"""Gene × dataset log2FC matrix — the data behind the interactive concordance heatmap.

Y-axis = the user's genes (or, if none supplied, the union of DE genes ranked by breadth).
X-axis = matched datasets, tagged with their match_tier so the UI can group columns.
Cells  = that gene's log2FC in each dataset's FULL DE table (not only where significant), so
concordance is visible even where a gene is sub-threshold — red-in-yours + red-in-theirs = a
consistent signal. The UI does the Top-10/20/50/100/All and Up/Down/All toggles by filtering
the returned `genes` (each carries rank + direction); the backend returns everything.

All ids are harmonized to HGNC symbols so datasets keyed by Ensembl vs symbol align.
"""
from __future__ import annotations

from typing import Iterable, Optional

import numpy as np

from .genes import map_to_symbols, map_symbols_to, detect_id_type, strip_version
from .types import DatasetMatch, DEResult


def _symbol_fc_map(table, species: str) -> dict[str, float]:
    """DE table -> {symbol: log2fc}, collapsing dup symbols by max |log2fc|.

    Forward-maps every id, so only call it on a SMALL table (e.g. the significant subset used to
    pick Y genes) — never a study's full 60k-row table. For full tables use ``_dataset_fc_for``,
    which maps only the handful of Y genes back into the table's id space instead.
    """
    ids = [str(g) for g in table["gene"]]
    mapping, _ = map_to_symbols(ids, detect_id_type(ids), species=species)
    out: dict[str, float] = {}
    for gid, fc in zip(ids, table["log2fc"]):
        try:
            fc = float(fc)
        except (TypeError, ValueError):
            continue
        sym = mapping.get(gid, gid.upper())
        if sym not in out or abs(fc) > abs(out[sym]):
            out[sym] = fc
    return out


def _dataset_fc_for(table, y_syms: list[str], species: str) -> dict[str, Optional[float]]:
    """{sym: log2fc} for just the Y genes, mapping only what's needed (bounded by |y_syms|).

    Builds one native-id -> max|log2fc| index over the study's own table (no mygene), then, for
    id-keyed studies, reverse-maps the Y symbols into that id space with a single small mygene
    call. This keeps heatmap cost proportional to the Y-gene count, not the 60k-row table size.
    """
    ids = [str(g) for g in table["gene"]]
    id_type = detect_id_type(ids)

    fc_by_key: dict[str, float] = {}
    if id_type == "symbol":
        for gid, fc in zip(ids, table["log2fc"]):
            try:
                fc = float(fc)
            except (TypeError, ValueError):
                continue
            k = gid.upper()
            if k not in fc_by_key or abs(fc) > abs(fc_by_key[k]):
                fc_by_key[k] = fc
        return {s: fc_by_key.get(s.upper()) for s in y_syms}

    for gid, fc in zip(ids, table["log2fc"]):
        try:
            fc = float(fc)
        except (TypeError, ValueError):
            continue
        k = strip_version(gid).upper()
        if k not in fc_by_key or abs(fc) > abs(fc_by_key[k]):
            fc_by_key[k] = fc
    sym2id = map_symbols_to(y_syms, id_type, species=species)
    out: dict[str, Optional[float]] = {}
    for s in y_syms:
        gid = sym2id.get(s.upper())
        out[s] = fc_by_key.get(strip_version(gid).upper()) if gid else None
    return out


def build_heatmap(de_results: Iterable[DEResult],
                  matches: Optional[list[DatasetMatch]] = None,
                  user_genes: Optional[list[str]] = None,
                  species: str = "human",
                  de_padj: float = 0.05, de_lfc: float = 1.0) -> Optional[dict]:
    ok = [d for d in de_results if d.status == "ok" and d.table is not None and len(d.table)]
    if not ok:
        return None

    tier_of = {m.accession: (m.match_tier or "") for m in (matches or [])}

    # choose Y genes FIRST (cheap): the user's list, else the union of significant genes.
    # Note the union maps only each study's SIGNIFICANT subset (small), not the full table.
    if user_genes:
        umap, _ = map_to_symbols([str(g) for g in user_genes], None, species=species)
        y_syms = list(dict.fromkeys(umap.values()))
    else:
        breadth: dict[str, int] = {}
        for d in ok:
            t = d.table
            sig = t[(t["padj"].fillna(1.0) < de_padj) & (t["log2fc"].abs() >= de_lfc)]
            for sym in _symbol_fc_map(sig, species):
                breadth[sym] = breadth.get(sym, 0) + 1
        y_syms = sorted(breadth, key=lambda s: breadth[s], reverse=True)

    # cells: for each dataset, look up ONLY the Y genes (bounded reverse-map), not the whole table
    ds_maps = {d.accession: _dataset_fc_for(d.table, y_syms, species) for d in ok}

    accs = [d.accession for d in ok]
    matrix, genes_meta = [], []
    for sym in y_syms:
        row = [ds_maps[a].get(sym) for a in accs]
        present = [v for v in row if v is not None]
        n_de = sum(1 for v in present if abs(v) >= de_lfc)
        mean_fc = round(float(np.mean(present)), 3) if present else None
        up = sum(1 for v in present if v > 0)
        down = sum(1 for v in present if v < 0)
        matrix.append([None if v is None else round(float(v), 3) for v in row])
        genes_meta.append({"gene": sym, "n_datasets_present": len(present),
                           "n_datasets_de": n_de, "mean_log2fc": mean_fc,
                           "up": up, "down": down,
                           "direction": "up" if (mean_fc or 0) > 0 else "down"})

    # rank Y genes for the Top-N toggle: broadest DE first, then strongest effect
    order = sorted(range(len(genes_meta)),
                   key=lambda i: (genes_meta[i]["n_datasets_de"],
                                  abs(genes_meta[i]["mean_log2fc"] or 0)), reverse=True)
    genes_meta = [dict(genes_meta[i], rank=r + 1) for r, i in enumerate(order)]
    matrix = [matrix[i] for i in order]

    datasets = [{"accession": a, "match_tier": tier_of.get(a, ""),
                 "source_db": d.source_db, "tier": d.tier}
                for a, d in zip(accs, ok)]
    # order columns by tier (direct, proxy, pathway, none, unknown)
    from .tiers import TIER_ORDER
    col_order = sorted(range(len(datasets)),
                       key=lambda j: TIER_ORDER.get(datasets[j]["match_tier"], 9))
    datasets = [datasets[j] for j in col_order]
    matrix = [[row[j] for j in col_order] for row in matrix]

    return {"genes": genes_meta, "datasets": datasets, "matrix": matrix,
            "value": "log2FC", "scope": "user_genes" if user_genes else "union_of_DE"}
