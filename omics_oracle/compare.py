"""Cross-dataset comparison — the structured data behind the Venn + gene-search UI.

Produces JSON-ready output:
  - datasets:  one row per reanalyzed study (accession, source_db, tier, n_de)
  - genes:     per-gene detection across studies — "in how many datasets is this gene DE"
               (drives the gene search bar); optionally scoped to the USER's own gene list
  - pairwise:  DE-gene overlap between every pair of studies (drives the interactive Venn)
  - user:      if a user gene list is supplied, its overlap with each study

All gene ids are harmonized to HGNC symbols first, so Ensembl-vs-symbol studies overlap.
"""
from __future__ import annotations

from typing import Iterable, Optional

from .genes import map_to_symbols, detect_id_type
from .types import DEResult


def build_comparison(de_results: Iterable[DEResult],
                     user_genes: Optional[list[str]] = None,
                     species: str = "human",
                     padj: float = 0.05, lfc: float = 1.0) -> Optional[dict]:
    ok = [d for d in de_results if d.status == "ok" and d.table is not None and len(d.table)]
    if not ok:
        return None

    # per-study DE gene set in symbol space: accession -> {symbol: log2fc}
    ds: dict[str, dict] = {}
    for d in ok:
        t = d.table
        sig = t[(t["padj"].fillna(1.0) < padj) & (t["log2fc"].abs() >= lfc)]
        ids = [str(g) for g in sig["gene"]]
        mapping, _ = map_to_symbols(ids, detect_id_type(ids), species=species)
        sym_fc: dict[str, float] = {}
        for g, fc in zip(ids, sig["log2fc"]):
            sym_fc[mapping.get(g, g.upper())] = float(fc)
        ds[d.accession] = {"map": sym_fc, "source_db": d.source_db, "tier": d.tier}

    accs = list(ds)
    sets = {a: set(ds[a]["map"]) for a in accs}

    # per-gene detection across studies
    gene_hits: dict[str, list[tuple[str, float]]] = {}
    for a in accs:
        for sym, fc in ds[a]["map"].items():
            gene_hits.setdefault(sym, []).append((a, fc))

    if user_genes:
        umap, _ = map_to_symbols([str(g) for g in user_genes], None, species=species)
        user_syms = set(umap.values())
        target = sorted(user_syms)
    else:
        user_syms = None
        target = sorted(gene_hits)

    genes_out = []
    for sym in target:
        hits = gene_hits.get(sym, [])
        up = sum(1 for _, fc in hits if fc > 0)
        down = sum(1 for _, fc in hits if fc < 0)
        mean_fc = round(sum(fc for _, fc in hits) / len(hits), 3) if hits else None
        genes_out.append({"gene": sym, "n_datasets": len(hits),
                          "datasets": [a for a, _ in hits], "up": up, "down": down,
                          "mean_log2fc": mean_fc})
    genes_out.sort(key=lambda g: (g["n_datasets"], abs(g["mean_log2fc"] or 0)), reverse=True)

    # pairwise overlap (Venn data)
    pairwise = []
    for i in range(len(accs)):
        for j in range(i + 1, len(accs)):
            a, b = accs[i], accs[j]
            sa, sb = sets[a], sets[b]
            inter, union = len(sa & sb), len(sa | sb)
            pairwise.append({"a": a, "b": b, "shared": inter,
                             "a_only": len(sa - sb), "b_only": len(sb - sa),
                             "jaccard": round(inter / union, 3) if union else 0.0})

    result = {
        "datasets": [{"accession": a, "source_db": ds[a]["source_db"], "tier": ds[a]["tier"],
                      "n_de": len(ds[a]["map"])} for a in accs],
        "genes": genes_out,
        "pairwise": pairwise,
        "thresholds": {"padj": padj, "abs_log2fc": lfc},
    }
    if user_syms is not None:
        result["user"] = {
            "n_genes": len(user_syms),
            "vs_datasets": [
                {"accession": a, "shared": len(user_syms & sets[a]),
                 "shared_genes": sorted(user_syms & sets[a])[:500],
                 "jaccard": round(len(user_syms & sets[a]) / len(user_syms | sets[a]), 3)
                 if (user_syms | sets[a]) else 0.0}
                for a in accs],
        }
    return result
