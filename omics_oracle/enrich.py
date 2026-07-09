"""Pathway enrichment per dataset + cross-study shared pathways.

Uses Enrichr (MaayanLab) — no bulk download, just the gene list per study. For each study's DE
genes we get enriched pathways (term, p, adjusted-p, overlapping genes). Then we aggregate:

  - shared_pathways: enriched (padj<cutoff) in >= 2 studies.
  - surfaced: the feature you asked for — a pathway SIGNIFICANT in some studies but only
    SUB-THRESHOLD in others *despite many of its genes moving there*. So a process that is
    "almost significant in ours but massively enriched in theirs" is flagged, not hidden.

Graceful offline: if Enrichr is unreachable, enrichment is skipped with a recorded reason.
"""
from __future__ import annotations

import hashlib
import json
import os
from typing import Iterable, Optional

import requests

from .genes import map_to_symbols, detect_id_type
from .types import DEResult

_ADDLIST = "https://maayanlab.cloud/Enrichr/addList"
_ENRICH = "https://maayanlab.cloud/Enrichr/enrich"
_CACHE_DIR = os.environ.get("OMICS_ORACLE_ENRICHCACHE", "data/enrichcache")

_DEFAULT_LIB = {"human": "KEGG_2021_Human", "mouse": "KEGG_2019_Mouse"}


def default_library(species: str) -> str:
    return _DEFAULT_LIB.get(species, "KEGG_2021_Human")


def run_enrichr(genes: list[str], library: str, timeout: float = 30.0) -> Optional[dict]:
    """Return {term: {p, padj, combined, genes:[...]}} for a gene list, or None if unavailable."""
    genes = sorted({str(g).upper() for g in genes if str(g).strip()})
    if len(genes) < 3:
        return {}
    cached = _cache_get(genes, library)
    if cached is not None:
        return cached
    try:
        r = requests.post(_ADDLIST, files={"list": (None, "\n".join(genes)),
                                            "description": (None, "omics-oracle")}, timeout=timeout)
        r.raise_for_status()
        uid = r.json()["userListId"]
        e = requests.get(_ENRICH, params={"userListId": uid, "backgroundType": library},
                         timeout=timeout)
        e.raise_for_status()
        rows = e.json().get(library, []) or []
    except Exception:
        return None
    out = {}
    for row in rows:
        # Enrichr row: [rank, term, p, zscore, combined, [overlap_genes], adjp, ...]
        try:
            term, p, combined, overlap, adjp = row[1], row[2], row[4], row[5], row[6]
        except (IndexError, TypeError):
            continue
        out[term] = {"p": float(p), "padj": float(adjp), "combined": float(combined),
                     "genes": [str(g).upper() for g in overlap]}
    _cache_put(genes, library, out)
    return out


def _cache_key(genes: list[str], library: str) -> str:
    h = hashlib.sha1(("|".join(genes) + "::" + library).encode()).hexdigest()
    return os.path.join(_CACHE_DIR, f"{library}_{h}.json")


def _cache_get(genes: list[str], library: str) -> Optional[dict]:
    p = _cache_key(genes, library)
    if os.path.exists(p):
        try:
            with open(p) as fh:
                return json.load(fh)
        except Exception:
            return None
    return None


def _cache_put(genes: list[str], library: str, data: dict) -> None:
    os.makedirs(_CACHE_DIR, exist_ok=True)
    p = _cache_key(genes, library)
    try:
        with open(p, "w") as fh:
            json.dump(data, fh)
    except Exception:
        pass


def enrich_across_datasets(de_results: Iterable[DEResult],
                           user_genes: Optional[list[str]] = None,
                           species: str = "human", library: Optional[str] = None,
                           padj_cutoff: float = 0.05, min_overlap: int = 3,
                           de_padj: float = 0.05, de_lfc: float = 1.0,
                           top: int = 60) -> Optional[dict]:
    library = library or default_library(species)
    sources: list[tuple[str, list[str], str, str]] = []   # (id, genes, tier, source_db)
    for d in de_results:
        if d.status != "ok" or d.table is None or not len(d.table):
            continue
        t = d.table
        sig = t[(t["padj"].fillna(1.0) < de_padj) & (t["log2fc"].abs() >= de_lfc)]
        ids = [str(g) for g in sig["gene"]]
        mapping, _ = map_to_symbols(ids, detect_id_type(ids), species=species)
        genes = sorted({mapping.get(g, g.upper()) for g in ids})
        if genes:
            sources.append((d.accession, genes, d.tier, d.source_db))
    if user_genes:
        umap, _ = map_to_symbols([str(g) for g in user_genes], None, species=species)
        sources.append(("USER", sorted(set(umap.values())), "user", "user"))
    if len(sources) < 1:
        return None

    # run enrichment per source
    per_source: dict[str, dict] = {}
    unavailable = 0
    for sid, genes, tier, sdb in sources:
        res = run_enrichr(genes, library)
        if res is None:
            unavailable += 1
            continue
        per_source[sid] = {"terms": res, "tier": tier, "source_db": sdb, "n_genes": len(genes)}
    if not per_source:
        return {"library": library, "note": "enrichment unavailable (Enrichr unreachable)",
                "sources": [], "shared_pathways": [], "surfaced": []}

    # term index across sources
    term_index: dict[str, dict[str, dict]] = {}
    for sid, info in per_source.items():
        for term, row in info["terms"].items():
            term_index.setdefault(term, {})[sid] = row

    shared, surfaced = [], []
    for term, by_src in term_index.items():
        sig_in = [s for s, r in by_src.items() if r["padj"] < padj_cutoff]
        if len(sig_in) >= 2:
            shared.append({
                "term": term, "n_significant": len(sig_in), "sources_significant": sig_in,
                "by_source": {s: {"padj": round(r["padj"], 4), "overlap_n": len(r["genes"]),
                                  "genes": r["genes"]} for s, r in by_src.items()},
            })
        # surfacing: significant somewhere, but present-yet-subthreshold elsewhere with real coverage
        if sig_in:
            sub = [{"source": s, "padj": round(r["padj"], 4), "overlap_n": len(r["genes"]),
                    "genes": r["genes"], "tier": per_source[s]["tier"]}
                   for s, r in by_src.items()
                   if r["padj"] >= padj_cutoff and len(r["genes"]) >= min_overlap]
            if sub:
                surfaced.append({
                    "term": term, "significant_in": sig_in,
                    "subthreshold_in": sub,
                    "note": (f"enriched in {len(sig_in)} study(ies) but only sub-threshold in "
                             f"{len(sub)} other(s) despite {max(x['overlap_n'] for x in sub)} "
                             "of its genes moving there — worth a look"),
                })

    shared.sort(key=lambda x: x["n_significant"], reverse=True)
    surfaced.sort(key=lambda x: max((s["overlap_n"] for s in x["subthreshold_in"]), default=0),
                  reverse=True)
    out = {
        "library": library,
        "sources": [{"id": s, "tier": per_source[s]["tier"], "source_db": per_source[s]["source_db"],
                     "n_genes": per_source[s]["n_genes"],
                     "n_terms_sig": sum(1 for r in per_source[s]["terms"].values()
                                        if r["padj"] < padj_cutoff)} for s in per_source],
        "shared_pathways": shared[:top],
        "surfaced": surfaced[:top],
    }
    if unavailable:
        out["note"] = f"{unavailable} source(s) could not be enriched (Enrichr unreachable)"
    return out
