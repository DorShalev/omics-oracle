"""Gene-identifier harmonization.

Studies deposit tables keyed by different identifiers — Ensembl (often with a version suffix),
HGNC symbols, or Entrez ids. Comparing datasets (concordance, overlap/Venn, shared pathways)
requires mapping them all to ONE space; otherwise Ensembl-vs-symbol studies show zero overlap
even when they agree. We map everything to uppercased HGNC symbols via mygene.info, cached on
disk. If mapping is unavailable (offline), we fall back to normalized raw ids and say so.
"""
from __future__ import annotations

import json
import os
import re
from typing import Iterable, Optional

_ENSEMBL = re.compile(r"^ENS[A-Z]*G\d{6,}(\.\d+)?$", re.IGNORECASE)
_ENTREZ = re.compile(r"^\d+$")

_CACHE_DIR = os.environ.get("OMICS_ORACLE_GENECACHE", "data/genecache")


def strip_version(gene_id: str) -> str:
    """ENSG00000141510.5 -> ENSG00000141510 (leave symbols/entrez untouched)."""
    s = str(gene_id).strip()
    if _ENSEMBL.match(s):
        return s.split(".")[0]
    return s


def detect_id_type(ids: Iterable[str]) -> str:
    """Return 'ensembl' | 'entrez' | 'symbol' | 'mixed' | 'unknown' for a set of ids."""
    ens = ent = sym = 0
    n = 0
    for g in ids:
        s = str(g).strip()
        if not s:
            continue
        n += 1
        if _ENSEMBL.match(s):
            ens += 1
        elif _ENTREZ.match(s):
            ent += 1
        else:
            sym += 1
        if n >= 200:
            break
    if n == 0:
        return "unknown"
    top = max(ens, ent, sym)
    if top / n < 0.7:
        return "mixed"
    return {ens: "ensembl", ent: "entrez", sym: "symbol"}[top]


def map_to_symbols(ids: list[str], id_type: Optional[str] = None,
                   species: str = "human") -> tuple[dict[str, str], bool]:
    """Map ids -> uppercased HGNC symbols. Returns (mapping, ok).

    ``ok`` is False if mapping could not run (offline / mygene missing); callers should then
    treat ids as-is and warn rather than silently comparing across identifier spaces.
    """
    ids = [str(i) for i in ids]
    id_type = id_type or detect_id_type(ids)
    if id_type == "symbol":
        return {i: i.upper() for i in ids}, True
    if id_type == "unknown":
        return {i: i.upper() for i in ids}, False

    stripped = {i: strip_version(i) for i in ids}
    cache = _load_cache(species)
    unknown = sorted({v for v in stripped.values() if v not in cache})
    if unknown:
        fetched = _query_mygene(unknown, id_type, species)
        if fetched is None:                        # mapping unavailable
            return {i: strip_version(i).upper() for i in ids}, False
        cache.update(fetched)
        _save_cache(species, cache)

    mapping = {}
    for original, key in stripped.items():
        sym = cache.get(key)
        mapping[original] = (sym or key).upper()   # fall back to the (version-stripped) id
    return mapping, True


# ---- mygene backend --------------------------------------------------------------------

def _query_mygene(ids: list[str], id_type: str, species: str) -> Optional[dict[str, str]]:
    try:
        import mygene
    except Exception:
        return None
    scopes = {"ensembl": "ensembl.gene", "entrez": "entrezgene"}.get(id_type, "ensembl.gene")
    try:
        mg = mygene.MyGeneInfo()
        res = mg.querymany(ids, scopes=scopes, fields="symbol", species=species,
                           as_dataframe=False, verbose=False)
    except Exception:
        return None
    out: dict[str, str] = {}
    for r in res:
        q = str(r.get("query", ""))
        sym = r.get("symbol")
        if q and sym:
            out[q] = sym
    # record misses too, so we don't re-query them every run
    for i in ids:
        out.setdefault(i, "")
    return out


def _cache_path(species: str) -> str:
    return os.path.join(_CACHE_DIR, f"{species}_symbols.json")


def _load_cache(species: str) -> dict[str, str]:
    p = _cache_path(species)
    if os.path.exists(p):
        try:
            with open(p) as fh:
                return json.load(fh)
        except Exception:
            return {}
    return {}


def _save_cache(species: str, cache: dict[str, str]) -> None:
    os.makedirs(_CACHE_DIR, exist_ok=True)
    tmp = _cache_path(species) + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(cache, fh)
    os.replace(tmp, _cache_path(species))
