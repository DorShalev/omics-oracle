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


def map_symbols_to(symbols: Iterable[str], target_type: str,
                   species: str = "human") -> dict[str, str]:
    """Reverse of :func:`map_to_symbols`: HGNC symbols -> ids of ``target_type``.

    Returns ``{SYMBOL_UPPER: id}`` for symbols that resolve (misses omitted). Bounded by the
    number of symbols asked for — used by the heatmap to look up only the handful of Y genes in
    each dataset's id space, instead of forward-mapping the study's whole 60k-row table. Cached
    per (species, target_type). Returns ``{}`` when mapping is unavailable (offline)."""
    syms = [str(s).strip() for s in symbols if str(s).strip()]
    if not syms:
        return {}
    if target_type == "symbol":
        return {s.upper(): s.upper() for s in syms}
    if target_type not in ("ensembl", "entrez"):
        return {}
    cache = _load_rev_cache(species, target_type)
    want = sorted({s.upper() for s in syms})
    unknown = [s for s in want if s not in cache]
    if unknown:
        fetched = _query_mygene_symbols(unknown, target_type, species)
        if fetched is None:                                 # offline / mygene missing
            return {s: cache[s] for s in want if cache.get(s)}
        cache.update(fetched)
        _save_rev_cache(species, target_type, cache)
    return {s: cache[s] for s in want if cache.get(s)}


def _query_mygene_symbols(symbols: list[str], target_type: str,
                          species: str) -> Optional[dict[str, str]]:
    try:
        import mygene
    except Exception:
        return None
    field = {"ensembl": "ensembl.gene", "entrez": "entrezgene"}[target_type]
    try:
        mg = mygene.MyGeneInfo()
        res = mg.querymany(symbols, scopes="symbol", fields=field, species=species,
                           as_dataframe=False, verbose=False)
    except Exception:
        return None
    out: dict[str, str] = {}
    for r in res:
        q = str(r.get("query", "")).upper()
        if not q:
            continue
        val = None
        if target_type == "entrez":
            v = r.get("entrezgene")
            val = str(v) if v is not None else None
        else:
            ens = r.get("ensembl")
            if isinstance(ens, dict):
                val = ens.get("gene")
            elif isinstance(ens, list) and ens:
                val = ens[0].get("gene")
        if val:
            out[q] = strip_version(str(val)).upper()
    for s in symbols:                                       # record misses so we don't re-query
        out.setdefault(s.upper(), "")
    return out


def _rev_cache_path(species: str, target_type: str) -> str:
    return os.path.join(_CACHE_DIR, f"{species}_to_{target_type}.json")


def _load_rev_cache(species: str, target_type: str) -> dict[str, str]:
    p = _rev_cache_path(species, target_type)
    if os.path.exists(p):
        try:
            with open(p) as fh:
                return json.load(fh)
        except Exception:
            return {}
    return {}


def _save_rev_cache(species: str, target_type: str, cache: dict[str, str]) -> None:
    os.makedirs(_CACHE_DIR, exist_ok=True)
    tmp = _rev_cache_path(species, target_type) + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(cache, fh)
    os.replace(tmp, _rev_cache_path(species, target_type))


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
