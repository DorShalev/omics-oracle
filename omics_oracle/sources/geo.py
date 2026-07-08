"""NCBI GEO finder via E-utilities (db=gds).

This is the default, key-free, remotely-queryable dataset finder. It returns rich metadata
(accession, taxon, dataset type, platform, sample count, summary) used for ranking and
mismatch flagging. It does not supply counts (see ARCHS4Source for that).
"""
from __future__ import annotations

import os
import time
from typing import Optional

import requests

from ..types import DatasetMatch, ExperimentSpec, CountMatrix

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Map our omics vocabulary to GEO "DataSet Type" filter values.
_GDS_TYPE = {
    "rna-seq": "expression profiling by high throughput sequencing",
    "scrna-seq": "expression profiling by high throughput sequencing",
    "atac-seq": "genome binding/occupancy profiling by high throughput sequencing",
    "proteomics": None,   # GEO is not the right index for proteomics
    "metabolomics": None,
    "unknown": None,
}

# Map our canonical organisms to GEO [Organism] tags.
_ORG_TAG = {
    "Homo sapiens": "Homo sapiens",
    "Mus musculus": "Mus musculus",
    "Rattus norvegicus": "Rattus norvegicus",
    "Saccharomyces cerevisiae": "Saccharomyces cerevisiae",
    "Danio rerio": "Danio rerio",
    "Drosophila melanogaster": "Drosophila melanogaster",
}


class GEOSource:
    name = "GEO"

    def __init__(self, api_key: Optional[str] = None, email: Optional[str] = None,
                 timeout: float = 20.0):
        # Optional NCBI credentials raise the rate limit; unused if absent.
        self.api_key = api_key or os.environ.get("NCBI_API_KEY")
        self.email = email or os.environ.get("NCBI_EMAIL")
        self.timeout = timeout

    # -- finder ---------------------------------------------------------------------------

    def search(self, spec: ExperimentSpec, limit: int = 10) -> list[DatasetMatch]:
        # Try progressively relaxed queries; Entrez ANDs bare words, so an over-stuffed
        # query returns nothing. Stop at the first term that yields hits.
        for term in self._candidate_terms(spec):
            ids = self._esearch(term, retmax=limit)
            if ids:
                return self._esummary(ids)
        return []

    # Generic tokens that over-constrain a GEO free-text AND without adding specificity.
    _FILLER = {"tumor", "tumour", "tissue", "adjacent", "normal", "cell", "cells",
               "profiling", "expression", "sequencing"}

    def _free_tokens(self, spec: ExperimentSpec, k: int) -> list[str]:
        bits: list[str] = []
        for b in [spec.tissue, spec.contrast_case, spec.contrast_control]:
            if b:
                bits += b.split()
        bits += spec.keywords
        toks: list[str] = []
        for t in bits:
            t = t.lower().strip(",.;:")
            if len(t) > 2 and t not in toks:
                toks.append(t)
        # prefer specific tokens over filler, but keep order otherwise
        specific = [t for t in toks if t not in self._FILLER]
        chosen = (specific or toks)[:k]
        return chosen

    def _candidate_terms(self, spec: ExperimentSpec):
        org = f'"{_ORG_TAG[spec.organism]}"[Organism]' if spec.organism in _ORG_TAG else None
        gds = _GDS_TYPE.get(spec.omics_type)
        gds_clause = f'"{gds}"[DataSet Type]' if gds else None

        def build(tokens, with_org, with_gds):
            parts = []
            if tokens:
                parts.append("(" + " ".join(f"{t}[All Fields]" for t in tokens) + ")")
            if with_org and org:
                parts.append(org)
            if with_gds and gds_clause:
                parts.append(gds_clause)
            return " AND ".join(parts) if parts else None

        seen = set()
        # most specific -> least: shrink free tokens, then drop dataset-type, then organism
        for k in (3, 2):
            for wg in (True, False):
                t = build(self._free_tokens(spec, k), with_org=True, with_gds=wg)
                if t and t not in seen:
                    seen.add(t); yield t
        t = build(self._free_tokens(spec, 2), with_org=False, with_gds=True)
        if t and t not in seen:
            seen.add(t); yield t
        # last resort: raw prompt words, organism only
        t = build(self._free_tokens(spec, 3), with_org=False, with_gds=False)
        if t and t not in seen:
            seen.add(t); yield t

    def _params(self, **kw) -> dict:
        p = {"retmode": "json", **kw}
        if self.api_key:
            p["api_key"] = self.api_key
        if self.email:
            p["email"] = self.email
        return p

    def _esearch(self, term: str, retmax: int) -> list[str]:
        r = requests.get(f"{_EUTILS}/esearch.fcgi",
                         params=self._params(db="gds", term=term, retmax=retmax),
                         timeout=self.timeout)
        r.raise_for_status()
        data = r.json().get("esearchresult", {})
        return data.get("idlist", []) or []

    def _esummary(self, ids: list[str]) -> list[DatasetMatch]:
        time.sleep(0.34 if not self.api_key else 0.0)  # respect default 3 req/s
        r = requests.get(f"{_EUTILS}/esummary.fcgi",
                         params=self._params(db="gds", id=",".join(ids)),
                         timeout=self.timeout)
        r.raise_for_status()
        result = r.json().get("result", {})
        out: list[DatasetMatch] = []
        for uid in result.get("uids", []):
            rec = result.get(uid, {})
            acc = rec.get("accession", "")
            out.append(DatasetMatch(
                accession=acc,
                source=self.name,
                title=rec.get("title", ""),
                organism=rec.get("taxon") or None,
                omics_type=_omics_from_gdstype(rec.get("gdstype", "")),
                platform=_platform_label(rec),
                n_samples=_as_int(rec.get("n_samples")),
                summary=rec.get("summary", ""),
                url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}" if acc else "",
            ))
        return out

    # -- counts (not supported here) ------------------------------------------------------

    def supports_counts(self) -> bool:
        return False

    def fetch_counts(self, match: DatasetMatch) -> Optional[CountMatrix]:
        return None


def _omics_from_gdstype(gdstype: str) -> str:
    g = (gdstype or "").lower()
    if "high throughput sequencing" in g and "expression" in g:
        return "rna-seq"
    if "binding/occupancy" in g:
        return "atac-seq"
    if "expression profiling by array" in g:
        return "microarray"
    return "unknown"


def _platform_label(rec: dict) -> Optional[str]:
    gpl = rec.get("gpl")
    return f"GPL{gpl}" if gpl else None


def _as_int(v) -> Optional[int]:
    try:
        return int(v)
    except (TypeError, ValueError):
        return None
