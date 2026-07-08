"""Dataset finder: query sources, then rank uniformly with an explainable score.

Similarity starts at 1.0 and is reduced by *attributed* penalties, so every point lost is
tied to a named Mismatch. A blocking mismatch (wrong organism / wrong omics) is recorded so
the pipeline can exclude that dataset from the counts/DE step even when its text overlaps.
"""
from __future__ import annotations

import re
from typing import Iterable, Optional

from .types import DatasetMatch, ExperimentSpec, Mismatch
from .sources.base import DatasetSource
from .sources.geo import GEOSource

# penalty, severity per field
_PENALTY = {
    "organism": (0.50, "blocking"),
    "omics_type": (0.40, "blocking"),
    "tissue": (0.25, "major"),
    "contrast": (0.20, "major"),
    "n": (0.05, "minor"),
}

_TOKEN = re.compile(r"[a-z0-9]+")


def find_datasets(spec: ExperimentSpec, sources: Optional[Iterable[DatasetSource]] = None,
                  limit: int = 10, per_source: int = 10) -> list[DatasetMatch]:
    sources = list(sources) if sources is not None else [GEOSource()]
    candidates: list[DatasetMatch] = []
    seen: set[str] = set()
    for src in sources:
        try:
            hits = src.search(spec, limit=per_source)
        except Exception:
            hits = []
        for h in hits:
            key = f"{h.source}:{h.accession}"
            if h.accession and key not in seen:
                seen.add(key)
                candidates.append(h)

    for m in candidates:
        score_match(spec, m)

    # rank: non-blocking first, then by similarity desc
    candidates.sort(key=lambda m: (m.blocking, -m.similarity))
    return candidates[:limit]


def score_match(spec: ExperimentSpec, match: DatasetMatch) -> DatasetMatch:
    """Populate match.similarity and match.mismatches in place; returns the match."""
    score = 1.0
    mismatches: list[Mismatch] = []
    haystack = f"{match.title} {match.summary}".lower()

    # organism (blocking)
    if spec.organism and match.organism:
        if not _same_organism(spec.organism, match.organism):
            pen, sev = _PENALTY["organism"]
            score -= pen
            mismatches.append(Mismatch("organism", spec.organism, match.organism, sev,
                                       f"{match.accession} is {match.organism}, not {spec.organism}"))

    # omics type
    if spec.omics_type != "unknown":
        if match.omics_type == "unknown":
            # assay could not be confirmed from source metadata — flag, but don't block
            score -= 0.10
            mismatches.append(Mismatch("omics_type", spec.omics_type, "unknown", "minor",
                                       f"assay type of {match.accession} could not be confirmed "
                                       f"from {match.source} metadata"))
        elif not _compatible_omics(spec.omics_type, match.omics_type):
            pen, sev = _PENALTY["omics_type"]
            score -= pen
            mismatches.append(Mismatch("omics_type", spec.omics_type, match.omics_type, sev,
                                       f"{match.accession} is {match.omics_type}, not {spec.omics_type}"))

    # tissue (major) — token overlap against title+summary
    if spec.tissue:
        if not _tokens_present(spec.tissue, haystack):
            pen, sev = _PENALTY["tissue"]
            score -= pen
            mismatches.append(Mismatch("tissue", spec.tissue, None, sev,
                                       f"requested tissue '{spec.tissue}' not evident in {match.accession} metadata"))

    # contrast (major) — at least one contrast term should appear
    contrast_terms = " ".join(t for t in [spec.contrast_case, spec.contrast_control] if t)
    if contrast_terms:
        if not _tokens_present(contrast_terms, haystack):
            pen, sev = _PENALTY["contrast"]
            score -= pen
            mismatches.append(Mismatch("contrast", contrast_terms, None, sev,
                                       f"contrast terms '{contrast_terms}' not evident in {match.accession} metadata"))

    # platform: informational only. The specific Illumina instrument (NovaSeq/HiSeq/NextSeq)
    # does not affect comparability — ARCHS4 reprocesses all of them uniformly — so we do NOT
    # penalize or filter on it. A genuine cross-technology gap (array vs seq) is already caught
    # by the omics_type check above. We record a note only when a non-Illumina platform is
    # explicitly named, as a heads-up, without touching the score.
    if spec.platform and "illumina" not in spec.platform.lower():
        mismatches.append(Mismatch("platform", spec.platform, match.platform, "minor",
                                   f"planned platform '{spec.platform}' is non-Illumina; most public "
                                   "RNA-seq (and ARCHS4) is Illumina — comparability is usually fine "
                                   "but library chemistry may differ"))

    # sample-size adequacy (minor)
    if spec.n_per_group and match.n_samples is not None:
        if match.n_samples < 2 * spec.n_per_group:
            pen, sev = _PENALTY["n"]
            score -= pen
            mismatches.append(Mismatch("n", f">= {2 * spec.n_per_group}", match.n_samples, sev,
                                       f"{match.accession} has {match.n_samples} samples; smaller than the planned design"))

    match.similarity = round(max(0.0, min(1.0, score)), 3)
    match.mismatches = mismatches
    return match


# ---- helpers ---------------------------------------------------------------------------

def _same_organism(a: str, b: str) -> bool:
    return a.strip().lower() == b.strip().lower()


def _compatible_omics(spec_omics: str, match_omics: str) -> bool:
    if spec_omics == match_omics:
        return True
    # bulk vs single-cell RNA are the same assay family for finding purposes
    rna = {"rna-seq", "scrna-seq"}
    return spec_omics in rna and match_omics in rna


def _tokens(text: str) -> set[str]:
    return set(_TOKEN.findall(text.lower()))


def _tokens_present(needle: str, haystack: str) -> bool:
    want = {t for t in _tokens(needle) if len(t) > 2}
    if not want:
        return True
    have = _tokens(haystack)
    # accept if any salient token (or a stem prefix) is present
    for w in want:
        if w in have:
            return True
        if any(h.startswith(w[:5]) for h in have if len(w) >= 5):
            return True
    return False
