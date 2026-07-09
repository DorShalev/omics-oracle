"""Fetcher: pull reprocessed counts for a matched dataset from a counts-capable source."""
from __future__ import annotations

from typing import Iterable, Optional

from .types import CountMatrix, DatasetMatch
from .sources.base import DatasetSource


def fetch_counts(match: DatasetMatch,
                 sources: Iterable[DatasetSource]) -> Optional[CountMatrix]:
    """Try counts-capable sources for this match, PREFERRING raw counts.

    Sources are tried in order (e.g. GEO per-study first, then ARCHS4/clinical). A raw-count
    matrix wins immediately; if only normalized/unknown matrices are found, the first of those
    is returned so the study is still reported honestly (DE will be skipped downstream). This is
    what lets ARCHS4 (raw) rescue a study for which GEO deposited only a TPM table.
    """
    fallback: Optional[CountMatrix] = None
    for src in sources:
        try:
            if not src.supports_counts():
                continue
            cm = src.fetch_counts(match)
        except Exception:
            continue
        if cm is None:
            continue
        if getattr(cm, "data_type", "raw") == "raw":
            return cm
        if fallback is None:
            fallback = cm
    return fallback
