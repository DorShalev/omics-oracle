"""Fetcher: pull reprocessed counts for a matched dataset from a counts-capable source."""
from __future__ import annotations

from typing import Iterable, Optional

from .types import CountMatrix, DatasetMatch
from .sources.base import DatasetSource


def fetch_counts(match: DatasetMatch,
                 sources: Iterable[DatasetSource]) -> Optional[CountMatrix]:
    """Try each counts-capable source for this match; return the first CountMatrix."""
    for src in sources:
        try:
            if src.supports_counts():
                cm = src.fetch_counts(match)
                if cm is not None:
                    return cm
        except Exception:
            continue
    return None
