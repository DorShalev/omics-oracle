"""The DatasetSource protocol — the seam every public-DB backend implements.

A source may supply metadata (for finding/ranking), counts (for reanalysis), or both.
Keeping these capabilities on one interface lets the finder mix a metadata-rich source
(GEO) with a counts source (ARCHS4) without special-casing either.
"""
from __future__ import annotations

from typing import Optional, Protocol, runtime_checkable

from ..types import DatasetMatch, ExperimentSpec, CountMatrix


@runtime_checkable
class DatasetSource(Protocol):
    name: str

    def search(self, spec: ExperimentSpec, limit: int = 10) -> list[DatasetMatch]:
        """Return candidate datasets (unranked or source-ranked); finder re-ranks uniformly."""
        ...

    def supports_counts(self) -> bool:
        """Whether fetch_counts can return reprocessed counts for a match."""
        ...

    def fetch_counts(self, match: DatasetMatch) -> Optional[CountMatrix]:
        """Return reprocessed counts for a match, or None if unavailable."""
        ...
