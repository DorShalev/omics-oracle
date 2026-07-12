"""ARCHS4 (Maayan) — a CO-EQUAL counts repository of uniformly-reprocessed RNA-seq.

ARCHS4 reprocesses GEO/SRA RNA-seq into one consistent raw-count matrix, so it is a peer of the
GEO per-study tables — not a fallback. For any accession the pipeline queries both repositories
and prefers raw counts regardless of which one supplied them; ARCHS4's uniform raw counts are
especially useful where a GEO study deposited only a normalized table (or nothing usable).

ARCHS4 distributes its matrix as a bulk H5 (>30 GB), which is not bundled. Point the source at a
local file with ``ARCHS4_H5`` or ``h5_path``. Without it the repository is INACTIVE:
``supports_counts()`` is False, the pipeline says so plainly, and it never invents counts. (The
former web API that allowed H5-free per-series access has been retired, so the H5 is required.)
"""
from __future__ import annotations

import os
from typing import Optional

from ..types import DatasetMatch, ExperimentSpec, CountMatrix


class ARCHS4Source:
    name = "ARCHS4"

    def __init__(self, h5_path: Optional[str] = None):
        self.h5_path = h5_path or os.environ.get("ARCHS4_H5")

    def _usable(self) -> bool:
        return bool(self.h5_path) and os.path.exists(self.h5_path)

    # -- peer finder over ARCHS4's own metadata (active only with the H5) ------------------

    def search(self, spec: ExperimentSpec, limit: int = 10) -> list[DatasetMatch]:
        """Find ARCHS4 series matching the spec. Active only with the H5; GEO E-utilities remains
        the richer index (ARCHS4 reprocesses the same GEO accessions), so this mainly confirms a
        study also has ARCHS4 uniform counts available."""
        if not self._usable():
            return []
        try:
            import archs4py as a4
            terms = [t for t in [spec.tissue, spec.contrast_case] if t] or (spec.keywords[:1] or [])
            if not terms:
                return []
            meta = a4.meta.meta(self.h5_path, terms[0], meta_fields=["series_id", "source_name_ch1"])
        except Exception:
            return []
        series: list[str] = []
        series_ids = meta.get("series_id", []) if hasattr(meta, "get") else []
        for sid in series_ids:
            if sid and sid not in series:
                series.append(sid)
        return [DatasetMatch(accession=s, source=self.name, organism=spec.organism,
                             omics_type="rna-seq") for s in series[:limit]]

    # -- counts ---------------------------------------------------------------------------

    def supports_counts(self) -> bool:
        return self._usable()

    def fetch_counts(self, match: DatasetMatch) -> Optional[CountMatrix]:
        if not self._usable():
            return None
        try:
            import archs4py as a4
            df = a4.data.series(self.h5_path, match.accession)
        except Exception:
            return None
        if df is None or getattr(df, "empty", True):
            return None
        df = _orient_genes_by_samples(df)
        return CountMatrix(
            accession=match.accession,
            genes=list(df.index),
            samples=list(df.columns),
            counts=df,
            sample_meta=None,
            source_db="ARCHS4", tier="clinical",
            provenance=f"ARCHS4 local H5 ({os.path.basename(self.h5_path)}), series {match.accession}",
        )


def _orient_genes_by_samples(df):
    """archs4py returns genes as rows; normalize defensively to genes x samples.

    Heuristic: gene symbols are strings; GSM sample ids start with 'GSM'. If the columns
    look like gene symbols and the index looks like samples, transpose.
    """
    cols = [str(c) for c in df.columns[:5]]
    idx = [str(i) for i in df.index[:5]]
    cols_are_samples = any(c.upper().startswith("GSM") for c in cols)
    idx_are_samples = any(i.upper().startswith("GSM") for i in idx)
    if idx_are_samples and not cols_are_samples:
        return df.T
    return df
