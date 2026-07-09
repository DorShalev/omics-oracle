"""Multi-source counts routing: prefer raw, let ARCHS4/clinical rescue normalized GEO."""
import numpy as np
import pandas as pd

from omics_oracle.types import CountMatrix, DatasetMatch
from omics_oracle.fetch import fetch_counts


def _cm(acc, data_type, source_db, tier):
    genes = [f"G{i}" for i in range(20)]
    cols = [f"{acc}_{i}" for i in range(4)]
    df = pd.DataFrame(np.ones((20, 4)), index=genes, columns=cols)
    return CountMatrix(accession=acc, genes=genes, samples=cols, counts=df,
                       data_type=data_type, source_db=source_db, tier=tier)


class _Src:
    def __init__(self, name, cm): self.name = name; self._cm = cm
    def supports_counts(self): return True
    def fetch_counts(self, m): return self._cm
    def search(self, spec, limit=10): return []


def test_raw_source_wins_over_normalized():
    geo = _Src("GEO", _cm("GSE1", "normalized", "GEO", "primary"))
    archs4 = _Src("ARCHS4", _cm("GSE1", "raw", "ARCHS4", "clinical"))
    cm = fetch_counts(DatasetMatch(accession="GSE1", source="GEO"), [geo, archs4])
    assert cm.data_type == "raw" and cm.tier == "clinical"  # ARCHS4 rescued it


def test_normalized_returned_when_no_raw_available():
    geo = _Src("GEO", _cm("GSE1", "normalized", "GEO", "primary"))
    cm = fetch_counts(DatasetMatch(accession="GSE1", source="GEO"), [geo])
    assert cm.data_type == "normalized" and cm.tier == "primary"


def test_first_raw_short_circuits():
    geo = _Src("GEO", _cm("GSE1", "raw", "GEO", "primary"))
    # a second source that would raise proves we stop at the first raw hit
    class Boom:
        name = "X"
        def supports_counts(self): raise AssertionError("should not be reached")
        def fetch_counts(self, m): raise AssertionError
        def search(self, s, limit=10): return []
    cm = fetch_counts(DatasetMatch(accession="GSE1", source="GEO"), [geo, Boom()])
    assert cm.tier == "primary" and cm.source_db == "GEO"
