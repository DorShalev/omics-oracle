"""Offline pipeline test: stub finder + counts source, no network, no LLM."""
import numpy as np
import pandas as pd

from omics_oracle.types import CountMatrix, DatasetMatch, ExperimentSpec
from omics_oracle.pipeline import project_experiment
from omics_oracle import OracleReport, SCHEMA_VERSION


class _StubFinder:
    name = "GEO"
    def __init__(self, matches): self._m = matches
    def search(self, spec, limit=10): return self._m
    def supports_counts(self): return False
    def fetch_counts(self, m): return None


class _StubCounts:
    name = "ARCHS4"
    def __init__(self, by_acc): self._by = by_acc
    def search(self, spec, limit=10): return []
    def supports_counts(self): return True
    def fetch_counts(self, m): return self._by.get(m.accession)


def _counts(acc, seed):
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(150)]
    ctrl = [f"{acc}_c{i}" for i in range(4)]
    case = [f"{acc}_t{i}" for i in range(4)]
    base = rng.integers(50, 400, size=len(genes))
    mat = np.vstack([rng.poisson(base) for _ in range(8)]).T.astype(float)  # genes x 8
    mat[:8, 4:] *= 4.0  # first 8 genes up in "case" columns
    df = pd.DataFrame(mat.astype(int), index=genes, columns=ctrl + case)
    # metadata labels must match the parsed contrast ("tumor" vs "control")
    meta = pd.DataFrame({"condition": (["control"] * 4) + (["tumor"] * 4)}, index=ctrl + case)
    return CountMatrix(accession=acc, genes=genes, samples=ctrl + case, counts=df,
                       sample_meta=meta, provenance="stub")


def test_pipeline_end_to_end_offline():
    matches = [
        DatasetMatch(accession="GSE_A", source="GEO", title="human tumor vs control rna-seq",
                     summary="tumor control", organism="Homo sapiens", omics_type="rna-seq",
                     n_samples=8),
        DatasetMatch(accession="GSE_B", source="GEO", title="human tumor vs control rna-seq",
                     summary="tumor control", organism="Homo sapiens", omics_type="rna-seq",
                     n_samples=8),
    ]
    counts = {"GSE_A": _counts("GSE_A", 1), "GSE_B": _counts("GSE_B", 2)}
    report = project_experiment(
        "RNA-seq of human tumor vs control, 4 per group",
        finder_sources=[_StubFinder(matches)],
        counts_sources=[_StubCounts(counts)],
    )
    assert isinstance(report, OracleReport)
    d = report.to_dict()
    assert d["schema_version"] == SCHEMA_VERSION
    assert len(d["matches"]) == 2
    ok = [r for r in report.de_results if r.status == "ok"]
    assert len(ok) == 2
    assert report.concordance is not None and report.concordance.n_datasets == 2
    # markdown renders without error and cites an accession
    md = report.to_markdown()
    assert "GSE_A" in md and "Projection" in md


def test_metadata_only_when_no_counts_source():
    matches = [DatasetMatch(accession="GSE_A", source="GEO", title="human rna-seq",
                            organism="Homo sapiens", omics_type="rna-seq", n_samples=8)]
    report = project_experiment(
        "RNA-seq of human tumor vs control",
        finder_sources=[_StubFinder(matches)],
        counts_sources=[],  # no counts backend
    )
    assert all(r.status != "ok" for r in report.de_results)
    assert report.concordance is None
    # with no counts source, each dataset is honestly skipped (not faked)
    assert all(r.status == "skipped" for r in report.de_results)
    assert report.projection.feasibility in ("cautionary", "thin-evidence")
