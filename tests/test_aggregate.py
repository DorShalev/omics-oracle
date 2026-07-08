import pandas as pd

from omics_oracle.types import DEResult
from omics_oracle.aggregate import aggregate


def _de(acc, rows):
    t = pd.DataFrame(rows, columns=["gene", "log2fc", "pvalue", "padj"])
    return DEResult(accession=acc, contrast="c", n_case=3, n_control=3, table=t, status="ok")


def test_concordant_gene_detected():
    d1 = _de("GSE1", [("A", 2.0, 1e-5, 1e-4), ("B", -1.5, 1e-4, 1e-3), ("C", 0.1, 0.5, 0.9)])
    d2 = _de("GSE2", [("A", 1.8, 1e-4, 1e-3), ("B", -1.2, 1e-3, 1e-2), ("D", 3.0, 1e-6, 1e-5)])
    agg = aggregate([d1, d2], min_datasets=2)
    genes = {g.gene: g for g in agg.genes}
    # A and B appear in both datasets with consistent direction
    assert "A" in genes and genes["A"].direction_agreement == 1.0 and genes["A"].n_up == 2
    assert "B" in genes and genes["B"].n_down == 2
    # C was not significant, D only in one dataset -> excluded from >=2 concordance
    assert "C" not in genes and "D" not in genes


def test_needs_two_datasets():
    d1 = _de("GSE1", [("A", 2.0, 1e-5, 1e-4)])
    agg = aggregate([d1], min_datasets=2)
    assert agg.genes == []
    assert "need >= 2" in agg.note


def test_discordant_gene_low_agreement():
    d1 = _de("GSE1", [("A", 2.0, 1e-5, 1e-4)])
    d2 = _de("GSE2", [("A", -2.0, 1e-5, 1e-4)])
    agg = aggregate([d1, d2], min_datasets=2)
    a = {g.gene: g for g in agg.genes}["A"]
    assert a.direction_agreement == 0.5  # one up, one down
