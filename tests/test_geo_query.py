"""GEO query builder: the sequencing dataset-type filter must never be relaxed away."""
from omics_oracle.types import ExperimentSpec
from omics_oracle.sources.geo import GEOSource


def _spec(omics):
    return ExperimentSpec(raw_prompt="x", organism="Homo sapiens", tissue="colorectal",
                          contrast_case="tumor", contrast_control="normal", omics_type=omics,
                          keywords=["colorectal", "tumor", "normal"])


def test_rnaseq_keeps_datasettype_filter_in_every_candidate():
    terms = list(GEOSource()._candidate_terms(_spec("rna-seq")))
    assert terms, "expected candidate query terms"
    # every relaxation step must still constrain to sequencing — else microarray leaks in
    assert all("DataSet Type" in t for t in terms)


def test_unknown_omics_may_relax_freely():
    terms = list(GEOSource()._candidate_terms(_spec("unknown")))
    # no dataset-type mapping for unknown omics, so no such clause is required
    assert terms
    assert all("DataSet Type" not in t for t in terms)
