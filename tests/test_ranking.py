from omics_oracle.types import ExperimentSpec, DatasetMatch
from omics_oracle.finder import score_match, find_datasets


def _spec():
    return ExperimentSpec(
        raw_prompt="x", organism="Homo sapiens", tissue="colorectal",
        contrast_case="tumor", contrast_control="normal", omics_type="rna-seq",
        n_per_group=10,
    )


def test_perfect_match_scores_high_no_blocking():
    m = DatasetMatch(accession="GSE1", source="GEO",
                     title="RNA-seq of human colorectal tumor vs normal",
                     summary="colorectal tumor and adjacent normal tissue",
                     organism="Homo sapiens", omics_type="rna-seq", n_samples=40)
    score_match(_spec(), m)
    assert not m.blocking
    assert m.similarity >= 0.9


def test_wrong_organism_is_blocking():
    m = DatasetMatch(accession="GSE2", source="GEO", title="mouse colorectal tumor rna-seq",
                     summary="colorectal tumor vs normal", organism="Mus musculus",
                     omics_type="rna-seq", n_samples=40)
    score_match(_spec(), m)
    assert m.blocking
    assert any(x.field == "organism" and x.severity == "blocking" for x in m.mismatches)
    # the caveat should name the found organism
    assert any("Mus musculus" in x.note for x in m.mismatches)


def test_wrong_cellline_flags_tissue_mismatch():
    m = DatasetMatch(accession="GSE3", source="GEO", title="HEK293 kidney cells rna-seq",
                     summary="HEK293 epithelial kidney cell line experiment",
                     organism="Homo sapiens", omics_type="rna-seq", n_samples=6)
    score_match(_spec(), m)
    assert any(x.field == "tissue" for x in m.mismatches)
    assert m.similarity < 0.9


def test_ranking_puts_nonblocking_first():
    good = DatasetMatch(accession="GSE_good", source="GEO",
                        title="human colorectal tumor vs normal rna-seq",
                        summary="colorectal tumor normal", organism="Homo sapiens",
                        omics_type="rna-seq", n_samples=40)
    bad = DatasetMatch(accession="GSE_bad", source="GEO", title="mouse colorectal tumor rna-seq",
                       summary="colorectal", organism="Mus musculus", omics_type="rna-seq",
                       n_samples=40)

    class Src:
        name = "GEO"
        def search(self, spec, limit=10): return [bad, good]
        def supports_counts(self): return False
        def fetch_counts(self, m): return None

    ranked = find_datasets(_spec(), sources=[Src()], limit=5)
    assert ranked[0].accession == "GSE_good"
    assert ranked[-1].accession == "GSE_bad"
