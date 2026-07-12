"""Match-tiers, heatmap matrix, narrative validation."""
import pandas as pd

from omics_oracle.types import DatasetMatch, DEResult, Mismatch, ExperimentSpec
from omics_oracle.tiers import classify_tiers
from omics_oracle.heatmap import build_heatmap
from omics_oracle.validate import validate_narrative, collect_allowed_genes


def _spec():
    return ExperimentSpec(raw_prompt="x", organism="Homo sapiens", tissue="colorectal",
                          contrast_case="tumor", contrast_control="normal", omics_type="rna-seq")


def test_tier_direct_proxy_pathway_none():
    direct = DatasetMatch(accession="GSE1", source="GEO", judge={"verdict": "relevant"})
    proxy = DatasetMatch(accession="GSE2", source="GEO", judge={"verdict": "partial"},
                         mismatches=[Mismatch("tissue", "colorectal", None, "major")])
    pathway = DatasetMatch(accession="GSE3", source="GEO", judge={"verdict": "irrelevant"},
                           mismatches=[Mismatch("tissue", "colorectal", None, "major")])
    none = DatasetMatch(accession="GSE4", source="GEO", judge={"verdict": "irrelevant"})
    enrichment = {"shared_pathways": [{"sources_significant": ["GSE3", "GSE1"]}], "surfaced": []}
    classify_tiers(_spec(), [direct, proxy, pathway, none], enrichment)
    assert direct.match_tier == "direct"
    assert proxy.match_tier == "proxy"
    assert pathway.match_tier == "pathway"    # irrelevant on metadata but a pathway match
    assert none.match_tier == "none"


def _de(acc, rows, source_db="GEO", tier="primary"):
    t = pd.DataFrame(rows, columns=["gene", "log2fc", "pvalue", "padj"])
    return DEResult(accession=acc, contrast="c", n_case=3, n_control=3, table=t, status="ok",
                    source_db=source_db, tier=tier)


def test_heatmap_matrix_shape_and_cells():
    d1 = _de("GSE1", [("TP53", 2.0, 1e-5, 1e-4), ("BRCA1", -1.5, 1e-4, 1e-3), ("EGFR", 0.2, 0.5, 0.9)])
    d2 = _de("GSE2", [("TP53", 1.8, 1e-4, 1e-3), ("MYC", 3.0, 1e-6, 1e-5)])
    matches = [DatasetMatch(accession="GSE1", source="GEO", match_tier="direct"),
               DatasetMatch(accession="GSE2", source="GEO", match_tier="proxy")]
    h = build_heatmap([d1, d2], matches, user_genes=["TP53", "BRCA1", "MYC"])
    assert h is not None
    assert {g["gene"] for g in h["genes"]} == {"TP53", "BRCA1", "MYC"}
    assert [d["accession"] for d in h["datasets"]] == ["GSE1", "GSE2"]  # direct before proxy
    # TP53 present in both; matrix cell holds the log2FC (full table, not only significant)
    tp53_row = h["matrix"][[g["gene"] for g in h["genes"]].index("TP53")]
    assert tp53_row[0] == 2.0 and tp53_row[1] == 1.8
    # MYC only in GSE2 -> None in GSE1 column
    myc_row = h["matrix"][[g["gene"] for g in h["genes"]].index("MYC")]
    assert myc_row[0] is None and myc_row[1] == 3.0


def test_heatmap_resolves_mixed_id_table(monkeypatch):
    """A study naming some genes by symbol and others by Ensembl id (a real GEO pattern:
    circRNA rows by symbol, linear-gene rows by Ensembl) must not silently drop the
    Ensembl-keyed genes from the heatmap. Reverse-map is stubbed to stay offline."""
    import omics_oracle.heatmap as hm
    known = {"MYC": "ENSG00000136997", "VEGFA": "ENSG00000112715"}

    def fake_map_symbols_to(syms, target, species="human"):
        if target != "ensembl":
            return {}
        return {s.upper(): known[s.upper()] for s in syms if s.upper() in known}

    monkeypatch.setattr(hm, "map_symbols_to", fake_map_symbols_to)
    d = _de("GSEMIX", [("FOO-AS1", 2.0, 1e-5, 1e-4),                # symbol-named row
                       ("ENSG00000136997.14", 1.4, 1e-4, 1e-3),     # MYC by Ensembl id
                       ("ENSG00000112715.20", -1.1, 1e-3, 1e-2)])   # VEGFA by Ensembl id
    matches = [DatasetMatch(accession="GSEMIX", source="GEO", match_tier="direct")]
    h = build_heatmap([d], matches, user_genes=["MYC", "VEGFA", "FOO-AS1"])
    cell = {g["gene"]: row[0] for g, row in zip(h["genes"], h["matrix"])}
    assert cell["MYC"] == 1.4        # resolved via Ensembl reverse-map despite the mixed table
    assert cell["VEGFA"] == -1.1
    assert cell["FOO-AS1"] == 2.0    # plain symbol row still resolves by direct hit


def test_narrative_validator_flags_invented_gene():
    allowed = {"TP53", "BRCA1", "MYC"}
    text = "TP53 and BRCA1 were consistent, and the invented gene ZZZ9 appeared strongly."
    flagged = validate_narrative(text, allowed)
    assert "ZZZ9" in flagged
    assert "TP53" not in flagged and "BRCA1" not in flagged


def test_validator_ignores_domain_acronyms():
    allowed = {"TP53"}
    text = "RNA-seq DE analysis via KEGG and MAPK signaling in GSE12345 showed TP53 up."
    flagged = validate_narrative(text, allowed)
    assert flagged == []  # RNA, DE, KEGG, MAPK, GSE12345 are not flagged as invented genes


def test_collect_allowed_from_comparison():
    comparison = {"genes": [{"gene": "TP53"}, {"gene": "MYC"}]}
    allowed = collect_allowed_genes(comparison)
    assert allowed == {"TP53", "MYC"}
