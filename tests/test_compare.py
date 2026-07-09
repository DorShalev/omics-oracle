"""Cross-dataset comparison: per-gene detection, pairwise overlap, user-set overlap."""
import pandas as pd

from omics_oracle.types import DEResult
from omics_oracle.compare import build_comparison


def _de(acc, rows, source_db="GEO", tier="primary"):
    t = pd.DataFrame(rows, columns=["gene", "log2fc", "pvalue", "padj"])
    return DEResult(accession=acc, contrast="c", n_case=3, n_control=3, table=t, status="ok",
                    source_db=source_db, tier=tier)


def test_detection_pairwise_and_datasets():
    d1 = _de("GSE1", [("TP53", 2.0, 1e-5, 1e-4), ("BRCA1", -1.5, 1e-4, 1e-3),
                      ("EGFR", 1.2, 1e-4, 1e-3)])
    d2 = _de("GSE2", [("TP53", 1.8, 1e-4, 1e-3), ("BRCA1", -1.2, 1e-3, 1e-2),
                      ("MYC", 3.0, 1e-6, 1e-5)], source_db="ARCHS4", tier="clinical")
    c = build_comparison([d1, d2])
    # datasets carry their source/tier
    tiers = {d["accession"]: d["tier"] for d in c["datasets"]}
    assert tiers == {"GSE1": "primary", "GSE2": "clinical"}
    # TP53 + BRCA1 detected in both; EGFR/MYC in one
    genes = {g["gene"]: g for g in c["genes"]}
    assert genes["TP53"]["n_datasets"] == 2 and genes["TP53"]["up"] == 2
    assert genes["BRCA1"]["n_datasets"] == 2 and genes["BRCA1"]["down"] == 2
    assert genes["EGFR"]["n_datasets"] == 1
    # pairwise overlap = {TP53, BRCA1} = 2 shared
    p = c["pairwise"][0]
    assert p["shared"] == 2 and p["a_only"] == 1 and p["b_only"] == 1


def test_user_gene_list_overlap():
    d1 = _de("GSE1", [("TP53", 2.0, 1e-5, 1e-4), ("BRCA1", -1.5, 1e-4, 1e-3)])
    d2 = _de("GSE2", [("TP53", 1.8, 1e-4, 1e-3), ("MYC", 3.0, 1e-6, 1e-5)])
    c = build_comparison([d1, d2], user_genes=["TP53", "KRAS"])
    assert c["user"]["n_genes"] == 2
    vs = {v["accession"]: v for v in c["user"]["vs_datasets"]}
    assert vs["GSE1"]["shared"] == 1 and "TP53" in vs["GSE1"]["shared_genes"]
    assert vs["GSE2"]["shared"] == 1
    # the gene list is scoped to the user's genes
    listed = {g["gene"] for g in c["genes"]}
    assert listed == {"TP53", "KRAS"}


def test_none_when_no_ok_results():
    assert build_comparison([DEResult(accession="GSE1", contrast="c", n_case=0, n_control=0,
                                      status="skipped")]) is None
