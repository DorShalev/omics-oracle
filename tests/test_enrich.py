"""Enrichment aggregation: shared pathways + the 'sub-threshold here / strong there' surfacing.

Enrichr is mocked so this is offline and deterministic.
"""
import pandas as pd

import omics_oracle.enrich as enrich
from omics_oracle.types import DEResult


def _de(acc, genes, tier="primary", source_db="GEO"):
    rows = [(g, 2.0, 1e-5, 1e-4) for g in genes]
    t = pd.DataFrame(rows, columns=["gene", "log2fc", "pvalue", "padj"])
    return DEResult(accession=acc, contrast="c", n_case=3, n_control=3, table=t, status="ok",
                    tier=tier, source_db=source_db)


def test_shared_and_surfaced(monkeypatch):
    # canned Enrichr results per gene-list content
    def fake_enrichr(genes, library, timeout=30.0):
        gs = set(genes)
        # "p53 pathway": significant in GSE1 & GSE2, present-but-subthreshold in USER
        if gs == {"CDKN1A", "MDM2", "BAX"}:                 # GSE1
            return {"p53 pathway": {"p": 1e-6, "padj": 1e-4, "combined": 50, "genes": ["CDKN1A", "MDM2", "BAX"]}}
        if gs == {"CDKN1A", "MDM2", "GADD45A"}:             # GSE2
            return {"p53 pathway": {"p": 1e-5, "padj": 1e-3, "combined": 40, "genes": ["CDKN1A", "MDM2", "GADD45A"]}}
        if gs == {"CDKN1A", "MDM2", "TP53", "ATM"}:         # USER: many p53 genes, NOT significant
            return {"p53 pathway": {"p": 0.2, "padj": 0.4, "combined": 5, "genes": ["CDKN1A", "MDM2", "TP53", "ATM"]}}
        return {}
    monkeypatch.setattr(enrich, "run_enrichr", fake_enrichr)

    des = [_de("GSE1", ["CDKN1A", "MDM2", "BAX"]),
           _de("GSE2", ["CDKN1A", "MDM2", "GADD45A"], tier="clinical", source_db="ARCHS4")]
    out = enrich.enrich_across_datasets(des, user_genes=["CDKN1A", "MDM2", "TP53", "ATM"])

    # shared: p53 pathway significant in GSE1 + GSE2
    shared = {p["term"]: p for p in out["shared_pathways"]}
    assert "p53 pathway" in shared and shared["p53 pathway"]["n_significant"] == 2

    # surfaced: significant in GSE1/GSE2 but sub-threshold in USER despite 4 genes moving
    surf = {p["term"]: p for p in out["surfaced"]}
    assert "p53 pathway" in surf
    subs = surf["p53 pathway"]["subthreshold_in"]
    assert any(s["source"] == "USER" and s["overlap_n"] == 4 for s in subs)


def test_offline_returns_note(monkeypatch):
    monkeypatch.setattr(enrich, "run_enrichr", lambda g, lib, timeout=30.0: None)
    out = enrich.enrich_across_datasets([_de("GSE1", ["A", "B", "C"])])
    assert out["shared_pathways"] == [] and "unreachable" in out["note"]
