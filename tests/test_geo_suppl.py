"""GeoSupplSource: URL derivation, file selection, split-merge, and a monkeypatched e2e."""
import gzip

import numpy as np
import pandas as pd

from omics_oracle.sources.geo_suppl import GeoSupplSource
from omics_oracle.types import DatasetMatch, ExperimentSpec
from omics_oracle.reanalyze import reanalyze


def test_suppl_url_derivation():
    s = GeoSupplSource()
    assert s._suppl_url("GSE313035").endswith("/GSE313nnn/GSE313035/suppl/")
    assert s._suppl_url("GSE100").endswith("/GSEnnn/GSE100/suppl/")
    assert s._suppl_url("GSE1").endswith("/GSEnnn/GSE1/suppl/")


def test_best_group_prefers_raw_over_normalized():
    s = GeoSupplSource()
    cands = [("GSE1_TPM.txt.gz", "u1"), ("GSE1_raw_counts.tsv.gz", "u2")]
    group = s._best_group(cands)
    assert group[0][0] == "GSE1_raw_counts.tsv.gz"


def test_best_group_merges_condition_splits():
    s = GeoSupplSource()
    cands = [("GSE1_TPM_ctrl.txt.gz", "uc"), ("GSE1_TPM_exp.txt.gz", "ue")]
    group = s._best_group(cands)
    assert len(group) == 2  # split set kept together for horizontal merge


def _gz_csv(df, sep=","):
    return gzip.compress(df.to_csv(sep=sep).encode())


def test_fetch_counts_end_to_end_monkeypatched(monkeypatch):
    # two per-condition split files sharing a gene index -> merge -> raw counts -> DE
    genes = [f"ENSG{i:04d}" for i in range(80)]
    rng = np.random.default_rng(0)
    ctrl = pd.DataFrame(rng.integers(50, 500, (80, 4)), index=genes,
                        columns=[f"ctrl_{i}" for i in range(4)])
    exp = pd.DataFrame(rng.integers(50, 500, (80, 4)), index=genes,
                       columns=[f"treat_{i}" for i in range(4)])
    exp.iloc[:8] *= 5  # a real up-signal in the treated arm

    files = [("GSE9_raw_counts_ctrl.tsv.gz", "http://x/ctrl"),
             ("GSE9_raw_counts_treat.tsv.gz", "http://x/treat")]
    blobs = {"http://x/ctrl": _gz_csv(ctrl, sep="\t"),
             "http://x/treat": _gz_csv(exp, sep="\t")}

    s = GeoSupplSource()
    monkeypatch.setattr(s, "_list_suppl", lambda acc: files)
    monkeypatch.setattr(s, "_download", lambda acc, fname, url: blobs[url])

    cm = s.fetch_counts(DatasetMatch(accession="GSE9", source="GEO"))
    assert cm is not None
    assert cm.data_type == "raw"
    assert cm.counts.shape == (80, 8)             # merged 4 + 4 samples
    assert "ctrl_0" in cm.samples and "treat_0" in cm.samples

    # grouping from column names works without an LLM (deterministic path)
    spec = ExperimentSpec(raw_prompt="x", omics_type="rna-seq",
                          contrast_case="treat", contrast_control="ctrl")
    res = reanalyze(cm, spec)
    assert res.status == "ok" and res.n_case == 4 and res.n_control == 4
    assert res.n_de() >= 5


def test_fetch_counts_normalized_is_analyzed_with_ttest(monkeypatch):
    rng = np.random.default_rng(1)
    genes = [f"G{i}" for i in range(60)]
    tpm = pd.DataFrame(np.abs(rng.random((60, 6))) * 3 + 0.1,
                       index=genes, columns=[f"s{i}" for i in range(6)])
    # plant a clean, low-variance signal in the first 5 genes: high in s3-s5, low in s0-s2
    tpm.iloc[:5, :3] = 0.5 + rng.random((5, 3)) * 0.05
    tpm.iloc[:5, 3:] = 5.0 + rng.random((5, 3)) * 0.05
    s = GeoSupplSource()
    monkeypatch.setattr(s, "_list_suppl", lambda acc: [("GSE8_TPM.tsv.gz", "http://x/tpm")])
    monkeypatch.setattr(s, "_download", lambda acc, fname, url: _gz_csv(tpm, sep="\t"))
    cm = s.fetch_counts(DatasetMatch(accession="GSE8", source="GEO"))
    assert cm is not None and cm.data_type == "normalized"
    assert any("normalized" in w for w in s.last_warnings)
    # normalized tables are NOT discarded: DE runs via the labelled t-test path, not DESeq2
    groups = {"s0": "control", "s1": "control", "s2": "control",
              "s3": "case", "s4": "case", "s5": "case"}
    res = reanalyze(cm, ExperimentSpec(raw_prompt="x", omics_type="rna-seq",
                                       contrast_case="case", contrast_control="control"),
                    groups=groups)
    assert res.status == "ok" and res.method == "welch-log2"
    assert res.n_de() >= 1 and "t-test" in res.note.lower()


def test_no_table_returns_none_with_reason(monkeypatch):
    s = GeoSupplSource()
    monkeypatch.setattr(s, "_list_suppl",
                        lambda acc: [("GSE7_RAW.bam", "u"), ("GSE7_clinical.txt.gz", "u2")])
    cm = s.fetch_counts(DatasetMatch(accession="GSE7", source="GEO"))
    assert cm is None
    assert any("no supplementary expression table" in w for w in s.last_warnings)
