"""Validate the pydeseq2 DE path on a small, deterministic count matrix.

We plant a handful of genes with a real fold change and confirm reanalyze() recovers them,
proving the counts -> DE seam works end-to-end (no network needed).
"""
import numpy as np
import pandas as pd

from omics_oracle.types import CountMatrix, ExperimentSpec
from omics_oracle.reanalyze import reanalyze, compute_qc


def _make_counts(seed=0):
    rng = np.random.default_rng(seed)
    n_genes, n_per = 200, 6
    genes = [f"G{i:03d}" for i in range(n_genes)]
    ctrl = [f"GSM_c{i}" for i in range(n_per)]
    case = [f"GSM_t{i}" for i in range(n_per)]
    base = rng.integers(50, 500, size=n_genes)
    def draw(mult):
        return np.vstack([rng.poisson(base * m) for m in mult]).T  # genes x samples
    ctrl_mat = draw(np.ones(n_per))
    up = np.ones(n_genes); up[:10] = 4.0     # 10 genes strongly up in case
    down = np.ones(n_genes); down[10:20] = 0.25  # 10 genes strongly down
    case_mat = draw(np.ones(n_per)) * (up * down)[:, None]
    counts = np.hstack([ctrl_mat, case_mat]).astype(int)
    df = pd.DataFrame(counts, index=genes, columns=ctrl + case)
    cm = CountMatrix(accession="GSE_TEST", genes=genes, samples=ctrl + case,
                     counts=df, provenance="synthetic fixture")
    groups = {**{s: "control" for s in ctrl}, **{s: "case" for s in case}}
    return cm, groups


def test_qc_summary_shape():
    cm, _ = _make_counts()
    qc = compute_qc(cm)
    assert qc.n_samples == 12
    assert qc.median_library_size > 0


def test_de_recovers_planted_signal():
    cm, groups = _make_counts()
    spec = ExperimentSpec(raw_prompt="x", omics_type="rna-seq",
                          contrast_case="case", contrast_control="control")
    res = reanalyze(cm, spec, groups=groups)
    assert res.status == "ok", res.note
    assert res.n_case == 6 and res.n_control == 6
    t = res.table.set_index("gene")
    # planted up genes G000..G009 should be significant & positive
    up_hits = t.loc[[f"G{i:03d}" for i in range(10)]]
    assert (up_hits["log2fc"] > 1).mean() >= 0.8
    assert (up_hits["padj"] < 0.05).mean() >= 0.8
    # a null gene should generally not be called
    assert res.n_de() >= 15  # ~20 planted, allow slack


def test_de_skipped_without_groups():
    cm, _ = _make_counts()
    spec = ExperimentSpec(raw_prompt="x", omics_type="rna-seq")  # no contrast, no meta
    res = reanalyze(cm, spec, groups=None)
    assert res.status == "skipped"
    assert res.qc is not None  # QC still reported
