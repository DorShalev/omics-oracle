"""Per-dataset reanalysis: QC summary + differential expression via pydeseq2.

DE runs only when samples can be assigned to case/control groups (explicitly, or derived
from sample metadata using the spec's contrast terms). If no grouping can be established we
return a QC-only result with ``status="skipped"`` — we never invent a contrast.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from .types import CountMatrix, DEResult, ExperimentSpec, QCSummary


def compute_qc(counts: CountMatrix) -> QCSummary:
    df = counts.counts  # genes x samples
    lib = df.sum(axis=0)
    detected = (df > 0).sum(axis=0)
    zero_frac = (df == 0).mean(axis=0)
    flags: list[str] = []
    med_lib = float(np.median(lib)) if len(lib) else 0.0
    if med_lib and med_lib < 1_000_000:
        flags.append("low median library size (<1M counts) — typical of low-input / EV / exosome preps")
    if len(detected) and float(np.median(detected)) < 10_000:
        flags.append("few detected genes (<10k) — possible degraded/low-complexity libraries")
    return QCSummary(
        n_samples=df.shape[1],
        n_genes_detected=int((df.sum(axis=1) > 0).sum()),
        median_library_size=med_lib,
        median_zero_fraction=round(float(np.median(zero_frac)) if len(zero_frac) else 0.0, 3),
        flags=flags,
    )


def reanalyze(counts: CountMatrix, spec: ExperimentSpec,
              groups: Optional[dict[str, str]] = None,
              min_group: int = 2, llm=None, group_consensus: int = 1,
              verify_samples: bool = True) -> DEResult:
    """QC always; DE when a verified case/control selection is available.

    Sample selection (which samples are case / control / excluded) is delegated to
    :func:`omics_oracle.selection.select_samples`, which proposes labels, VERIFIES them with an
    independent judge, and records the full audit trail on ``DEResult.sample_selection``. Pass
    ``groups`` to fix the selection explicitly (skips the LLM). DE method is chosen by data type:
    raw counts -> pydeseq2, normalized -> Welch-on-log; unknown -> skipped.

    The DE result inherits the counts' repository labels (source_db / tier) so the report can
    show whether a study's DE came from the GEO per-study table (primary) or ARCHS4 (clinical).
    """
    res = _reanalyze_impl(counts, spec, groups, min_group, llm, group_consensus, verify_samples)
    res.source_db = getattr(counts, "source_db", "GEO")
    res.tier = getattr(counts, "tier", "primary")
    return res


def _reanalyze_impl(counts: CountMatrix, spec: ExperimentSpec,
                    groups: Optional[dict[str, str]], min_group: int, llm,
                    group_consensus: int = 1, verify_samples: bool = True) -> DEResult:
    qc = compute_qc(counts)

    # P0 honesty gate: pick the DE method by data type, never run the WRONG method.
    #   raw counts     -> pydeseq2 (Wald)               [method="deseq2"]
    #   normalized      -> Welch t-test on log2 values   [method="welch-log2"]
    #   unknown/unconfirmed -> skip (we won't guess which method is valid)
    # Normalized studies are NOT discarded: analysing each study case-vs-its-own-control still
    # yields a comparable DE gene list — we just label the method so cross-study comparison is
    # honest about mixing DESeq2 and t-test effect sizes.
    dtype = getattr(counts, "data_type", "raw")
    if dtype not in ("raw", "normalized"):
        return DEResult(accession=counts.accession, contrast="(none)", n_case=0, n_control=0,
                        qc=qc, status="skipped",
                        note=f"{counts.accession} table data-type is '{dtype}' "
                             f"({counts.source_file or 'source n/a'}); could not confirm raw vs "
                             "normalized, so DE was skipped rather than run with the wrong method")

    # propose + verify which samples are case / control / excluded (full trail on the result)
    from .selection import select_samples
    audit = select_samples(counts, spec, llm=llm, consensus=group_consensus,
                           verify=verify_samples, explicit_groups=groups)
    sel = audit.to_dict()
    sel_groups = audit.groups()

    if not sel_groups:
        return DEResult(accession=counts.accession, contrast="(none)", n_case=0, n_control=0,
                        qc=qc, status="skipped", sample_selection=sel,
                        note=_drop_note(audit))

    case = [s for s, g in sel_groups.items() if g == "case" and s in counts.counts.columns]
    control = [s for s, g in sel_groups.items() if g == "control" and s in counts.counts.columns]
    if len(case) < min_group or len(control) < min_group:
        return DEResult(accession=counts.accession,
                        contrast=f"{spec.contrast_case or 'case'} vs {spec.contrast_control or 'control'}",
                        n_case=len(case), n_control=len(control), qc=qc, status="skipped",
                        sample_selection=sel,
                        note=f"insufficient replicates after sample selection (case={len(case)}, "
                             f"control={len(control)}, need >= {min_group} each)")

    contrast = f"{spec.contrast_case or 'case'} vs {spec.contrast_control or 'control'}"
    try:
        if dtype == "raw":
            table, method = _run_pydeseq2(counts.counts, case, control), "deseq2"
        else:
            table, method = _run_normalized_de(counts.counts, case, control), "welch-log2"
    except Exception as e:  # keep the pipeline honest about failure rather than faking a table
        return DEResult(accession=counts.accession, contrast=contrast,
                        n_case=len(case), n_control=len(control), qc=qc, status="error",
                        sample_selection=sel, note=f"DE failed ({dtype}): {type(e).__name__}: {e}")

    note = ("DE via pydeseq2 (Wald test, case vs control)" if method == "deseq2" else
            "DE via Welch t-test on log2 normalized values (case vs control) — this table is "
            "normalized, not raw counts, so effect sizes are log2 ratios of normalized "
            "expression and significance uses a t-test, not DESeq2")
    return DEResult(accession=counts.accession, contrast=contrast,
                    n_case=len(case), n_control=len(control), table=table, qc=qc, status="ok",
                    method=method, sample_selection=sel, note=note)


def _drop_note(audit) -> str:
    """Honest skip reason built from the sample-selection audit."""
    if audit.outcome == "dropped" and audit.verdict is not None:
        concerns = "; ".join(audit.verdict.concerns) if audit.verdict.concerns else audit.verdict.reason
        return ("sample selection was rejected by the verification judge and dropped rather than "
                f"run on a questionable grouping — {concerns or 'no case/control contrast confirmed'}")
    return ("no case/control selection could be established for this study; QC reported, "
            "DE skipped (no contrast invented)")


def _run_pydeseq2(counts_gxs: pd.DataFrame, case: list[str], control: list[str]) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    samples = control + case  # control first so it is the reference level
    counts_sxg = counts_gxs[samples].T                      # samples x genes
    counts_sxg = counts_sxg.round().astype(int)
    counts_sxg = counts_sxg.loc[:, counts_sxg.sum(axis=0) >= 10]   # drop near-zero genes

    metadata = pd.DataFrame(
        {"condition": (["control"] * len(control)) + (["case"] * len(case))},
        index=samples,
    )

    # The explicit contrast below fixes numerator/denominator, so a reference level is not
    # needed (ref_level is deprecated and inert in pydeseq2 0.5.x).
    dds = DeseqDataSet(counts=counts_sxg, metadata=metadata, design="~condition", quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "case", "control"], quiet=True)
    stats.summary()
    res = stats.results_df.rename(columns={"log2FoldChange": "log2fc"})
    out = res[["log2fc", "pvalue", "padj"]].copy()
    out.index.name = "gene"
    return out.reset_index()


def _run_normalized_de(expr_gxs: pd.DataFrame, case: list[str], control: list[str]) -> pd.DataFrame:
    """DE for NORMALIZED tables (TPM/FPKM/log) — Welch t-test per gene on log2 values.

    DESeq2 is invalid here (it models raw-count dispersions), so we use a distribution-light
    per-gene test on the log scale: log2FC = mean(case) - mean(control) in log2 space, a Welch
    (unequal-variance) t-test for significance, and Benjamini-Hochberg FDR. Effect sizes are
    log2 ratios of *normalized* expression, not shrunk DESeq2 estimates — labelled as such.
    """
    from scipy import stats
    from scipy.stats import false_discovery_control

    sub = expr_gxs[control + case].astype(float)
    # Work in log2 space. If the table looks linear (large dynamic range) log-transform it;
    # if it already looks log-scaled (compact range) use it as-is, so we never double-log.
    hi = float(np.nanpercentile(sub.values, 99)) if sub.size else 0.0
    logged = sub if hi <= 30 else np.log2(sub + 1.0)
    keep = logged.notna().all(axis=1) & (logged.std(axis=1) > 0)   # drop flat/absent genes
    logged = logged[keep]
    if logged.empty:
        raise ValueError("no genes with usable variation after log-transform")

    ctrl_m = logged[control].to_numpy(dtype=float)
    case_m = logged[case].to_numpy(dtype=float)
    _, p = stats.ttest_ind(case_m, ctrl_m, axis=1, equal_var=False)
    p = np.where(np.isnan(p), 1.0, p)
    log2fc = case_m.mean(axis=1) - ctrl_m.mean(axis=1)
    padj = false_discovery_control(p, method="bh")
    return pd.DataFrame({"gene": logged.index.astype(str), "log2fc": log2fc,
                         "pvalue": p, "padj": padj}).reset_index(drop=True)


# Sample-selection helpers (heuristic + LLM propose/verify/consensus) live in selection.py.
