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
              min_group: int = 2) -> DEResult:
    """QC always; DE when a case/control grouping is available.

    ``groups`` maps sample id -> "case"|"control". If None, we try to derive it from
    ``counts.sample_meta`` using the spec's contrast terms.
    """
    qc = compute_qc(counts)

    if groups is None:
        groups = _derive_groups(counts, spec)

    if not groups:
        return DEResult(accession=counts.accession, contrast="(none)", n_case=0, n_control=0,
                        qc=qc, status="skipped",
                        note="no case/control grouping could be established from metadata; "
                             "QC reported, DE skipped (no contrast invented)")

    case = [s for s, g in groups.items() if g == "case" and s in counts.counts.columns]
    control = [s for s, g in groups.items() if g == "control" and s in counts.counts.columns]
    if len(case) < min_group or len(control) < min_group:
        return DEResult(accession=counts.accession,
                        contrast=f"{spec.contrast_case or 'case'} vs {spec.contrast_control or 'control'}",
                        n_case=len(case), n_control=len(control), qc=qc, status="skipped",
                        note=f"insufficient replicates for DE (case={len(case)}, control={len(control)}, "
                             f"need >= {min_group} each)")

    try:
        table = _run_pydeseq2(counts.counts, case, control)
    except Exception as e:  # keep the pipeline honest about failure rather than faking a table
        return DEResult(accession=counts.accession,
                        contrast=f"{spec.contrast_case or 'case'} vs {spec.contrast_control or 'control'}",
                        n_case=len(case), n_control=len(control), qc=qc, status="error",
                        note=f"pydeseq2 DE failed: {type(e).__name__}: {e}")

    return DEResult(accession=counts.accession,
                    contrast=f"{spec.contrast_case or 'case'} vs {spec.contrast_control or 'control'}",
                    n_case=len(case), n_control=len(control), table=table, qc=qc, status="ok",
                    note="DE via pydeseq2 (Wald test, case vs control)")


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


def _derive_groups(counts: CountMatrix, spec: ExperimentSpec) -> Optional[dict[str, str]]:
    """Best-effort case/control assignment from sample metadata using contrast terms."""
    meta = counts.sample_meta
    if meta is None or spec.contrast_case is None or spec.contrast_control is None:
        return None
    case_terms = _terms(spec.contrast_case)
    ctrl_terms = _terms(spec.contrast_control)
    groups: dict[str, str] = {}
    for sample in counts.counts.columns:
        if sample not in meta.index:
            continue
        blob = " ".join(str(v) for v in meta.loc[sample].values).lower()
        if any(t in blob for t in case_terms):
            groups[sample] = "case"
        elif any(t in blob for t in ctrl_terms):
            groups[sample] = "control"
    return groups or None


def _terms(s: str) -> list[str]:
    return [t for t in s.lower().replace("/", " ").split() if len(t) > 2]
