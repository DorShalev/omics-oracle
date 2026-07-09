"""Projector: turn matched studies + concordance into a forward-looking sanity check.

Everything here is grounded in the matched datasets or in well-established assay caveats.
Where evidence is thin, the projection says so (feasibility="thin-evidence") rather than
manufacturing a confident number.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from .types import Concordance, DatasetMatch, DEResult, ExperimentSpec, Projection

# Assay-specific QC pitfalls surfaced regardless of matches (well-established priors).
_ASSAY_PITFALLS = {
    "rna-seq": [
        "Check for rRNA / globin carryover, especially in blood/plasma inputs.",
        "Verify library size and 3' bias; low-input EV/exosome RNA often yields shallow, "
        "high-duplication libraries.",
        "Normalization choice (median-of-ratios vs TMM vs CPM) changes DE calls — fix one and report it.",
        "Batch/lane confounding with condition is the most common false-signal source; block or randomize.",
    ],
    "scrna-seq": [
        "Ambient RNA and doublets inflate apparent DE; run doublet detection and empty-drop filtering.",
        "Cell-type composition shifts can masquerade as DE — analyze within matched cell types.",
    ],
    "proteomics": [
        "Missing-value structure (MNAR vs MCAR) drives DE; pick an imputation strategy deliberately.",
        "Normalize for total intensity / median; watch for batch effects across MS runs.",
    ],
    "atac-seq": [
        "TSS enrichment and fragment-size periodicity are the key QC gates before differential accessibility.",
    ],
}

_TISSUE_PITFALLS = {
    "plasma": "Plasma/EV RNA is low-abundance and hemolysis-sensitive; expect high technical variance.",
    "exosome": "Exosome/EV preps are frequently low-quality; RNA integrity and yield vary widely.",
    "ev": "Extracellular-vesicle RNA is low-input; shallow libraries and high dropout are common.",
    "serum": "Serum RNA is even lower-input than plasma; platelet contamination confounds signal.",
    "ffpe": "FFPE RNA is degraded; use appropriate library prep and expect elevated duplication.",
}


def project(spec: ExperimentSpec, matches: list[DatasetMatch],
            de_results: list[DEResult], concordance: Optional[Concordance]) -> Projection:
    proj = Projection()

    # --- expected QC pitfalls (assay + tissue priors) ---
    pitfalls = list(_ASSAY_PITFALLS.get(spec.omics_type, []))
    if spec.tissue:
        low = spec.tissue.lower()
        for key, msg in _TISSUE_PITFALLS.items():
            if key in low and msg not in pitfalls:
                pitfalls.append(msg)
    # QC flags actually observed in the reanalyzed prior studies carry the most weight
    for d in de_results:
        if d.qc:
            for f in d.qc.flags:
                obs = f"Observed in {d.accession}: {f}"
                if obs not in pitfalls:
                    pitfalls.append(obs)
    proj.qc_pitfalls = pitfalls

    # --- expected #DE and effect-size prior from reanalyzed studies ---
    n_de_counts = [d.n_de() for d in de_results if d.status == "ok"]
    if n_de_counts:
        lo, hi = int(min(n_de_counts)), int(max(n_de_counts))
        proj.expected_n_de = (lo, hi)
        proj.effect_size_prior = _effect_prior(de_results)
    else:
        # explain WHY no study yielded DE, from what actually happened per dataset
        reasons = []
        norm = sum(1 for d in de_results if "normalized" in (d.note or ""))
        notab = sum(1 for d in de_results if "no usable per-study counts" in (d.note or ""))
        nogrp = sum(1 for d in de_results if "grouping" in (d.note or "").lower()
                    or "replicates" in (d.note or "").lower())
        if norm:
            reasons.append(f"{norm} had only normalized tables (TPM/FPKM/log), unusable for pydeseq2")
        if notab:
            reasons.append(f"{notab} had no usable per-study counts table")
        if nogrp:
            reasons.append(f"{nogrp} had no assignable case/control grouping")
        detail = ("; ".join(reasons)) if reasons else "no reanalyzable raw-count table was found"
        proj.caveats.append(f"No prior study yielded raw counts for DE ({detail}); DE-count and "
                            "effect-size priors are not estimated. Try raising --reanalyze-top to "
                            "reach studies deeper in the ranking that may deposit raw counts.")

    # --- power / feasibility read ---
    proj.power_read = _power_read(spec, de_results)

    # --- overall feasibility label + evidence-strength caveats ---
    usable = [m for m in matches if not m.blocking]
    reanalyzed = sum(1 for d in de_results if d.status == "ok")
    if reanalyzed >= 2 and concordance and concordance.genes:
        proj.feasibility = "encouraging"
    elif usable:
        proj.feasibility = "cautionary"
    else:
        proj.feasibility = "thin-evidence"

    if not usable:
        proj.caveats.append("No matched dataset cleared the blocking-mismatch filter "
                            "(organism/omics); treat any projection as weakly supported.")
    if matches and all(m.blocking for m in matches):
        proj.caveats.append("Every candidate had a blocking mismatch — the closest public "
                            "studies differ in organism or assay from the plan.")
    return proj


def _effect_prior(de_results: list[DEResult]) -> str:
    absfc = []
    for d in de_results:
        if d.status == "ok" and d.table is not None and len(d.table):
            t = d.table
            sig = t[(t["padj"].fillna(1.0) < 0.05) & (t["log2fc"].abs() >= 1.0)]
            if len(sig):
                absfc.append(float(np.median(sig["log2fc"].abs())))
    if not absfc:
        return "no significant genes in reanalyzed studies at padj<0.05, |log2fc|>=1"
    med = float(np.median(absfc))
    return f"median |log2FC| of significant genes ≈ {med:.2f} across reanalyzed studies"


def _power_read(spec: ExperimentSpec, de_results: list[DEResult]) -> str:
    n = spec.n_per_group
    if not n:
        return ("Planned N per group not specified; specify it for a power read. As a rule of "
                "thumb, bulk RNA-seq needs ≈6–12 per group to detect moderate effects reliably.")
    reanalyzed = [d for d in de_results if d.status == "ok"]
    prior_ns = [min(d.n_case, d.n_control) for d in reanalyzed if d.n_case and d.n_control]
    base = f"Planned N = {n} per group. "
    if n < 3:
        base += "Below the practical floor for DESeq2 dispersion estimation (need ≥3, ideally ≥6). "
    elif n < 6:
        base += "Adequate only for large effects; underpowered for subtle transcriptional changes. "
    else:
        base += "Reasonable for moderate effects in bulk RNA-seq. "
    if prior_ns:
        base += (f"Comparable prior studies used as few as {min(prior_ns)} per group, so a "
                 f"detectable signal at N={n} is plausible if effect sizes match.")
    return base.strip()
