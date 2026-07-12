# Omics Oracle — experiment sanity check

> _Human pancreatic islet RNA-seq, type 2 diabetes vs non-diabetic donors_

**Parsed spec** (`llm` parser): organism=**Homo sapiens**, tissue=**pancreatic islet**, omics=**rna-seq**, contrast=**type 2 diabetes donors vs non-diabetic donors**, platform=**—**, N/group=**—**.

**Warnings:** LLM backend: claude-cli (spec parsing + narrative phrasing only; numbers come from data, not the model).; GSE335464: no supplementary expression table found (files: GSE335464_README.md.txt, GSE335464_processed_data.tar.gz); GSE285253: assembled a 26870×12 matrix from 12 per-sample files inside GSE285253_RAW.tar (raw); GSE333721: single-cell RNA-seq deposition (.mtx/.h5/.rds) — set aside; bulk DE not applicable (single-cell analysis is planned); GSE305932: single-cell RNA-seq (10x matrix / .h5 / .rds in RAW.tar) — set aside; bulk DE not applicable (single-cell analysis is planned); GSE310965: supplementary table is normalized (GSE310965_Salmon_Batch-Corrected_Normalized_Counts_PIC_and_Vehicle.xlsx); DE will use a Welch t-test on log values (not DESeq2, which needs raw counts); GSE332956: single-cell RNA-seq (10x matrix / .h5 / .rds in RAW.tar) — set aside; bulk DE not applicable (single-cell analysis is planned); GSE301598: single-cell RNA-seq (10x matrix / .h5 / .rds in RAW.tar) — set aside; bulk DE not applicable (single-cell analysis is planned)

## Narrative
Prior public data offers only weak support for this plan so far: several human pancreatic-islet RNA-seq studies rank as very close matches to the proposed type 2 diabetes vs. non-diabetic donor comparison, but none could actually be reanalyzed into that contrast. The closest datasets by accession are GSE335464 (n=28), GSE308718 (n=6), GSE312607 (n=12), GSE285253 (n=12), and GSE333721 (n=6). The main mismatches are substantive: GSE335464 had no usable per-study counts table; GSE285253 contains only iPSC-derived HNF1A-MODY (monogenic diabetes, not T2D donor islets); GSE310965 is a vehicle-vs-poly(I:C) treatment design; and GSE307751 is a MYT1/MYT3 knockdown experiment — none carry donor T2D-vs-non-diabetic labels. Additional candidates (GSE333721, GSE305932, GSE332956, GSE301598) are single-cell studies set aside from bulk differential expression. Because no study yielded raw counts for a valid contrast, no differential-expression count or effect-size priors could be estimated, and general QC risks (normalization choice, batch confounding, and the low library sizes flagged in GSE310965) remain unaddressed. Note also that no planned N per group was specified, so no power read is possible beyond the general rule that bulk RNA-seq typically needs roughly 6–12 samples per group. Overall feasibility is cautionary: the biological question is well-represented in principle, but the retrieved datasets do not yet provide a directly matching, reanalyzable contrast, so deeper searching for donor-status–annotated islet cohorts is warranted before drawing conclusions.

## Matched datasets (ranked)
| Rank | Accession | Sim | Tier | Organism | Omics | N | Blocking? | Judge | Mismatches |
|---|---|---|---|---|---|---|---|---|---|
| 1 | [GSE335464](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335464) | 1.0 | proxy | Homo sapiens | rna-seq | 28 | — | partial | — |
| 2 | [GSE308718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308718) | 1.0 | none | Homo sapiens | rna-seq | 6 | — | irrelevant | — |
| 3 | [GSE312607](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312607) | 1.0 | none | Homo sapiens | rna-seq | 12 | — | irrelevant | — |
| 4 | [GSE285253](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285253) | 1.0 | proxy | Homo sapiens | rna-seq | 12 | — | partial | — |
| 5 | [GSE333721](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333721) | 1.0 | proxy | Homo sapiens | rna-seq | 6 | — | — | — |
| 6 | [GSE305932](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305932) | 1.0 | proxy | Homo sapiens | rna-seq | 2 | — | — | — |
| 7 | [GSE310965](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310965) | 1.0 | proxy | Homo sapiens | rna-seq | 60 | — | — | — |
| 8 | [GSE307751](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307751) | 1.0 | proxy | Homo sapiens | rna-seq | 27 | — | — | — |
| 9 | [GSE332956](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE332956) | 0.8 | proxy | Homo sapiens | rna-seq | 8 | — | — | contrast(major) |
| 10 | [GSE301598](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301598) | 0.8 | proxy | Homo sapiens | rna-seq | 1 | — | — | contrast(major) |

**Relevance verification (LLM judge):**
- GSE335464: **partial** (conf 0.72) — Organism (human) and tissue (pancreatic islet) match, but the disease and contrast differ substantively: this study targets type 1 diabetes (GWAS effector-gene mechanisms), not type 2 diabetes, and its RNA-seq contrast is gene knockout vs wild-type in hESC-derived beta-like cells (plus single-cell baseline vs T1D-stimulating conditions) — not a T2D-vs-non-diabetic donor comparison. Overlaps in sys
- GSE308718: **irrelevant** (conf 0.9) — Although human RNA-seq with diabetes/beta-cell themes, the actual system is stem-cell-derived beta cells (SC-β) in vitro, not primary pancreatic islets; the contrast is macrophage polarization states (unpolarized/inflammatory/regenerative) during differentiation, not T2D vs non-diabetic donors; and the disease context is type 1 diabetes, not T2D. No shared biological system, contrast axis, or dono
- GSE312607: **irrelevant** (conf 0.86) — Right organism, tissue (islet) and assay (RNA-seq), but the biological system and contrast are wrong: it is a stem-cell-derived islet HNF1A-MODY (monogenic diabetes) mutant-vs-control comparison, not primary human pancreatic islets from T2D vs non-diabetic donors. Different disease entity (MODY, not type 2 diabetes) and different contrast axis (in vitro genotype vs donor disease status).
- GSE285253: **partial** (conf 0.72) — Matches organism (human), assay (RNA-seq), and islet context, but the biological system and contrast differ: it uses stem-cell-derived islets (not primary donor islets) and contrasts an HNF1A-MODY monogenic-diabetes mutation vs control, not type 2 diabetes vs non-diabetic donors. MODY is a distinct diabetes subtype, so this offers only tangential, not direct, prior evidence.

**Mismatch caveats:**
- contrast terms 'type 2 diabetes donors non-diabetic donors' not evident in GSE332956 metadata
- contrast terms 'type 2 diabetes donors non-diabetic donors' not evident in GSE301598 metadata

## Per-dataset reanalysis
- **GSE335464** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE285253** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case set is empty (0 samples) and control set is empty (0 samples) — a differential-expression contrast requires non-empty groups on both sides; nothing can be computed.; All 12 samples are iPSC-derived in-vitro cells (non-carrier or HNF1A-MODY, some purmorphamine-treated), none are primary T2D or non-diabetic donor islets — the exclusions are individually justified, b; The dataset as described appears to contain no primary pancreatic-islet donor samples matching either the T2D or non-diabetic-donor definition; the requested contrast does not exist in these GSMs, so ; HNF1A-MODY is monogenic diabetes, correctly not equated with T2D; but again this only confirms there is no case arm available.
  - QC: 12 samples, median lib size 28,105,480, median zero-fraction 0.218.
  - Sample selection (dropped): case=0, control=0, excluded=12 — judge ⚠️ disagreed (Case set is empty (0 samples) and control set is empty (0 samples) — a differential-expression contrast requires non-empty groups on both sides; nothing can be computed.; All 12 samples are iPSC-derived in-vitro cells (non-carrier or HNF1A-MODY, some purmorphamine-treated), none are primary T2D or non-diabetic donor islets — the exclusions are individually justified, b)
- **GSE333721** _[primary: GEO]_ — SKIPPED: single-cell RNA-seq study — set aside from bulk DE (single-cell analysis is planned)
- **GSE305932** _[primary: GEO]_ — SKIPPED: single-cell RNA-seq study — set aside from bulk DE (single-cell analysis is planned)
- **GSE310965** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case (0) and control (0) are both empty — a differential-expression contrast requires non-empty case and control groups.; All 40 samples belong to a VEH-vs-PIC (vehicle vs poly I:C dsRNA-mimic) islet stimulation experiment, i.e. a treatment-arm design, not a donor diabetes-status cohort.; The proposer's per-sample reasons are individually correct (none of these labels encode T2D vs non-diabetic status), which actually confirms this dataset does not contain the requested contrast.; This dataset is unsuitable for the T2D-vs-non-diabetic islet question; a different dataset with donor diabetes-status annotation is needed. (If donor T2D status exists as separate metadata not shown h
  - QC: 40 samples, median lib size 151,948, median zero-fraction 0.0. Flags: low median library size (<1M counts) — typical of low-input / EV / exosome preps
  - Sample selection (dropped): case=0, control=0, excluded=40 — judge ⚠️ disagreed (Case (0) and control (0) are both empty — a differential-expression contrast requires non-empty case and control groups.; All 40 samples belong to a VEH-vs-PIC (vehicle vs poly I:C dsRNA-mimic) islet stimulation experiment, i.e. a treatment-arm design, not a donor diabetes-status cohort.)
- **GSE307751** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case set is empty (0 samples) and control set is empty (0 samples); a differential-expression contrast requires both arms to be non-empty.; All 26 samples are MYT1/MYT3 shRNA-knockdown vs Sham experiments — a gene-perturbation design, not donor disease status. The dataset contains no T2D or non-diabetic donor arm, so the requested contras; The per-sample exclusion reasons are individually correct (knockdown/Sham samples genuinely lack a T2D/non-diabetic label), but excluding 100% of samples confirms this dataset is the wrong one for the; Recommendation: do not run this contrast on this dataset; source pancreatic-islet samples that actually carry donor T2D vs non-diabetic status.
  - QC: 26 samples, median lib size 56,412,368, median zero-fraction 0.536.
  - Sample selection (dropped): case=0, control=0, excluded=26 — judge ⚠️ disagreed (Case set is empty (0 samples) and control set is empty (0 samples); a differential-expression contrast requires both arms to be non-empty.; All 26 samples are MYT1/MYT3 shRNA-knockdown vs Sham experiments — a gene-perturbation design, not donor disease status. The dataset contains no T2D or non-diabetic donor arm, so the requested contras)
- **GSE332956** _[primary: GEO]_ — SKIPPED: single-cell RNA-seq study — set aside from bulk DE (single-cell analysis is planned)
- **GSE301598** _[primary: GEO]_ — SKIPPED: single-cell RNA-seq study — set aside from bulk DE (single-cell analysis is planned)

## Projection
- **Power / feasibility:** Planned N per group not specified; specify it for a power read. As a rule of thumb, bulk RNA-seq needs ≈6–12 per group to detect moderate effects reliably.
- **Overall feasibility:** **cautionary**.
- **Expected QC pitfalls:**
  - Check for rRNA / globin carryover, especially in blood/plasma inputs.
  - Verify library size and 3' bias; low-input EV/exosome RNA often yields shallow, high-duplication libraries.
  - Normalization choice (median-of-ratios vs TMM vs CPM) changes DE calls — fix one and report it.
  - Batch/lane confounding with condition is the most common false-signal source; block or randomize.
  - Observed in GSE310965: low median library size (<1M counts) — typical of low-input / EV / exosome preps
- **Caveats:**
  - No prior study yielded raw counts for DE (1 had no usable per-study counts table; 3 had no assignable case/control grouping); DE-count and effect-size priors are not estimated. Try raising --reanalyze-top to reach studies deeper in the ranking that may deposit raw counts.

<sub>schema 0.1.0 · real data only · projections are cited to the datasets above and are not a substitute for formal power analysis.</sub>