# Omics Oracle — experiment sanity check

> _RNA-seq of human Alzheimer's disease prefrontal cortex vs cognitively normal controls, post-mortem brain_

**Parsed spec** (`llm` parser): organism=**Homo sapiens**, tissue=**prefrontal cortex (post-mortem brain)**, omics=**rna-seq**, contrast=**Alzheimer's disease vs cognitively normal controls**, platform=**—**, N/group=**—**.

**Warnings:** LLM backend: claude-cli (spec parsing + narrative phrasing only; numbers come from data, not the model).; GSE310554: GSE310554_raw_counts_matrix.tsv.gz is 1093 MB > cap 52 MB; skipped; GSE310554: found candidate table(s) but none parsed; GSE291605: GSE291605_RAW.tar is 962 MB > cap 419 MB; skipped; GSE291605: no supplementary expression table found (files: GSE291605_RAW.tar, filelist.txt); GSE297551: GSE297551_RAW.tar had no per-sample count tables to assemble; GSE297551: no supplementary expression table found (files: GSE297551_RAW.tar, GSE297551_WB_LJ_A10-WT_S29_Aligned.sortedByCoord.rm_dup.bed.gz, GSE297551_WB_LJ_A11-WT_S30_Aligned.sortedByCoord.rm_dup.bed.gz, GSE297551_WB_LJ_A12-WT_S31_Aligned.sortedByCoord.rm_dup.bed.gz, GSE297551_WB_LJ_A13-WT_S32_Aligned.sortedByCoord.rm_dup.bed.gz, GSE297551_WB_LJ_A15-WT_S34_Aligned.sortedByCoord.rm_dup.bed.gz); GSE296943: GSE296943_RAW.tar is 5290 MB > cap 419 MB; skipped; GSE296943: single-cell RNA-seq deposition (.mtx/.h5/.rds) — set aside; bulk DE not applicable (single-cell analysis is planned); GSE247416: could not parse GSE247416_old_new_numberings.xlsx (no numeric sample columns); GSE247416: found candidate table(s) but none parsed; GSE212630: single-cell RNA-seq (10x matrix / .h5 / .rds in RAW.tar) — set aside; bulk DE not applicable (single-cell analysis is planned)

## Narrative
**Sanity-check summary — human AD prefrontal cortex RNA-seq vs. cognitively normal controls**

Public data broadly supports the premise: several human, post-mortem brain RNA-seq studies matched the plan closely, with the nearest being GSE313307 (n=202), GSE282551 (n=30), GSE310554 (n=39), GSE291605 (n=31), and GSE282552 (n=18). However, none of these could be turned into a working Alzheimer's-vs-control differential-expression check — four of the ranked studies (including GSE310554 and GSE291605) had no usable per-study counts table, and GSE313307, despite passing QC on all 202 samples, had no case/control grouping that could be established. Two matches (GSE296943, GSE212630) are single-cell studies set aside from bulk DE, and one candidate, GSE275676, turned out to be a major-depressive-disorder cohort with no Alzheimer's samples at all, so its 39 samples were correctly not repurposed. Because no prior study yielded raw counts, expected DE-gene counts and effect-size priors could not be estimated here. The plan also did not specify N per group; as a rough guide, bulk RNA-seq typically needs about 6–12 per group to detect moderate effects, and batch/condition confounding plus normalization choice are the usual sources of false signal to guard against. Overall the feasibility read is **cautionary**: the disease, tissue, organism, and assay are all well-represented in public data, but this particular check could not be quantitatively grounded, so consider reaching deeper into the ranking for studies that deposit raw counts before drawing conclusions.

## Matched datasets (ranked)
| Rank | Accession | Sim | Tier | Organism | Omics | N | Blocking? | Judge | Mismatches |
|---|---|---|---|---|---|---|---|---|---|
| 1 | [GSE313307](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313307) | 1.0 | direct | Homo sapiens | rna-seq | 202 | — | relevant | — |
| 2 | [GSE282552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282552) | 1.0 | none | Homo sapiens | rna-seq | 18 | — | irrelevant | — |
| 3 | [GSE282551](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282551) | 1.0 | none | Homo sapiens | rna-seq | 30 | — | irrelevant | — |
| 4 | [GSE310554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310554) | 1.0 | proxy | Homo sapiens | rna-seq | 39 | — | partial | — |
| 5 | [GSE291605](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291605) | 1.0 | proxy | Homo sapiens | rna-seq | 31 | — | — | — |
| 6 | [GSE297551](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297551) | 1.0 | proxy | Homo sapiens | rna-seq | 60 | — | — | — |
| 7 | [GSE296943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296943) | 1.0 | proxy | Homo sapiens | rna-seq | 8 | — | — | — |
| 8 | [GSE247416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247416) | 1.0 | proxy | Homo sapiens | rna-seq | 73 | — | — | — |
| 9 | [GSE212630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212630) | 1.0 | proxy | Homo sapiens | rna-seq | 52 | — | — | — |
| 10 | [GSE275676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275676) | 0.8 | proxy | Homo sapiens | rna-seq | 78 | — | — | contrast(major) |

**Relevance verification (LLM judge):**
- GSE282552: **irrelevant** (conf 0.9) — Matches organism (human), tissue (post-mortem DLPFC/prefrontal cortex), and assay (short-read RNA-seq), but the biological system and contrast are schizophrenia cases vs neurotypical controls with an isoform-usage focus — Alzheimer's disease is never mentioned. The disease and case-control axis differ entirely, so it cannot serve as prior evidence for an AD-vs-control study despite tissue/assay ov
- GSE282551: **irrelevant** (conf 0.82) — Same organism (human), same tissue (post-mortem DLPFC/prefrontal cortex), and same broad assay (RNA-seq, though long-read isoform-focused) — but the disease contrast is schizophrenia vs neurotypical controls, not Alzheimer's disease vs cognitively normal. The planned AD case-control axis is entirely absent; a schizophrenia isoform-usage study cannot serve as prior evidence for AD-associated expres
- GSE310554: **partial** (conf 0.72) — Matches organism (human), tissue (prefrontal cortex, post-mortem), and contrast (AD vs cognitively normal controls). However the assay is single-nuclei RNA-seq (snRNA-seq) focused on the neurovascular unit/vascular cell types, not the bulk prefrontal-cortex RNA-seq the plan implies; usable as prior evidence but not a direct assay match.

**Mismatch caveats:**
- contrast terms 'Alzheimer's disease cognitively normal controls' not evident in GSE275676 metadata

## Per-dataset reanalysis
- **GSE313307** _[primary: GEO]_ — SKIPPED: no case/control selection could be established for this study; QC reported, DE skipped (no contrast invented)
  - QC: 202 samples, median lib size 18,771,214, median zero-fraction 0.554.
- **GSE310554** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE291605** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE297551** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE296943** _[primary: GEO]_ — SKIPPED: single-cell RNA-seq study — set aside from bulk DE (single-cell analysis is planned)
- **GSE247416** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE212630** _[primary: GEO]_ — SKIPPED: single-cell RNA-seq study — set aside from bulk DE (single-cell analysis is planned)
- **GSE275676** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case set is empty (0 samples) — no Alzheimer's disease samples exist in this dataset.; Control set is empty (0 samples) — the 20 'ctrl' samples are MDD-study controls, not cognitively-normal AD controls, and were (reasonably) not repurposed.; The dataset is a major depressive disorder cohort (mdd cases + study-matched ctrls); the AD-vs-normal contrast does not exist in these 39 samples.; With both arms empty, no pydeseq2/limma contrast can be fit; the correct action is to reject the analysis for this dataset, not to proceed with an empty selection.
  - QC: 39 samples, median lib size 19,730,705, median zero-fraction 0.428.
  - Sample selection (dropped): case=0, control=0, excluded=39 — judge ⚠️ disagreed (Case set is empty (0 samples) — no Alzheimer's disease samples exist in this dataset.; Control set is empty (0 samples) — the 20 'ctrl' samples are MDD-study controls, not cognitively-normal AD controls, and were (reasonably) not repurposed.)

## Projection
- **Power / feasibility:** Planned N per group not specified; specify it for a power read. As a rule of thumb, bulk RNA-seq needs ≈6–12 per group to detect moderate effects reliably.
- **Overall feasibility:** **cautionary**.
- **Expected QC pitfalls:**
  - Check for rRNA / globin carryover, especially in blood/plasma inputs.
  - Verify library size and 3' bias; low-input EV/exosome RNA often yields shallow, high-duplication libraries.
  - Normalization choice (median-of-ratios vs TMM vs CPM) changes DE calls — fix one and report it.
  - Batch/lane confounding with condition is the most common false-signal source; block or randomize.
- **Caveats:**
  - No prior study yielded raw counts for DE (4 had no usable per-study counts table; 1 had no assignable case/control grouping); DE-count and effect-size priors are not estimated. Try raising --reanalyze-top to reach studies deeper in the ranking that may deposit raw counts.

<sub>schema 0.1.0 · real data only · projections are cited to the datasets above and are not a substitute for formal power analysis.</sub>