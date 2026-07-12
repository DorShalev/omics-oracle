# Omics Oracle — experiment sanity check

> _Human pancreatic islet RNA-seq, type 2 diabetes vs non-diabetic donors_

**Parsed spec** (`llm` parser): organism=**Homo sapiens**, tissue=**pancreatic islet**, omics=**rna-seq**, contrast=**type 2 diabetes donors vs non-diabetic donors**, platform=**—**, N/group=**—**.

**Warnings:** LLM backend: claude-cli (spec parsing + narrative phrasing only; numbers come from data, not the model).; GSE285253: no supplementary expression table found (files: GSE285253_RAW.tar, filelist.txt); GSE333721: no supplementary expression table found (files: GSE333721_stressor_screen_6conditions.rds); GSE305932: no supplementary expression table found (files: GSE305932_RAW.tar, filelist.txt); GSE332956: no supplementary expression table found (files: GSE332956_RAW.tar, filelist.txt); GSE301598: no supplementary expression table found (files: GSE301598_RAW.tar, filelist.txt)

## Narrative
Prior public data offers only weak support for this plan as it stands: several human pancreatic-islet RNA-seq studies rank as close matches by keyword, but none could actually be reanalyzed for the requested type 2 diabetes vs non-diabetic contrast. The closest datasets by accession are GSE335464 (n=28), GSE312607 (n=12), GSE285253 (n=12), GSE308718 (n=6), and GSE333721 (n=6), all Homo sapiens RNA-seq. The main caveat is that six of the reanalyzed studies (including GSE285253 and GSE333721) had no usable per-study counts table, so no differential-expression numbers or effect-size priors could be estimated. A further mismatch: GSE307751's 26 samples were all shMYT1/shMYT3 knockdown-vs-Sham perturbation experiments with no donor diabetic status annotated, so the case and control arms came out empty and the contrast could not be run. Note too that no planned N per group was specified, and as a rough guide bulk RNA-seq typically needs about 6–12 samples per group to detect moderate effects. Because no prior study yielded raw counts for this contrast, the feasibility read is cautionary rather than confident. It may be worth reaching deeper into the ranking (e.g. raising the reanalysis depth) to find studies that deposit raw counts before committing to this design.

## Matched datasets (ranked)
| Rank | Accession | Sim | Tier | Organism | Omics | N | Blocking? | Judge | Mismatches |
|---|---|---|---|---|---|---|---|---|---|
| 1 | [GSE335464](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335464) | 1.0 | none | Homo sapiens | rna-seq | 28 | — | irrelevant | — |
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
- GSE335464: **irrelevant** (conf 0.82) — Organism (human), tissue (pancreatic islet), and assay (RNA-seq) match, but the dataset studies type 1 diabetes (T1D GWAS mechanisms, baseline vs T1D-stimulating conditions, and hESC-derived beta-like knockout cells), not type 2 diabetes. There is no T2D-vs-non-diabetic donor contrast; the biological question and contrast axis differ entirely, and much of the RNA-seq is on isogenic hESC-derived ce
- GSE308718: **irrelevant** (conf 0.9) — Both are human RNA-seq mentioning beta cells/diabetes, but the actual system and contrast don't match: this dataset studies stem cell-derived beta cells (SC-β) co-aggregated with macrophages in different polarization states, aimed at type 1 diabetes cell therapy — not primary pancreatic islet tissue. The contrast is macrophage polarization (SC-MReg/Inf/Unp), not T2D vs non-diabetic donors. Wrong b
- GSE312607: **irrelevant** (conf 0.72) — Same organism/assay (human RNA-seq) and islet/diabetes theme, but the actual biology and contrast diverge on every axis that matters: it studies HNF1A-MODY (monogenic diabetes), not type 2 diabetes; the contrast is mutant-vs-control genotype during development, not diabetic-vs-non-diabetic donor state; and the samples are stem-cell-derived islets, not primary donor pancreatic islet tissue. Keyword
- GSE285253: **partial** (conf 0.8) — Matches organism (human), assay (RNA-seq), and islet cell-type, but diverges on both key axes: the contrast is an isogenic HNF1A-MODY mutation (monogenic diabetes) vs control, not T2D-vs-nondiabetic donors, and the system is stem-cell-derived islets rather than primary donor pancreatic islet tissue. Related islet/diabetes biology but wrong contrast and model system, so only weak prior evidence.

**Mismatch caveats:**
- contrast terms 'type 2 diabetes donors non-diabetic donors' not evident in GSE332956 metadata
- contrast terms 'type 2 diabetes donors non-diabetic donors' not evident in GSE301598 metadata

## Per-dataset reanalysis
- **GSE285253** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE333721** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE305932** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE310965** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE307751** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case set (T2D donors) is empty and control set (non-diabetic donors) is empty; a DE contrast requires both arms to be non-empty.; The 26 samples are all gene-perturbation experiments (shMYT1/shMYT3 knockdowns and Sham knockdown controls across donors 8695, 9201-GG, 9508, 9994) — the dataset carries a knockdown-vs-Sham design, no; The T2D-vs-non-diabetic contrast is simply not present in these samples: donor diabetic status is not annotated, only perturbation arm. Excluding all 26 is defensible on its face, but it leaves nothin; Note: the Sham samples are perturbation controls, NOT 'non-diabetic donor' controls — so they correctly should not be repurposed as the control arm for this contrast; doing so would be an off-contrast
  - QC: 26 samples, median lib size 56,412,368, median zero-fraction 0.536.
  - Sample selection (dropped): case=0, control=0, excluded=26 — judge ⚠️ disagreed (Case set (T2D donors) is empty and control set (non-diabetic donors) is empty; a DE contrast requires both arms to be non-empty.; The 26 samples are all gene-perturbation experiments (shMYT1/shMYT3 knockdowns and Sham knockdown controls across donors 8695, 9201-GG, 9508, 9994) — the dataset carries a knockdown-vs-Sham design, no)
- **GSE332956** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE301598** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset

## Projection
- **Power / feasibility:** Planned N per group not specified; specify it for a power read. As a rule of thumb, bulk RNA-seq needs ≈6–12 per group to detect moderate effects reliably.
- **Overall feasibility:** **cautionary**.
- **Expected QC pitfalls:**
  - Check for rRNA / globin carryover, especially in blood/plasma inputs.
  - Verify library size and 3' bias; low-input EV/exosome RNA often yields shallow, high-duplication libraries.
  - Normalization choice (median-of-ratios vs TMM vs CPM) changes DE calls — fix one and report it.
  - Batch/lane confounding with condition is the most common false-signal source; block or randomize.
- **Caveats:**
  - No prior study yielded raw counts for DE (6 had no usable per-study counts table; 1 had no assignable case/control grouping); DE-count and effect-size priors are not estimated. Try raising --reanalyze-top to reach studies deeper in the ranking that may deposit raw counts.

<sub>schema 0.1.0 · real data only · projections are cited to the datasets above and are not a substitute for formal power analysis.</sub>