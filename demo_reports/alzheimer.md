# Omics Oracle — experiment sanity check

> _RNA-seq of human Alzheimer's disease prefrontal cortex vs cognitively normal controls, post-mortem brain_

**Parsed spec** (`llm` parser): organism=**Homo sapiens**, tissue=**prefrontal cortex (post-mortem brain)**, omics=**rna-seq**, contrast=**Alzheimer's disease vs cognitively normal controls**, platform=**—**, N/group=**—**.

**Warnings:** LLM backend: claude-cli (spec parsing + narrative phrasing only; numbers come from data, not the model).; GSE310554: GSE310554_raw_counts_matrix.tsv.gz is 1093 MB > cap 52 MB; skipped; GSE310554: found candidate table(s) but none parsed; GSE291605: no supplementary expression table found (files: GSE291605_RAW.tar, filelist.txt); GSE297551: no supplementary expression table found (files: GSE297551_RAW.tar, GSE297551_WB_LJ_A10-WT_S29_Aligned.sortedByCoord.rm_dup.bed.gz, GSE297551_WB_LJ_A11-WT_S30_Aligned.sortedByCoord.rm_dup.bed.gz, GSE297551_WB_LJ_A12-WT_S31_Aligned.sortedByCoord.rm_dup.bed.gz, GSE297551_WB_LJ_A13-WT_S32_Aligned.sortedByCoord.rm_dup.bed.gz, GSE297551_WB_LJ_A15-WT_S34_Aligned.sortedByCoord.rm_dup.bed.gz); GSE296943: no supplementary expression table found (files: GSE296943_Control1_filtered_feature_bc_matrix.h5, GSE296943_Control2_filtered_feature_bc_matrix.h5, GSE296943_Hand1_filtered_feature_bc_matrix.h5, GSE296943_Hand2_filtered_feature_bc_matrix.h5, GSE296943_RAW.tar, filelist.txt); GSE212630: no supplementary expression table found (files: GSE212630_RAW.tar, filelist.txt)

## Narrative
Prior public data offers only weak support for this plan: several closely matching human RNA-seq datasets of Alzheimer's disease brain tissue exist, but none could be turned into a usable case-vs-control comparison, so the overall feasibility read is cautionary. The closest datasets by accession are GSE313307 (n=202), GSE310554 (n=39), and GSE291605 (n=31), along with GSE282552 (n=18) and GSE282551 (n=30), all human RNA-seq. The main caveats are concrete: six of the reanalyzed studies had no usable per-study counts table, and two others could not be split into non-empty Alzheimer's and cognitively-normal arms — in GSE313307 all 202 samples were dropped for lacking per-sample diagnosis metadata, and GSE275676 turned out to be a major depressive disorder study with no Alzheimer's samples at all. Because no dataset yielded raw counts for a real contrast, no differential-expression count or effect-size prior could be estimated. The planned number of samples per group was also not specified, so no power assessment is possible; as a rough guide, bulk RNA-seq typically needs about 6–12 samples per group to detect moderate effects. Watch for the usual pitfalls flagged here — normalization choice affecting DE calls, and batch or lane effects confounded with condition, which is the most common source of false signal. In short, the biological match is strong on paper, but until a study with recoverable raw counts and assignable diagnosis groups is reached, this sanity-check cannot confirm the plan is supported by existing data.

## Matched datasets (ranked)
| Rank | Accession | Sim | Tier | Organism | Omics | N | Blocking? | Judge | Mismatches |
|---|---|---|---|---|---|---|---|---|---|
| 1 | [GSE313307](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313307) | 1.0 | direct | Homo sapiens | rna-seq | 202 | — | relevant | — |
| 2 | [GSE282552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282552) | 1.0 | none | Homo sapiens | rna-seq | 18 | — | irrelevant | — |
| 3 | [GSE282551](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282551) | 1.0 | none | Homo sapiens | rna-seq | 30 | — | irrelevant | — |
| 4 | [GSE310554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310554) | 1.0 | direct | Homo sapiens | rna-seq | 39 | — | relevant | — |
| 5 | [GSE291605](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291605) | 1.0 | proxy | Homo sapiens | rna-seq | 31 | — | — | — |
| 6 | [GSE297551](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE297551) | 1.0 | proxy | Homo sapiens | rna-seq | 60 | — | — | — |
| 7 | [GSE296943](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296943) | 1.0 | proxy | Homo sapiens | rna-seq | 8 | — | — | — |
| 8 | [GSE247416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247416) | 1.0 | proxy | Homo sapiens | rna-seq | 73 | — | — | — |
| 9 | [GSE212630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212630) | 1.0 | proxy | Homo sapiens | rna-seq | 52 | — | — | — |
| 10 | [GSE275676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275676) | 0.8 | proxy | Homo sapiens | rna-seq | 78 | — | — | contrast(major) |

**Relevance verification (LLM judge):**
- GSE282552: **irrelevant** (conf 0.9) — Right organism, tissue (post-mortem DLPFC/prefrontal cortex), and short-read RNA-seq assay, but the disease and contrast axis are entirely different: this study is schizophrenia cases vs neurotypical controls, with no mention of Alzheimer's disease. It cannot serve as prior evidence for an AD-vs-control contrast despite the anatomical/assay overlap.
- GSE282551: **irrelevant** (conf 0.82) — Organism (human), tissue (postmortem DLPFC/prefrontal cortex), and assay (RNA-seq) all match, but the disease contrast is schizophrenia vs neurotypical controls — not Alzheimer's vs cognitively normal. The dataset never involves AD, so it cannot serve as prior evidence for an AD-specific contrast despite the shared tissue and assay.

**Mismatch caveats:**
- contrast terms 'Alzheimer's disease cognitively normal controls' not evident in GSE275676 metadata

## Per-dataset reanalysis
- **GSE313307** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case group (Alzheimer's disease) is empty (0 samples) — a DE contrast requires a non-empty case arm.; Control group (cognitively normal) is empty (0 samples) — a DE contrast requires a non-empty control arm.; All 202 samples were dumped into 'exclude' for identical boilerplate reason ('no diagnosis metadata; ambiguous arm'), which is implausible for a whole brain-bank cohort and suggests the diagnosis meta; If per-sample diagnosis truly cannot be recovered for any sample, the correct action is to report that the contrast cannot be built and fix the metadata linkage (e.g., match BBN accessions to the brai; The result is a null selection: nothing to compare, so the proposal cannot represent the requested AD vs cognitively-normal prefrontal-cortex contrast.
  - QC: 202 samples, median lib size 18,771,214, median zero-fraction 0.554.
  - Sample selection (dropped): case=0, control=0, excluded=202 — judge ⚠️ disagreed (Case group (Alzheimer's disease) is empty (0 samples) — a DE contrast requires a non-empty case arm.; Control group (cognitively normal) is empty (0 samples) — a DE contrast requires a non-empty control arm.)
- **GSE310554** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE291605** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE297551** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE296943** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE247416** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE212630** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE275676** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case (0) and control (0) are both empty; a differential-expression contrast requires non-empty case AND control arms.; The dataset (GSE275676) is a major depressive disorder study containing only MDD cases and MDD-study controls — there are zero Alzheimer's disease samples, so the AD-vs-control contrast simply does no; The MDD-study controls are NOT cognitively-normal AD-cohort controls; using them as the requested control arm would be off-contrast, so they cannot rescue the comparison either.; The per-sample exclusion reasons are individually correct (right call to drop MDD cases and non-AD controls), but the net result is a null selection that cannot answer the research question — the corr
  - QC: 39 samples, median lib size 19,730,705, median zero-fraction 0.428.
  - Sample selection (dropped): case=0, control=0, excluded=39 — judge ⚠️ disagreed (Case (0) and control (0) are both empty; a differential-expression contrast requires non-empty case AND control arms.; The dataset (GSE275676) is a major depressive disorder study containing only MDD cases and MDD-study controls — there are zero Alzheimer's disease samples, so the AD-vs-control contrast simply does no)

## Projection
- **Power / feasibility:** Planned N per group not specified; specify it for a power read. As a rule of thumb, bulk RNA-seq needs ≈6–12 per group to detect moderate effects reliably.
- **Overall feasibility:** **cautionary**.
- **Expected QC pitfalls:**
  - Check for rRNA / globin carryover, especially in blood/plasma inputs.
  - Verify library size and 3' bias; low-input EV/exosome RNA often yields shallow, high-duplication libraries.
  - Normalization choice (median-of-ratios vs TMM vs CPM) changes DE calls — fix one and report it.
  - Batch/lane confounding with condition is the most common false-signal source; block or randomize.
- **Caveats:**
  - No prior study yielded raw counts for DE (6 had no usable per-study counts table; 2 had no assignable case/control grouping); DE-count and effect-size priors are not estimated. Try raising --reanalyze-top to reach studies deeper in the ranking that may deposit raw counts.

<sub>schema 0.1.0 · real data only · projections are cited to the datasets above and are not a substitute for formal power analysis.</sub>