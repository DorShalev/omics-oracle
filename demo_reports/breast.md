# Omics Oracle — experiment sanity check

> _Human breast cancer tumor vs normal breast tissue RNA-seq, triple-negative, ~8 per group_

**Parsed spec** (`llm` parser): organism=**Homo sapiens**, tissue=**breast**, omics=**rna-seq**, contrast=**triple-negative breast cancer tumor vs normal breast tissue**, platform=**—**, N/group=**8**.

**Warnings:** LLM backend: claude-cli (spec parsing + narrative phrasing only; numbers come from data, not the model).; GSE273078: supplementary table is normalized (GSE273078_DESeq2_output_26to30.csv.gz); DE requires raw counts and will be skipped (not faked); GSE309543: no supplementary expression table found (files: GSE309543_RAW.tar, filelist.txt); GSE303509: no supplementary expression table found (files: GSE303509_RAW.tar, filelist.txt); GSE312563: no supplementary expression table found (files: GSE312563_Hs578_T_vehicle_vs_Hs578_T_24h_6403_DEG_by_DESeq2.csv.gz, GSE312563_MDA_231_vehicle_vs_MDA_231_24h_6403_DEG_by_DESeq2.csv.gz, GSE312563_MDA_436_vehicle_vs_MDA_436_24h_6403_DEG_by_DESeq2.csv.gz, GSE312563_RAW.tar, filelist.txt); GSE334884: supplementary table is normalized (GSE334884_PSI_MDA-MB231_MCF10a_Stiffness.csv.gz); DE requires raw counts and will be skipped (not faked); Narrative mentioned gene-like tokens not found in the analyzed data — verify: MDA, MCF

## Narrative
Prior public data offers only weak, indirect support for this plan. The keyword search surfaced several high-similarity human breast RNA-seq studies — the closest by accession being GSE261033, GSE273078, GSE337023, GSE333496, and GSE309543 — but on reanalysis none of these actually contained a triple-negative-tumor-versus-normal-breast contrast: they turned out to be shRNA-knockdown ribosome profiling, cell-line co-cultures, drug/knockout perturbations, or lacked usable per-study counts tables. The only dataset that yielded a differential-expression result was GSE334884, and even there the "tumor vs normal" arms were cell-line proxies (MDA-MB-231 vs MCF-10A), analyzed 2-vs-2 by a Welch t-test on normalized (not raw) values, so its 23 significant genes and effect-size prior (median |log2FC| ≈ 3.60) should be read with real caution. Key mismatches to flag: no true primary-tissue tumor-vs-normal cohort was recovered, replicates were effectively n=1 per condition in the one usable study, and GSE334884 carried a low-library-size QC flag typical of low-input preps. Because only one dataset produced a DE table (two are needed), no cross-experiment concordance could be computed. On the positive side, a planned N=8 per group is reasonable for detecting moderate effects in bulk RNA-seq, and comparable priors suggest signal is plausible if effect sizes hold. Overall the feasibility read is **cautionary** — the design is sound in principle, but the public evidence base here is thin and proxy-based rather than a confirmed tissue-level precedent.

## Matched datasets (ranked)
| Rank | Accession | Sim | Tier | Organism | Omics | N | Blocking? | Judge | Mismatches |
|---|---|---|---|---|---|---|---|---|---|
| 1 | [GSE261033](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261033) | 1.0 | none | Homo sapiens | rna-seq | 28 | — | irrelevant | — |
| 2 | [GSE273078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273078) | 1.0 | proxy | Homo sapiens | rna-seq | 18 | — | partial | — |
| 3 | [GSE337023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337023) | 1.0 | proxy | Homo sapiens | rna-seq | 24 | — | partial | — |
| 4 | [GSE333496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333496) | 0.95 | none | Homo sapiens | rna-seq | 11 | — | irrelevant | n(minor) |
| 5 | [GSE309543](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309543) | 0.95 | proxy | Homo sapiens | rna-seq | 4 | — | — | n(minor) |
| 6 | [GSE303509](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303509) | 0.95 | proxy | Homo sapiens | rna-seq | 14 | — | — | n(minor) |
| 7 | [GSE302168](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302168) | 0.95 | proxy | Homo sapiens | rna-seq | 12 | — | — | n(minor) |
| 8 | [GSE312563](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312563) | 0.95 | proxy | Homo sapiens | rna-seq | 12 | — | — | n(minor) |
| 9 | [GSE334884](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334884) | 0.5 | proxy | Homo sapiens | rna-seq | 4 | — | — | tissue(major); contrast(major); n(minor) |

**Relevance verification (LLM judge):**
- GSE261033: **irrelevant** (conf 0.95) — Dataset profiles human Tr1 regulatory CD4+ T cells (immune cell differentiation, TF regulators via scMultiome/single-cell), not bulk RNA-seq of breast tumor vs normal tissue. TNBC appears only in passing as one of several cancers where a Tr1 signature was checked in tumor-resident T cells; there is no tumor-vs-normal breast tissue contrast and the biological system, assay, and contrast axis all di
- GSE273078: **partial** (conf 0.72) — Shares organism (human) and disease (TNBC) with the plan, but the actual contrast axis is fibrillarin depletion vs control (translational control), not tumor vs normal breast tissue. The assay is Ribo-seq (translation efficiency) rather than standard RNA-seq DE, and there is no indication of normal breast tissue samples or a tumor-vs-normal comparison. Keyword overlap on TNBC only; cannot serve as
- GSE337023: **partial** (conf 0.6) — Same organism (human), disease (TNBC), and assay (RNA-seq/TagSeq), but the biological system is cell lines (HCC1806, MDA-MB-231), not breast tissue, and there is no normal-tissue arm. The contrast axis is whole-genome-doubling mechanism vs parental aneuploid controls, not tumor vs normal. It shares the disease context but cannot serve as tumor-vs-normal prior evidence.
- GSE333496: **irrelevant** (conf 0.94) — The dataset profiles human Tr1 (CD4+ regulatory T) cells and their transcriptional regulators (IRF4/BATF/MAF), not breast tumor vs normal tissue. TNBC is mentioned only in passing as one of several tumor types where a Tr1 signature was later applied to CD4+ T cells. Wrong biological system (immune cells, not breast epithelium/tissue), no tumor-vs-normal-breast contrast, and no per-group breast coh

**Mismatch caveats:**
- requested tissue 'breast' not evident in GSE334884 metadata
- contrast terms 'triple-negative breast cancer tumor normal breast tissue' not evident in GSE334884 metadata

## Per-dataset reanalysis
- **GSE273078** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case (n=0) and control (n=0) are both empty — a DE contrast requires two non-empty groups; nothing here can be compared.; The 49 items are all from a shRNA-knockdown ribosome-profiling study (GSE273078): condSH2/condSH3 knockdowns vs condSC scrambled control, plus IN/RPF and derived DESeq2 normalized/statistic columns — ; The individual exclusion reasons are factually correct (knockdown/control ribo-seq samples and derived stat columns are legitimately off-contrast), but the dataset simply lacks any TNBC-tumor or norma; This dataset appears to be the wrong resource for the stated research question; a TNBC-vs-normal-breast contrast needs samples from a tumor/normal cohort (e.g., TCGA-BRCA TNBC vs adjacent-normal), not
  - QC: 49 samples, median lib size 7,984,315, median zero-fraction 0.134.
  - Sample selection (dropped): case=0, control=0, excluded=49 — judge ⚠️ disagreed (Case (n=0) and control (n=0) are both empty — a DE contrast requires two non-empty groups; nothing here can be compared.; The 49 items are all from a shRNA-knockdown ribosome-profiling study (GSE273078): condSH2/condSH3 knockdowns vs condSC scrambled control, plus IN/RPF and derived DESeq2 normalized/statistic columns — )
- **GSE337023** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case set is empty (0 samples) and control set is empty (0 samples); a DE contrast requires both arms non-empty.; The dataset is entirely in-vitro (HCC1806 and MDA-MB-231 cell-line monocultures, cancer+fibroblast co-cultures, fibroblast monocultures, and 'MS' conditions) — it contains neither primary TNBC tumor t; The proposer's exclusion reasoning is biologically correct (none of the 24 samples match either requested arm), but the correct verdict is that the contrast is not answerable here, not that a valid ca; If a cell-line-based proxy contrast were acceptable (e.g., TNBC cancer monoculture vs fibroblast/normal-like control), that would need to be explicitly restated as a different research question — it d
  - QC: 24 samples, median lib size 3,458,930, median zero-fraction 0.504.
  - Sample selection (dropped): case=0, control=0, excluded=24 — judge ⚠️ disagreed (Case set is empty (0 samples) and control set is empty (0 samples); a DE contrast requires both arms non-empty.; The dataset is entirely in-vitro (HCC1806 and MDA-MB-231 cell-line monocultures, cancer+fibroblast co-cultures, fibroblast monocultures, and 'MS' conditions) — it contains neither primary TNBC tumor t)
- **GSE309543** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE303509** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE302168** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case (TNBC tumor) = 0 samples and control (normal breast tissue) = 0 samples — no contrast can be formed.; All 12 samples are cell-line perturbation conditions (GLI2-KO/WT × DMSO/PTX), i.e. a drug/genetic-perturbation experiment, not tumor vs normal tissue; the dataset simply does not contain the requested; The per-sample exclusion reasoning is individually correct (none are tumor or normal tissue), but the net result is a null selection — the right conclusion is that this contrast is unanswerable in thi; If the actual biological interest here is drug response or GLI2 knockout effect (e.g. WT_DMSO vs WT_PTX, or WT vs GLI2_KO), the research question/contrast should be redefined to match the samples rath
  - QC: 12 samples, median lib size 89,941,280, median zero-fraction 0.021.
  - Sample selection (dropped): case=0, control=0, excluded=12 — judge ⚠️ disagreed (Case (TNBC tumor) = 0 samples and control (normal breast tissue) = 0 samples — no contrast can be formed.; All 12 samples are cell-line perturbation conditions (GLI2-KO/WT × DMSO/PTX), i.e. a drug/genetic-perturbation experiment, not tumor vs normal tissue; the dataset simply does not contain the requested)
- **GSE312563** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE334884** _[primary: GEO; welch-log2]_ — OK: DE via Welch t-test on log2 normalized values (case vs control) — this table is normalized, not raw counts, so effect sizes are log2 ratios of normalized expression and significance uses a t-test, not DESeq2
  - QC: 7 samples, median lib size 603,308, median zero-fraction 0.42. Flags: low median library size (<1M counts) — typical of low-input / EV / exosome preps
  - Sample selection (used): case=2, control=2, excluded=3 — judge ✅ agreed
  - DE: 23 genes at padj<0.05 & |log2FC|>=1 (case=2, control=2).

## Cross-experiment concordance
only 1 dataset(s) yielded DE tables; need >= 2 for cross-experiment concordance

## Cross-dataset comparison (Venn + gene search)
Reanalyzed studies: GSE334884 [primary] (23 DE)

**Your gene list (6 genes) overlap:** GSE334884 0 shared (J=0.0)

## Concordance heatmap (gene × dataset log2FC)
Matrix: 6 genes (user_genes) × 1 datasets (GSE334884[proxy]). Full data in JSON `heatmap` (UI toggles Top-N / Up-Down over `heatmap.genes`).

## Pathway enrichment (KEGG_2021_Human)

## Projection
- **Expected #DE genes:** ~23–23 (from reanalyzed priors).
- **Effect-size prior:** median |log2FC| of significant genes ≈ 3.60 across reanalyzed studies.
- **Power / feasibility:** Planned N = 8 per group. Reasonable for moderate effects in bulk RNA-seq. Comparable prior studies used as few as 2 per group, so a detectable signal at N=8 is plausible if effect sizes match.
- **Overall feasibility:** **cautionary**.
- **Expected QC pitfalls:**
  - Check for rRNA / globin carryover, especially in blood/plasma inputs.
  - Verify library size and 3' bias; low-input EV/exosome RNA often yields shallow, high-duplication libraries.
  - Normalization choice (median-of-ratios vs TMM vs CPM) changes DE calls — fix one and report it.
  - Batch/lane confounding with condition is the most common false-signal source; block or randomize.
  - Observed in GSE334884: low median library size (<1M counts) — typical of low-input / EV / exosome preps

<sub>schema 0.1.0 · real data only · projections are cited to the datasets above and are not a substitute for formal power analysis.</sub>