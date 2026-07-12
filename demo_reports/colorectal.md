# Omics Oracle — experiment sanity check

> _RNA-seq of human colorectal cancer tumor vs adjacent normal tissue, ~10 per group, to find differentially expressed genes_

**Parsed spec** (`llm` parser): organism=**Homo sapiens**, tissue=**colorectal tissue**, omics=**rna-seq**, contrast=**colorectal cancer tumor vs adjacent normal tissue**, platform=**—**, N/group=**10**.

**Warnings:** LLM backend: claude-cli (spec parsing + narrative phrasing only; numbers come from data, not the model).; GSE313035: supplementary table is normalized (GSE313035_TPM_ctrl.txt.gz); DE will use a Welch t-test on log values (not DESeq2, which needs raw counts); GSE312260: single-cell RNA-seq (10x matrix / .h5 / .rds in RAW.tar) — set aside; bulk DE not applicable (single-cell analysis is planned); GSE285319: supplementary table is normalized (GSE285319_TPM.xlsx); DE will use a Welch t-test on log values (not DESeq2, which needs raw counts)

## Narrative
Public data does offer partial support for this plan: several human colorectal RNA-seq studies of tumor versus adjacent normal tissue already exist, the closest by match being GSE261033, GSE313035, GSE333496, and GSE245885, with GSE245885 sitting right at your intended ~10-per-group scale. However, the reanalyses came with important caveats. The datasets that were actually reprocessed (GSE313035 with 36 tumor vs 31 normal, and GSE285319 with 3 vs 3) provided only normalized expression tables rather than raw counts, so differential expression was estimated with a Welch t-test on log2 values instead of DESeq2, which weakens how directly those effect sizes translate to your design. QC also flagged a low median library size (under 1M counts) in GSE313035, consistent with a low-input/EV/exosome-style prep, and group assignments rested on file-arm labels rather than inspected clinical metadata. Notably, no genes came out as differentially expressed in two or more datasets under the concordance criteria, and the projection expects roughly zero significant genes at the stricter threshold. On the positive side, prior studies used as few as 3 per group, so N=10 is a reasonable size for moderate effects if your true effect sizes are comparable. Overall the feasibility read is cautionary rather than a clear green light: the biological question is well-precedented, but the available reanalyzed evidence is thin and method-limited, so I would confirm sample grouping and library quality before drawing strong conclusions.

## Matched datasets (ranked)
| Rank | Accession | Sim | Tier | Organism | Omics | N | Blocking? | Judge | Mismatches |
|---|---|---|---|---|---|---|---|---|---|
| 1 | [GSE261033](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261033) | 1.0 | none | Homo sapiens | rna-seq | 28 | — | irrelevant | — |
| 2 | [GSE313035](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313035) | 1.0 | proxy | Homo sapiens | rna-seq | 67 | — | partial | — |
| 3 | [GSE333496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333496) | 0.95 | none | Homo sapiens | rna-seq | 11 | — | irrelevant | n(minor) |
| 4 | [GSE245885](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245885) | 0.95 | none | Homo sapiens | rna-seq | 10 | — | irrelevant | n(minor) |
| 5 | [GSE312260](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312260) | 0.95 | proxy | Homo sapiens | rna-seq | 2 | — | — | n(minor) |
| 6 | [GSE285319](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285319) | 0.95 | proxy | Homo sapiens | rna-seq | 6 | — | — | n(minor) |

**Relevance verification (LLM judge):**
- GSE261033: **irrelevant** (conf 0.93) — The dataset profiles human Tr1 regulatory CD4+ T cells (scMultiome: transcriptomic/epigenomic/TCR) to study transcription factors driving Tr1 differentiation. Colorectal cancer appears only in passing — as one of three tumor types where a Tr1 signature is checked in tumor-resident CD4+ T cells. There is no colorectal tumor-vs-adjacent-normal tissue contrast, no bulk tissue DEG design; the biologic
- GSE313035: **partial** (conf 0.6) — Same organism (human), disease (colorectal cancer), and assay (RNA-seq), but the contrast axis differs fundamentally: this is a metastatic CRC (mCRC) clinical trial comparing treatment arms (FLOX vs FLOX+nivolumab) for immune modulation, not tumor vs adjacent normal tissue. Metadata describes no adjacent-normal samples, so the planned tumor-vs-normal DE contrast cannot be constructed; only shared 
- GSE333496: **irrelevant** (conf 0.95) — Dataset profiles human Tr1 regulatory CD4+ T cells (differentiation/transcription factors), an immunology system with no tumor-vs-adjacent-normal contrast. Colorectal cancer is mentioned only in passing as one setting where a Tr1 signature can be detected in scRNA-seq CD4+ T cells; the biological system, contrast axis, and sample composition do not match colorectal tumor vs adjacent normal tissue.
- GSE245885: **irrelevant** (conf 0.95) — Same organism (human), tissue (colon), and assay (RNA-seq), but the biological system and contrast differ entirely: the dataset profiles slow transit constipation vs non-STC colon tissue with no cancer component, whereas the plan needs colorectal tumor vs adjacent-normal. No cancer contrast is present, so it cannot serve as prior evidence for a tumor-vs-normal DE experiment.

## Per-dataset reanalysis
- **GSE313035** _[primary: GEO; welch-log2]_ — OK: DE via Welch t-test on log2 normalized values (case vs control) — this table is normalized, not raw counts, so effect sizes are log2 ratios of normalized expression and significance uses a t-test, not DESeq2
  - QC: 67 samples, median lib size 680,348, median zero-fraction 0.089. Flags: low median library size (<1M counts) — typical of low-input / EV / exosome preps
  - Sample selection (used): case=36, control=31, excluded=0 — judge ✅ agreed
  - DE: 0 genes at padj<0.05 & |log2FC|>=1 (case=36, control=31).
- **GSE312260** _[primary: GEO]_ — SKIPPED: single-cell RNA-seq study — set aside from bulk DE (single-cell analysis is planned)
- **GSE285319** _[primary: GEO; welch-log2]_ — OK: DE via Welch t-test on log2 normalized values (case vs control) — this table is normalized, not raw counts, so effect sizes are log2 ratios of normalized expression and significance uses a t-test, not DESeq2
  - QC: 6 samples, median lib size 1,000,000, median zero-fraction 0.598.
  - Sample selection (used): case=3, control=3, excluded=0 — judge ✅ agreed
  - DE: 0 genes at padj<0.05 & |log2FC|>=1 (case=3, control=3).

## Cross-experiment concordance
0 gene(s) differentially expressed in >= 2 datasets (padj<0.1, |log2fc|>=0.5), compared in HGNC-symbol space.

## Cross-dataset comparison (Venn + gene search)
Reanalyzed studies: GSE313035 [primary] (0 DE), GSE285319 [primary] (0 DE)

**Pairwise DE-gene overlap:**
- GSE313035 ∩ GSE285319: **0** shared (Jaccard 0.0; 0/0 unique)

**Your gene list (6 genes) overlap:** GSE313035 0 shared (J=0.0); GSE285319 0 shared (J=0.0)

## Concordance heatmap (gene × dataset log2FC)
Matrix: 6 genes (user_genes) × 2 datasets (GSE313035[proxy] × GSE285319[proxy]). Full data in JSON `heatmap` (UI toggles Top-N / Up-Down over `heatmap.genes`).

## Pathway enrichment (KEGG_2021_Human)

## Projection
- **Expected #DE genes:** ~0–0 (from reanalyzed priors).
- **Effect-size prior:** no significant genes in reanalyzed studies at padj<0.05, |log2fc|>=1.
- **Power / feasibility:** Planned N = 10 per group. Reasonable for moderate effects in bulk RNA-seq. Comparable prior studies used as few as 3 per group, so a detectable signal at N=10 is plausible if effect sizes match.
- **Overall feasibility:** **cautionary**.
- **Expected QC pitfalls:**
  - Check for rRNA / globin carryover, especially in blood/plasma inputs.
  - Verify library size and 3' bias; low-input EV/exosome RNA often yields shallow, high-duplication libraries.
  - Normalization choice (median-of-ratios vs TMM vs CPM) changes DE calls — fix one and report it.
  - Batch/lane confounding with condition is the most common false-signal source; block or randomize.
  - Observed in GSE313035: low median library size (<1M counts) — typical of low-input / EV / exosome preps

<sub>schema 0.1.0 · real data only · projections are cited to the datasets above and are not a substitute for formal power analysis.</sub>