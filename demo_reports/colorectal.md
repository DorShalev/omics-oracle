# Omics Oracle — experiment sanity check

> _RNA-seq of human colorectal cancer tumor vs adjacent normal tissue, ~10 per group, to find differentially expressed genes_

**Parsed spec** (`llm` parser): organism=**Homo sapiens**, tissue=**colorectal (tumor and adjacent normal)**, omics=**rna-seq**, contrast=**colorectal cancer tumor vs adjacent normal tissue**, platform=**—**, N/group=**10**.

**Warnings:** LLM backend: claude-cli (spec parsing + narrative phrasing only; numbers come from data, not the model).; GSE312260: no supplementary expression table found (files: GSE312260_RAW.tar, filelist.txt); GSE335811: no supplementary expression table found (files: GSE335811_RAW.tar, filelist.txt); GSE324345: no supplementary expression table found (files: GSE324345_RAW.tar, filelist.txt)

## Narrative
Prior public data offers only partial, unconfirmed support for this plan: several human RNA-seq datasets match the colorectal tumor-versus-adjacent-normal design by metadata similarity, but none could be reanalyzed to actually confirm differential-expression yields. The closest datasets by accession are GSE261033 and GSE319664 (highest metadata similarity) and GSE245885, a colon-specimen expression-profiling study whose ~10-sample size mirrors your planned design; GSE333496 is also a near match. The main caveat is that no candidate study produced a usable per-study counts table — five were skipped for lacking raw counts, and one (GSE289055) was dropped because it contained only cultured cell lines from a knockout experiment rather than patient tumor/normal tissue, so no expected DE-gene count or effect-size prior could be estimated. Note too that one top-similarity hit (GSE319664) is a xenograft/drug-treatment study, which does not match a straightforward tumor-versus-normal tissue contrast. On the general RNA-seq caveats, watch for batch/lane confounding with condition, and fix and report a single normalization choice, since these are common sources of false signal. Your planned N of 10 per group is reasonable for detecting moderate effects in bulk RNA-seq. Overall feasibility is cautionary: the design is plausible and metadata-similar studies exist, but because no prior dataset could be reanalyzed into usable DE estimates, this should be treated as unvalidated — consider searching deeper in the ranking for studies that deposit raw counts before relying on any quantitative prior.

## Matched datasets (ranked)
| Rank | Accession | Sim | Tier | Organism | Omics | N | Blocking? | Judge | Mismatches |
|---|---|---|---|---|---|---|---|---|---|
| 1 | [GSE261033](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261033) | 1.0 | none | Homo sapiens | rna-seq | 28 | — | irrelevant | — |
| 2 | [GSE319664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE319664) | 1.0 | none | Homo sapiens | rna-seq | 29 | — | irrelevant | — |
| 3 | [GSE333496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333496) | 0.95 | none | Homo sapiens | rna-seq | 11 | — | irrelevant | n(minor) |
| 4 | [GSE245885](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245885) | 0.95 | none | Homo sapiens | rna-seq | 10 | — | irrelevant | n(minor) |
| 5 | [GSE312260](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312260) | 0.95 | proxy | Homo sapiens | rna-seq | 2 | — | — | n(minor) |
| 6 | [GSE289055](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289055) | 0.95 | proxy | Homo sapiens | rna-seq | 18 | — | — | n(minor) |
| 7 | [GSE285319](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285319) | 0.95 | proxy | Homo sapiens | rna-seq | 6 | — | — | n(minor) |
| 8 | [GSE335811](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335811) | 0.95 | proxy | Homo sapiens | rna-seq | 12 | — | — | n(minor) |
| 9 | [GSE337988](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337988) | 0.8 | proxy | Homo sapiens | rna-seq | 48 | — | — | contrast(major) |
| 10 | [GSE324345](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324345) | 0.8 | proxy | Homo sapiens | rna-seq | 66 | — | — | contrast(major) |

**Relevance verification (LLM judge):**
- GSE261033: **irrelevant** (conf 0.94) — This dataset studies human Tr1 regulatory CD4+ T-cell differentiation (scMultiome/epigenomic profiling of antigen-induced T cells, focus on IRF4/BATF/MAF). Colorectal cancer is mentioned only in passing as one of several tumor types where a Tr1 signature was detected in tumor-resident CD4+ T cells. There is no tumor-vs-adjacent-normal colorectal tissue contrast, no bulk colorectal tumor/normal RNA
- GSE319664: **irrelevant** (conf 0.9) — Both are human RNA-seq, but the dataset profiles tumors from a xenograft model treated with KDM5B inhibitors (JB-157/JB-161) — a drug-treatment vs vehicle contrast in a xenograft, not colorectal tumor vs adjacent normal tissue. No indication the tumors are colorectal or that adjacent-normal tissue is included; the biological system and contrast axis do not match the plan.
- GSE333496: **irrelevant** (conf 0.95) — Despite RNA-seq/Homo sapiens/keyword match, the actual system is human Tr1 regulatory CD4+ T cells and IRF4/BATF/MAF-driven Tr1 differentiation. Colorectal cancer is mentioned only in passing as one tissue where a Tr1 signature was projected onto external scRNA-seq. There is no CRC tumor-vs-adjacent-normal contrast, no colorectal tissue samples, and a different biological system/assay focus, so it
- GSE245885: **irrelevant** (conf 0.9) — Same organism (human), same assay (RNA-seq), and same organ region (colon), but the biological system and contrast are entirely different: the dataset compares slow transit constipation (STC) vs non-STC colon tissue, with no cancer/tumor component. There is no tumor-vs-adjacent-normal contrast, so it cannot serve as prior evidence for colorectal cancer DEGs. Overlap is limited to tissue keyword on

**Mismatch caveats:**
- contrast terms 'colorectal cancer tumor adjacent normal tissue' not evident in GSE337988 metadata
- contrast terms 'colorectal cancer tumor adjacent normal tissue' not evident in GSE324345 metadata

## Per-dataset reanalysis
- **GSE312260** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE289055** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case set is empty (0 samples) — no colorectal cancer tumor tissue samples exist here.; Control set is empty (0 samples) — no adjacent normal tissue samples exist here.; All 18 samples are cultured cell lines from a gene-knockout experiment (BxPC3 pancreatic, SW480 CRC, KBM7 leukemia), not patient tissue; the tumor-vs-adjacent-normal design does not exist in this data; Individual exclusion reasons are correct, but the net result is a non-runnable DE contrast — this dataset cannot answer the stated research question and should be flagged as no valid samples rather th
  - QC: 18 samples, median lib size 43,187,570, median zero-fraction 0.639.
  - Sample selection (dropped): case=0, control=0, excluded=18 — judge ⚠️ disagreed (Case set is empty (0 samples) — no colorectal cancer tumor tissue samples exist here.; Control set is empty (0 samples) — no adjacent normal tissue samples exist here.)
- **GSE285319** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE335811** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE337988** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE324345** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset

## Projection
- **Power / feasibility:** Planned N = 10 per group. Reasonable for moderate effects in bulk RNA-seq.
- **Overall feasibility:** **cautionary**.
- **Expected QC pitfalls:**
  - Check for rRNA / globin carryover, especially in blood/plasma inputs.
  - Verify library size and 3' bias; low-input EV/exosome RNA often yields shallow, high-duplication libraries.
  - Normalization choice (median-of-ratios vs TMM vs CPM) changes DE calls — fix one and report it.
  - Batch/lane confounding with condition is the most common false-signal source; block or randomize.
- **Caveats:**
  - No prior study yielded raw counts for DE (5 had no usable per-study counts table; 1 had no assignable case/control grouping); DE-count and effect-size priors are not estimated. Try raising --reanalyze-top to reach studies deeper in the ranking that may deposit raw counts.

<sub>schema 0.1.0 · real data only · projections are cited to the datasets above and are not a substitute for formal power analysis.</sub>