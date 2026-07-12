# Omics Oracle — experiment sanity check

> _Human breast cancer tumor vs normal breast tissue RNA-seq, triple-negative, ~8 per group_

**Parsed spec** (`llm` parser): organism=**Homo sapiens**, tissue=**breast**, omics=**rna-seq**, contrast=**triple-negative breast cancer tumor vs normal breast tissue**, platform=**—**, N/group=**8**.

**Warnings:** LLM backend: claude-cli (spec parsing + narrative phrasing only; numbers come from data, not the model).; GSE273078: supplementary table is normalized (GSE273078_DESeq2_output_26to30.csv.gz); DE will use a Welch t-test on log values (not DESeq2, which needs raw counts); GSE309543: single-cell RNA-seq (10x matrix / .h5 / .rds in RAW.tar) — set aside; bulk DE not applicable (single-cell analysis is planned); GSE303509: assembled a 44391×14 matrix from 14 per-sample files inside GSE303509_RAW.tar (raw); GSE312563: assembled a 162473×12 matrix from 12 per-sample files inside GSE312563_RAW.tar (normalized); GSE334884: supplementary table is normalized (GSE334884_PSI_MDA-MB231_MCF10a_Stiffness.csv.gz); DE will use a Welch t-test on log values (not DESeq2, which needs raw counts); Narrative mentioned gene-like tokens not found in the analyzed data — verify: MDA, TNBC, MCF10A

## Narrative
Prior public data offers only weak, indirect support for this plan: although several Homo sapiens breast RNA-seq datasets matched the query with high similarity (GSE261033, GSE273078, and GSE337023 at sim=1.0, plus GSE333496 and GSE309543 at sim=0.95), the reanalysis could not form the requested triple-negative-tumor-vs-normal-breast contrast in almost all of them. The closest dataset that actually yielded a differential-expression table was GSE334884, and even there the "case" and "control" arms are cell-line proxies (MDA-MB-231 as TNBC, MCF10A as normal breast epithelium) rather than primary tumor and normal tissue. That single usable comparison came from a substrate-stiffness experiment with only 2 samples per arm and was run on normalized values with a Welch t-test rather than DESeq2 on raw counts, so its ~23 significant genes should be read as preliminary. The other candidates were set aside for legitimate reasons — ribosome-profiling knockdowns (GSE273078), tumor-only cohorts with no normal arm (GSE303509), in-vitro knockout/drug designs (GSE302168, GSE312563), and single-cell data (GSE309543) — meaning no true tumor-vs-normal-tissue design was recovered. A further limitation is that only one dataset produced a DE table, so no cross-experiment concordance check (which needs at least two) was possible, and several studies carried low-library-size QC flags. The planned N=8 per group is reasonable for detecting moderate effects in bulk RNA-seq, and prior studies detected signal with as few as 2 per arm. Overall the feasibility read is cautionary: the experiment is plausible at the proposed sample size, but the available public evidence is thin, indirect, and cell-line-based, so I would not overstate how well it is de-risked.

## Matched datasets (ranked)
| Rank | Accession | Sim | Tier | Organism | Omics | N | Blocking? | Judge | Mismatches |
|---|---|---|---|---|---|---|---|---|---|
| 1 | [GSE261033](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261033) | 1.0 | none | Homo sapiens | rna-seq | 28 | — | irrelevant | — |
| 2 | [GSE273078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273078) | 1.0 | proxy | Homo sapiens | rna-seq | 18 | — | partial | — |
| 3 | [GSE337023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337023) | 1.0 | none | Homo sapiens | rna-seq | 24 | — | irrelevant | — |
| 4 | [GSE333496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333496) | 0.95 | none | Homo sapiens | rna-seq | 11 | — | irrelevant | n(minor) |
| 5 | [GSE309543](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309543) | 0.95 | proxy | Homo sapiens | rna-seq | 4 | — | — | n(minor) |
| 6 | [GSE303509](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303509) | 0.95 | proxy | Homo sapiens | rna-seq | 14 | — | — | n(minor) |
| 7 | [GSE302168](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302168) | 0.95 | proxy | Homo sapiens | rna-seq | 12 | — | — | n(minor) |
| 8 | [GSE312563](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312563) | 0.95 | proxy | Homo sapiens | rna-seq | 12 | — | — | n(minor) |
| 9 | [GSE334884](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334884) | 0.5 | proxy | Homo sapiens | rna-seq | 4 | — | — | tissue(major); contrast(major); n(minor) |

**Relevance verification (LLM judge):**
- GSE261033: **irrelevant** (conf 0.94) — Dataset studies human Tr1 (CD4+ regulatory T) cell differentiation and its transcriptional regulators via scMultiome/single-cell profiling; TNBC is mentioned only as one tissue context where a Tr1 signature is checked among tumor-resident CD4+ T cells. The biological system (immune T cells, not breast tissue), the contrast axis (Tr1 differentiation/TF regulation, not tumor vs normal), and the assa
- GSE273078: **partial** (conf 0.72) — Same organism (human) and disease context (TNBC) as planned, but the actual contrast is a fibrillarin-depletion perturbation (knockdown vs control) profiled by ribosome profiling/translation efficiency — not tumor vs normal breast tissue, and there is no normal-tissue arm. The disease keyword matches, but the contrast axis and likely sample type (cell-based Ribo-seq, not tumor/normal tissue RNA-se
- GSE337023: **irrelevant** (conf 0.85) — Both TNBC and RNA-seq match on keywords, but the actual system is TNBC cell lines (HCC1806, MDA-MB-231), not breast tissue, and there is no normal breast tissue at all. The contrast axis is whole-genome-doubling mechanism vs parental aneuploid controls, not tumor vs normal. No tumor-vs-normal comparison can be derived, so it cannot serve as prior evidence for the planned contrast.
- GSE333496: **irrelevant** (conf 0.95) — Dataset profiles human antigen-induced Tr1 (CD4+ regulatory T) cells; TNBC is only mentioned in passing as one context where a Tr1 signature can be detected in tumor-resident CD4+ T cells. It has no breast tumor-vs-normal-tissue contrast, no breast tissue samples, and studies an immune-cell system rather than the planned TNBC tumor vs normal breast epithelium comparison. Only keyword overlap, not 

**Mismatch caveats:**
- requested tissue 'breast' not evident in GSE334884 metadata
- contrast terms 'triple-negative breast cancer tumor normal breast tissue' not evident in GSE334884 metadata

## Per-dataset reanalysis
- **GSE273078** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case set is empty (0 samples) and control set is empty (0 samples) — no differential-expression contrast is possible.; The data is a ribosome-profiling/SeRP knockdown experiment (condSC scrambled control, condSH2/condSH3 shRNA knockdowns, RPF/IN pairs) — none of these represent TNBC tumor or normal breast tissue, so t; Correctly, the exclusions themselves are justified: knockdown vs scrambled-control samples are off-contrast, and the *_DESeq2_norm / lfcMLE / p / padj entries are value/statistic columns, not samples.; The right conclusion is that this dataset cannot answer the TNBC-vs-normal question at all; approving a selection with two empty arms would misrepresent that no valid samples exist.
  - QC: 49 samples, median lib size 7,984,315, median zero-fraction 0.134.
  - Sample selection (dropped): case=0, control=0, excluded=49 — judge ⚠️ disagreed (Case set is empty (0 samples) and control set is empty (0 samples) — no differential-expression contrast is possible.; The data is a ribosome-profiling/SeRP knockdown experiment (condSC scrambled control, condSH2/condSH3 shRNA knockdowns, RPF/IN pairs) — none of these represent TNBC tumor or normal breast tissue, so t)
- **GSE309543** _[primary: GEO]_ — SKIPPED: single-cell RNA-seq study — set aside from bulk DE (single-cell analysis is planned)
- **GSE303509** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case set is empty (0 samples): no confirmed TNBC tumor was assigned to the case arm.; Control set is empty (0 samples): the dataset contains no normal breast tissue, so the control arm of the requested contrast cannot exist.; All 14 samples are tumor specimens; even if some are genuinely TNBC, excluding every one leaves no case arm, and there is nothing to contrast them against.; The exclusion reasons themselves confirm the contrast is infeasible here (no normal-tissue arm, TNBC subtype unconfirmed) — the correct conclusion is that this dataset does not support the requested c; If TNBC status could be established from associated metadata, the tumor samples might form a case arm, but a normal-tissue control would still have to come from another source; as proposed, the contra
  - QC: 14 samples, median lib size 748,738, median zero-fraction 0.411. Flags: low median library size (<1M counts) — typical of low-input / EV / exosome preps
  - Sample selection (dropped): case=0, control=0, excluded=14 — judge ⚠️ disagreed (Case set is empty (0 samples): no confirmed TNBC tumor was assigned to the case arm.; Control set is empty (0 samples): the dataset contains no normal breast tissue, so the control arm of the requested contrast cannot exist.)
- **GSE302168** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case set is empty (0 samples): no 'triple-negative breast cancer tumor' samples exist here.; Control set is empty (0 samples): no 'normal breast tissue' samples exist here.; All 12 samples belong to an in vitro GLI2-knockout/WT cell-line study with DMSO vs paclitaxel treatments (project 13051-JM) — a drug-response/knockout design, not a tumor-vs-normal-tissue design.; The exclusions are individually correct (none of these cell-line samples is a tumor biopsy or normal tissue), but the consequence is that the contrast simply is not represented in these samples; DE an; This dataset would instead support contrasts like WT_PTX vs WT_DMSO or GLI2_KO vs WT — not the requested clinical tumor/normal comparison.
  - QC: 12 samples, median lib size 89,941,280, median zero-fraction 0.021.
  - Sample selection (dropped): case=0, control=0, excluded=12 — judge ⚠️ disagreed (Case set is empty (0 samples): no 'triple-negative breast cancer tumor' samples exist here.; Control set is empty (0 samples): no 'normal breast tissue' samples exist here.)
- **GSE312563** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case set is empty (0 samples) and control set is empty (0 samples) — a DE contrast requires both arms to be non-empty, so this selection produces no runnable comparison.; The 12 samples are all in-vitro TNBC cell lines (MDA-MB-231, MDA-MB-436, Hs578T) under 24h drug treatment vs vehicle — none are primary tumor tissue and none are normal breast tissue, so the requested; The individual exclusions are correctly reasoned (cell line ≠ primary tumor; vehicle arm ≠ normal tissue), but the honest conclusion is that the contrast is unanswerable with these samples, not that a; If any contrast were to be salvaged from this dataset it would be a different question entirely (drug-treated vs vehicle within each TNBC cell line), which does not match the case=TNBC tumor / control
  - QC: 12 samples, median lib size 1,000,000, median zero-fraction 0.376. Flags: low median library size (<1M counts) — typical of low-input / EV / exosome preps
  - Sample selection (dropped): case=0, control=0, excluded=12 — judge ⚠️ disagreed (Case set is empty (0 samples) and control set is empty (0 samples) — a DE contrast requires both arms to be non-empty, so this selection produces no runnable comparison.; The 12 samples are all in-vitro TNBC cell lines (MDA-MB-231, MDA-MB-436, Hs578T) under 24h drug treatment vs vehicle — none are primary tumor tissue and none are normal breast tissue, so the requested)
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
  - Observed in GSE303509: low median library size (<1M counts) — typical of low-input / EV / exosome preps
  - Observed in GSE312563: low median library size (<1M counts) — typical of low-input / EV / exosome preps
  - Observed in GSE334884: low median library size (<1M counts) — typical of low-input / EV / exosome preps

<sub>schema 0.1.0 · real data only · projections are cited to the datasets above and are not a substitute for formal power analysis.</sub>