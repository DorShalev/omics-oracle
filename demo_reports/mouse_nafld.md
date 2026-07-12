# Omics Oracle — experiment sanity check

> _Mouse liver RNA-seq, high-fat diet vs normal chow, to study NAFLD/steatosis_

**Parsed spec** (`llm` parser): organism=**Mus musculus**, tissue=**liver**, omics=**rna-seq**, contrast=**high-fat diet vs normal chow**, platform=**—**, N/group=**—**.

**Warnings:** LLM backend: claude-cli (spec parsing + narrative phrasing only; numbers come from data, not the model).; GSE289033: no supplementary expression table found (files: GSE289033_RAW.tar, filelist.txt); GSE313258: supplementary table is normalized (GSE313258_AllSamples_Countsmatrix_TPM.csv.gz); DE requires raw counts and will be skipped (not faked); GSE302447: could not list supplementary files (HTTPError); GSE302445: no supplementary expression table found (files: GSE302445_MUP_uPA_HFD_TOFA_liver_FPKM.tab.gz); GSE302444: no supplementary expression table found (files: GSE302444_WT_CDAA-HFD_TOFA_liver_FPKM.tab.gz); GSE302442: no supplementary expression table found (files: GSE302442_PPARa-KO_HFD_TOFA_liver_FPKM.tab.gz); GSE302440: no supplementary expression table found (files: GSE302440_WT_HFD_TOFA_liver_FPKM.tab.gz); GSE334406: no supplementary expression table found (files: GSE334406_RAW.tar, filelist.txt)

## Narrative
Prior public data only weakly supports this plan as specified. The closest matches by title/organism/omics were mouse RNA-seq datasets GSE306578, GSE289033, and GSE313258 (all top-similarity), with GSE302447 and GSE302445 ranking lower due to a major tissue mismatch. However, none of these could be turned into a high-fat-diet vs normal-chow differential-expression check: seven of the reanalyzed studies had no usable per-study counts table, and the one liver-lipid-metabolism dataset that did (GSE313258) turned out to be a dose-response design (vehicle plus dose arms) with no chow-vs-HFD contrast present, so its samples were dropped rather than force-fit to your comparison. As a result, no differential-expression count or effect-size prior could be estimated from prior data. One concrete QC caveat carried over from GSE313258 is a low median library size (<1M counts), and general pitfalls to watch include normalization choice and batch/condition confounding. Note also that no per-group sample size was specified, so no power read is possible beyond the rough guideline of ~6–12 samples per group for moderate effects. Overall feasibility is cautionary: the biological topic is well represented in public mouse liver RNA-seq, but no retrieved dataset provided a runnable matching contrast to anchor expectations, so treat any projections as unconfirmed until a dataset with raw counts and a true HFD-vs-chow grouping is located.

## Matched datasets (ranked)
| Rank | Accession | Sim | Tier | Organism | Omics | N | Blocking? | Judge | Mismatches |
|---|---|---|---|---|---|---|---|---|---|
| 1 | [GSE306578](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306578) | 1.0 | none | Mus musculus | rna-seq | 9 | — | irrelevant | — |
| 2 | [GSE289033](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289033) | 1.0 | proxy | Mus musculus | rna-seq | 9 | — | partial | — |
| 3 | [GSE313258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313258) | 1.0 | proxy | Mus musculus | rna-seq | 54 | — | partial | — |
| 4 | [GSE302447](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302447) | 0.75 | proxy | Mus musculus | rna-seq | 67 | — | partial | tissue(major) |
| 5 | [GSE302445](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302445) | 0.75 | proxy | Mus musculus | rna-seq | 16 | — | — | tissue(major) |
| 6 | [GSE302444](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302444) | 0.75 | proxy | Mus musculus | rna-seq | 9 | — | — | tissue(major) |
| 7 | [GSE302442](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302442) | 0.75 | proxy | Mus musculus | rna-seq | 18 | — | — | tissue(major) |
| 8 | [GSE302440](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302440) | 0.75 | proxy | Mus musculus | rna-seq | 8 | — | — | tissue(major) |
| 9 | [GSE334406](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334406) | 0.55 | proxy | Mus musculus | rna-seq | 6 | — | — | tissue(major); contrast(major) |

**Relevance verification (LLM judge):**
- GSE306578: **irrelevant** (conf 0.9) — Same organism (mouse) and tissue (liver) and broad assay (RNA-seq, though single-cell), but the biological system and contrast are entirely different: oncogenic β-catenin in hepatocytes / normal vs chronically-damaged liver, a hepatic oncogenesis model. There is no dietary intervention (high-fat vs chow) and no NAFLD/steatosis axis, so it cannot serve as prior evidence for the planned diet-induced
- GSE289033: **partial** (conf 0.6) — Matches organism (mouse), assay (RNA-seq), tissue (liver/hepatocytes), and disease (MAFLD/NAFLD, HFD-induced). However, the study's actual contrast axis is a bone-liver intervention (osteocyte-derived EV ablation / miR-483-5p antagonism), not high-fat diet vs normal chow. HFD is only the background disease model, and the 9-sample design likely compares EV/miRNA perturbation rather than diet, so it
- GSE313258: **partial** (conf 0.72) — Matches organism (mouse), tissue (liver), and assay (RNA-seq) with a metabolic/hepatic-lipid focus relevant to NAFLD/steatosis. However, the actual contrast axis is DINP phthalate dose (0/1.5/15/150 mg/kg), not high-fat diet vs normal chow — all mice are on a 45% HFD background, so the planned HFD-vs-chow contrast is not directly present (the 0 mg/kg control is still HFD). Useful as related prior 
- GSE302447: **partial** (conf 0.55) — Right organism (mouse) and assay (RNA-seq), and the diet-induced obesity / lipid-metabolism theme is biologically adjacent to NAFLD/steatosis. But the primary contrast axis is an oral drug treatment vs control in DIO mice, not the requested high-fat-diet vs normal-chow diet contrast, and the tissue is never stated (liver unconfirmed — could be adipose, whole-body, or multiple). As a thin SuperSeri

**Mismatch caveats:**
- requested tissue 'liver' not evident in GSE302447 metadata
- requested tissue 'liver' not evident in GSE302445 metadata
- requested tissue 'liver' not evident in GSE302444 metadata
- requested tissue 'liver' not evident in GSE302442 metadata
- requested tissue 'liver' not evident in GSE302440 metadata
- requested tissue 'liver' not evident in GSE334406 metadata
- contrast terms 'high-fat diet normal chow' not evident in GSE334406 metadata

## Per-dataset reanalysis
- **GSE289033** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE313258** _[primary: GEO]_ — SKIPPED: sample selection was rejected by the verification judge and dropped rather than run on a questionable grouping — Case set is empty (0 samples) and control set is empty (0 samples) — a DE contrast requires both arms to be non-empty; nothing is runnable.; The proposer's per-sample exclusion reasons are individually correct: 'CONTROL_REP*' are vehicle/untreated controls of a dose-response study (not normal-chow dietary controls), and '1.5/15/150_REP*' a; However, the correct conclusion is that the requested contrast does not exist in these samples — there is neither a high-fat-diet arm nor a normal-chow arm. This should be flagged as 'contrast not pre; No liver tissue is confirmed for any retained sample either, but this is moot since no samples are retained.; Do not attempt to salvage this by mapping vehicle control to 'normal chow' or a dose arm to 'high-fat diet' — that would fabricate a contrast that mixes an unrelated dose-response design.
  - QC: 54 samples, median lib size 974,444, median zero-fraction 0.369. Flags: low median library size (<1M counts) — typical of low-input / EV / exosome preps
  - Sample selection (dropped): case=0, control=0, excluded=54 — judge ⚠️ disagreed (Case set is empty (0 samples) and control set is empty (0 samples) — a DE contrast requires both arms to be non-empty; nothing is runnable.; The proposer's per-sample exclusion reasons are individually correct: 'CONTROL_REP*' are vehicle/untreated controls of a dose-response study (not normal-chow dietary controls), and '1.5/15/150_REP*' a)
- **GSE302447** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE302445** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE302444** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE302442** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE302440** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset
- **GSE334406** _[primary: GEO]_ — SKIPPED: no usable per-study counts table for this dataset

## Projection
- **Power / feasibility:** Planned N per group not specified; specify it for a power read. As a rule of thumb, bulk RNA-seq needs ≈6–12 per group to detect moderate effects reliably.
- **Overall feasibility:** **cautionary**.
- **Expected QC pitfalls:**
  - Check for rRNA / globin carryover, especially in blood/plasma inputs.
  - Verify library size and 3' bias; low-input EV/exosome RNA often yields shallow, high-duplication libraries.
  - Normalization choice (median-of-ratios vs TMM vs CPM) changes DE calls — fix one and report it.
  - Batch/lane confounding with condition is the most common false-signal source; block or randomize.
  - Observed in GSE313258: low median library size (<1M counts) — typical of low-input / EV / exosome preps
- **Caveats:**
  - No prior study yielded raw counts for DE (7 had no usable per-study counts table; 1 had no assignable case/control grouping); DE-count and effect-size priors are not estimated. Try raising --reanalyze-top to reach studies deeper in the ranking that may deposit raw counts.

<sub>schema 0.1.0 · real data only · projections are cited to the datasets above and are not a substitute for formal power analysis.</sub>