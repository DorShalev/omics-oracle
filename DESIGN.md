# Omics Oracle — Design

> Given a natural-language description of a *planned* omics experiment, mine public
> databases for prior similar experiments and return an **evidence-based projection**
> that sanity-checks the experiment **before it is run**.

Status: **v0.1 — working vertical slice** (RNA-seq prompt → GEO/ARCHS4 matched datasets →
concordance report). Correctness and honest uncertainty are prioritized over breadth.

---

## 1. Problem & scope

A researcher is about to spend money and time on an experiment (e.g. *"RNA-seq of plasma
EVs from 20 colorectal-cancer patients vs 20 healthy controls, Illumina, human"*). Before
they run it, we want to answer:

- **Has something like this been done?** Which public datasets are genuinely comparable,
  and where do they *differ* (wrong tissue, wrong organism, wrong contrast)?
- **What should we expect?** Rough number of DE genes, effect-size priors, which signals
  are *reproducible* across prior studies, and known QC pitfalls for this assay.
- **Is it powered / feasible?** A rough read on N vs. expected effect sizes.

Non-goals (v0.1): we do not run the experiment, do not design primers/probes, and do not
claim to replace a statistician. Every projection is **cited** and **caveated**.

### Honesty contract (hard rule)
Real data only. The engine never fabricates a dataset, a count, a p-value, or a
projection. When evidence is thin or a step cannot run (e.g. reprocessed counts not
locally available), it says so explicitly and degrades to metadata-only rather than
inventing numbers. Every quantitative claim in a report carries its source accession(s).

---

## 2. Pipeline

```
prompt ─▶ (1) spec parser ─▶ (2) dataset finder ─▶ (3) fetcher ─▶ (4) reanalysis ─▶ (5) aggregator ─▶ (6) projector ─▶ (7) report
             ExperimentSpec     ranked+flagged        counts        per-DS DE          concordance       priors/power     cited md/json
```

1. **Spec parser** (`spec.py`) — LLM extracts a structured `ExperimentSpec`
   (organism, tissue/cell-type, perturbation/contrast, omics type, platform, expected N).
   Falls back to a deterministic keyword heuristic when no LLM key is configured; the
   result is labelled with `parser="llm"|"heuristic"` so downstream honesty is preserved.

2. **Dataset finder** (`finder.py` + `sources/`) — query public DBs, then **rank by
   similarity** and **explicitly flag mismatches** (organism, omics type, platform,
   tissue). Ranking is transparent: a weighted, additive score over comparable fields,
   with every deduction attributed to a named mismatch. No black-box embedding is required
   for v0.1 (optional later).

3. **Fetcher** (`fetch.py`) — pull *reprocessed* counts/intensities for the top 3–5
   matches. RNA-seq path uses **ARCHS4** uniformly-reprocessed counts.

4. **Per-dataset reanalysis** (`reanalyze.py`) — QC summary (library size, detected genes,
   zero-fraction) + **DE via pydeseq2** for the closest available contrast.

5. **Cross-experiment aggregator** (`aggregate.py`) — effect-size **meta-analysis** and
   per-gene **concordance** (how consistently a gene moves in the same direction across
   datasets), each finding annotated with its per-dataset context.

6. **Projector** (`project.py`) — assay-specific expected QC pitfalls, expected #DE and
   effect-size priors derived from the matched studies, and a rough power/feasibility read.

7. **Report** (`report.py`) — honest, cited Markdown/JSON with mismatch caveats.

---

## 3. Data sources & the reprocessed-counts reality

| Source | Role | Access in this build |
|---|---|---|
| **NCBI GEO** (E-utilities `esearch`/`esummary`, db=`gds`) | dataset finder (RNA-seq & others) | **remote, live** — the default finder. Rich metadata: accession, taxon, dataset type, platform (GPL), sample count, summary. |
| **ARCHS4** (uniformly reprocessed RNA-seq counts) | fetcher / reanalysis | via `archs4py` against a **local H5** (`ARCHS4_H5` / `--h5`). The full human H5 is >30 GB, so it is *not* bundled; when absent, reanalysis degrades to metadata-only and says so. |
| pysradb / SRA | finder enrichment (later) | optional |
| Expression Atlas, PRIDE (proteomics), ENCODE | later omics/DBs | interface-ready, not wired in v0.1 |

**Why GEO for finding, ARCHS4 for counts:** GEO's E-utilities give a reliable, key-free,
remotely queryable metadata index for *finding and ranking* studies. ARCHS4 gives
*uniformly reprocessed* counts (so cross-dataset DE is comparable) but is distributed as a
bulk H5. The `DatasetSource` interface (below) keeps these decoupled: a source can supply
metadata, counts, or both.

---

## 4. The stable function boundary (this must plug into another app)

The entire engine is reachable through **one importable entry point** returning a
JSON-serializable dataclass. Callers never touch internals.

```python
from omics_oracle import project_experiment, OracleReport

report: OracleReport = project_experiment(
    "RNA-seq of colorectal tumor vs adjacent normal, human, ~20 per arm",
    max_datasets=5,
    h5_path=None,          # ARCHS4 local H5; None => metadata-only reanalysis
    llm_model=None,        # None => heuristic parser + templated narrative
)
report.to_dict()           # fully serializable; stable schema (report.SCHEMA_VERSION)
report.to_markdown()
```

Internal seams, each independently testable and swappable:

```python
ExperimentSpec  = parse_spec(prompt, llm=...)                 # spec.py
matches         = find_datasets(spec, sources=[...])          # finder.py -> list[DatasetMatch]
counts          = fetch_counts(match, backend=...)            # fetch.py
de              = reanalyze(counts, spec)                      # reanalyze.py -> DEResult
agg             = aggregate([de, ...])                         # aggregate.py -> Concordance
proj            = project(spec, matches, agg)                  # project.py -> Projection
report          = build_report(spec, matches, [de...], agg, proj)   # report.py
```

`DatasetSource` protocol (`sources/base.py`):

```python
class DatasetSource(Protocol):
    name: str
    def search(self, spec: ExperimentSpec, limit: int) -> list[DatasetMatch]: ...
    def supports_counts(self) -> bool: ...
    def fetch_counts(self, match: DatasetMatch) -> CountMatrix | None: ...
```

### Core data types (dataclasses, all `to_dict()`-able)
- `ExperimentSpec` — organism, tissue, contrast(case/control terms), omics_type, platform,
  n_per_group, raw_prompt, parser.
- `DatasetMatch` — accession, source, title, organism, omics_type, platform, n_samples,
  summary, **similarity score (0–1)**, **mismatches: list[Mismatch]**.
- `Mismatch` — field, expected, found, severity (`blocking|major|minor`), note.
- `CountMatrix` — genes × samples (pandas), sample metadata, provenance.
- `DEResult` — per-dataset table (gene, log2fc, pvalue, padj), contrast used, QC summary.
- `Concordance` — per-gene direction agreement, mean/median effect, n_supporting, context.
- `Projection` — expected_n_de (range), effect_size_prior, qc_pitfalls, power_read.
- `OracleReport` — everything above + provenance + `SCHEMA_VERSION`.

---

## 5. Ranking & mismatch model

Similarity starts at 1.0 and is reduced by attributed penalties so the score is always
explainable:

| Field | Match test | Penalty if mismatch | Severity |
|---|---|---|---|
| organism | taxon equality | −0.5 | blocking |
| omics type | dataset-type / assay equality | −0.4 | blocking |
| tissue / cell-type | token overlap of tissue terms in title+summary | −0.25 | major |
| perturbation/contrast | token overlap of contrast terms | −0.2 | major |
| platform | GPL / instrument family equality | −0.05 | minor |
| N adequacy | n_samples ≥ 2·n_per_group | −0.05 | minor |

A **blocking** mismatch (wrong organism/omics) is surfaced first and prevents a dataset
from being used for the counts/DE step even if its text overlaps. This is the mechanism
behind caveats like *"GSE##### is HEK293, not the requested colorectal cancer."*

---

## 6. LLM usage & fallback

- **Provider:** OpenRouter via **LiteLLM** (`OPENROUTER_API_KEY`, `OMICS_ORACLE_MODEL`).
- **Where used:** (a) spec parsing (prompt→structured JSON, validated against the
  dataclass; malformed → retry once → heuristic), (b) report narrative phrasing only —
  **never** to produce numbers.
- **No key present:** the engine still runs end-to-end using the deterministic heuristic
  parser and a templated (non-LLM) narrative. Numbers are identical either way, because
  numbers come from data, not the model.

---

## 7. Testing strategy

- **Unit, offline:** heuristic parser, ranking/mismatch logic, concordance math, projector
  thresholds, report serialization — all run without network (fixtures).
- **Integration, live (opt-in via `OMICS_ORACLE_LIVE=1`):** GEO finder round-trip.
- **DE correctness:** pydeseq2 path validated on a small real RNA-seq count matrix fixture.
- The vertical slice is considered "working" when
  `project_experiment(<rna-seq prompt>)` returns a populated, cited `OracleReport` with
  real GEO matches and correct mismatch flags, with reanalysis either computed (H5 present)
  or honestly skipped.

---

## 8. Roadmap after the slice
1. ARCHS4 counts + pydeseq2 DE wired through when an H5 is available (production counts).
2. Proteomics via PRIDE API; Expression Atlas differential API as a second effect-size source.
3. Optional embedding-based similarity to complement the transparent additive score.
4. Pathway-level concordance (feeds naturally into multi-omics enrichment work).
5. Minimal web UI over the same `project_experiment` boundary.
