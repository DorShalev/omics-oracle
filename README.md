# Omics Oracle

**Sanity-check a *planned* omics experiment against public data — before you run it.**

Give Omics Oracle a natural-language description of an experiment you're about to run. It
parses the spec, mines public databases for prior similar experiments, ranks them (flagging
mismatches like *"GSE##### is HEK293, not colorectal"*), reanalyzes the closest studies, and
returns an **evidence-based, cited projection**: expected DE-gene counts, effect-size priors,
cross-study concordance, assay-specific QC pitfalls, and a rough power/feasibility read.

> **Honesty contract:** real data only. The engine never fabricates a dataset, a count, a
> p-value, or a projection. When evidence is thin or a step can't run, it says so and
> degrades gracefully. Every quantitative claim is tied to a source accession.

See [`DESIGN.md`](DESIGN.md) for the full architecture.

## Status — v0.1 working vertical slice
`prompt → GEO/ARCHS4 matched datasets → concordance report`, end-to-end, for RNA-seq.
Correctness and honest uncertainty are prioritized over breadth.

## Install
```bash
python -m venv .venv && . .venv/bin/activate
pip install -e .          # core engine
# optional extras:
pip install -e ".[llm]"      # LLM spec-parsing + narrative via OpenRouter/LiteLLM
pip install -e ".[archs4]"   # ARCHS4 reprocessed counts (needs a local H5)
```

## Use — CLI
```bash
python -m omics_oracle "RNA-seq of human colorectal tumor vs adjacent normal, Illumina, 20/arm"
python -m omics_oracle "RNA-seq of mouse liver, high-fat diet vs chow, n=8" --json
python -m omics_oracle "<prompt>" --no-judge      # skip the LLM relevance judge
```

### How counts work (DE + concordance)
By default, for each matched study the engine downloads **that study's own supplementary
expression table** from GEO (`*_counts.tsv.gz`, `*_TPM.txt.gz`, etc. — a few MB each), parses
it, and reanalyzes it. **No bulk database download.** It:
- prefers a **raw-counts** table over a normalized one, and **refuses to run DE on normalized
  data** (TPM/FPKM/log) — reporting it honestly instead of rounding-and-faking;
- groups samples into case/control from column names / clinical files, with an **LLM fallback**
  (via the configured backend) for messy naming;
- caps each download (`GEO_SUPPL_MAX_BYTES`, default 50 MB) and caches under `data/geo_suppl/`.

An ARCHS4 whole-database backend remains available as an explicit opt-in (uniformly reprocessed
counts, for heavy/offline use), but it is no longer the default:
```bash
python -m omics_oracle "<prompt>" --h5 data/archs4/human_gene_v2.latest.h5
python -m omics_oracle --download-h5 human --yes   # one 58 GB file, resumable
```

## Use — importable API (the stable boundary for embedding in another app)
```python
from omics_oracle import project_experiment

report = project_experiment(
    "RNA-seq of colorectal tumor vs adjacent normal, human, ~20 per arm",
    max_datasets=5,
    h5_path=None,      # ARCHS4 local H5; None => metadata-only reanalysis
    llm_model=None,    # None => heuristic parser + templated narrative
)
report.to_dict()       # JSON-serializable, schema-versioned (report.schema_version)
report.to_markdown()   # human-readable, cited
```

Lower-level stages are individually importable and swappable: `parse_spec`, `find_datasets`,
`reanalyze`, `aggregate`, `project`. New databases plug in by implementing the
`DatasetSource` protocol (`omics_oracle/sources/base.py`).

## LLM backend
The LLM is used only for spec parsing and narrative phrasing — **never to produce numbers**.
Backends are auto-selected (override with `OMICS_ORACLE_LLM` = `auto`|`claude-cli`|`litellm`|`none`):

1. **Claude CLI** — if the local `claude` binary is present (Claude Code subscription, **no API
   key, no per-token cost**). This is the default when no API key is set.
2. **LiteLLM / OpenRouter** — used when `OPENROUTER_API_KEY` (or `OPENAI_API_KEY`) is set.
3. **none** — no model available → deterministic heuristic parser + templated narrative. The
   engine still runs end-to-end; the numbers are identical because they come from data.

## Configuration (all optional)
| Env var | Purpose |
|---|---|
| `OMICS_ORACLE_LLM` | backend: `auto` (default), `claude-cli`, `litellm`, `none` |
| `OMICS_ORACLE_CLAUDE_MODEL` | Claude CLI model alias/id (e.g. `sonnet`, `opus`; else CLI default) |
| `OPENROUTER_API_KEY` / `OPENAI_API_KEY` | switch to the LiteLLM backend |
| `OMICS_ORACLE_MODEL` | LiteLLM model id (default `openrouter/anthropic/claude-3.5-sonnet`) |
| `ARCHS4_H5` | path to a local ARCHS4 H5 for reprocessed counts (>30 GB; else metadata-only) |
| `NCBI_API_KEY`, `NCBI_EMAIL` | raise NCBI E-utilities rate limits |

## Data sources
- **NCBI GEO** (E-utilities) — the live, key-free dataset finder (RNA-seq & more).
- **GEO supplementary tables** — each matched study's own expression table (the default counts
  path; small per-study downloads, raw-counts-aware, no bulk data).
- **ARCHS4** (opt-in) — uniformly-reprocessed RNA-seq counts from a local H5, for heavy/offline use.
- Expression Atlas, PRIDE (proteomics), ENCODE, pysradb — interface-ready, on the roadmap.

## Relevance judge (anti-fabrication)
On top of the deterministic mismatch flags, an **LLM judge** reads each top match's real
metadata and rules whether it *genuinely* fits the question — `relevant` / `partial` /
`irrelevant` / `uncertain`, with a reason. It catches semantic mismatches the keyword ranker
misses (disease mentioned only in passing, wrong contrast axis, cell line vs tissue), can only
*lower* trust in a match (never invent one), and answers `uncertain` when metadata is thin.
A dataset judged `irrelevant` is excluded from reanalysis. Disable with `--no-judge`.

## Does the sequencing platform matter?
For finding comparable RNA-seq, **no** — the specific Illumina instrument (NovaSeq / HiSeq /
NextSeq) does not gate or penalize matches. ARCHS4 reprocesses all Illumina RNA-seq uniformly,
so we gather every relevant RNA-seq regardless of model. Platform is recorded only as an
informational note, and flagged solely when a genuinely different *technology* (e.g. a
non-Illumina platform) is named. The real assay boundary — sequencing vs microarray vs
proteomics — is enforced by the omics-type check, not the platform.

## Download-size safeguards
- The only large asset is a **single** ARCHS4 H5 (~58 GB), and you usually skip it entirely
  via remote streaming (tens of MB per dataset). It is never "many files".
- `--download-h5` discloses the size, requires `--yes`, checks free disk (+5 GB headroom), and
  resumes partial downloads.
- The number of datasets reanalyzed is capped (`--reanalyze-top`, default 3), and each remote
  fetch is capped at `ARCHS4_MAX_SAMPLES` (default 200) samples.

## Tests
```bash
pytest                      # offline: parser, ranking, DE (pydeseq2 on a real fixture), concordance, report
OMICS_ORACLE_LIVE=1 pytest  # also runs the live NCBI GEO integration test
```

## What it is not
Not a replacement for a formal power analysis or a statistician. It is a fast, cited
gut-check that surfaces prior evidence and likely pitfalls before you spend the money.
