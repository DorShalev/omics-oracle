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
python -m omics_oracle "<prompt>" --h5 /path/to/human_matrix.h5   # enables counts + DE + concordance
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

## Configuration (all optional)
| Env var | Purpose |
|---|---|
| `OPENROUTER_API_KEY` | enable the LLM spec parser + narrative (else heuristic) |
| `OMICS_ORACLE_MODEL` | LiteLLM model id (default `openrouter/anthropic/claude-3.5-sonnet`) |
| `ARCHS4_H5` | path to a local ARCHS4 H5 for reprocessed counts (>30 GB; else metadata-only) |
| `NCBI_API_KEY`, `NCBI_EMAIL` | raise NCBI E-utilities rate limits |

## Data sources
- **NCBI GEO** (E-utilities) — the live, key-free dataset finder (RNA-seq & more).
- **ARCHS4** — uniformly-reprocessed RNA-seq counts via a local H5 (for DE + concordance).
- Expression Atlas, PRIDE (proteomics), ENCODE, pysradb — interface-ready, on the roadmap.

## Tests
```bash
pytest                      # offline: parser, ranking, DE (pydeseq2 on a real fixture), concordance, report
OMICS_ORACLE_LIVE=1 pytest  # also runs the live NCBI GEO integration test
```

## What it is not
Not a replacement for a formal power analysis or a statistician. It is a fast, cited
gut-check that surfaces prior evidence and likely pitfalls before you spend the money.
