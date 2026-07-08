"""Omics Oracle — evidence-based projections that sanity-check a planned omics experiment.

Stable public API (this is the boundary other apps import against):

    from omics_oracle import project_experiment, OracleReport
    report = project_experiment("RNA-seq of colorectal tumor vs adjacent normal, human, ~20/arm")
    report.to_dict()        # JSON-serializable, schema-versioned
    report.to_markdown()

Lower-level seams (spec/finder/reanalyze/aggregate/project) are also importable for callers
that want to run stages individually.
"""
from .pipeline import project_experiment
from .types import (
    OracleReport, ExperimentSpec, DatasetMatch, Mismatch, CountMatrix, DEResult,
    Concordance, GeneConcordance, Projection, QCSummary, SCHEMA_VERSION,
)
from .spec import parse_spec
from .finder import find_datasets, score_match
from .reanalyze import reanalyze, compute_qc
from .aggregate import aggregate
from .project import project

__version__ = "0.1.0"

__all__ = [
    "project_experiment", "OracleReport", "ExperimentSpec", "DatasetMatch", "Mismatch",
    "CountMatrix", "DEResult", "Concordance", "GeneConcordance", "Projection", "QCSummary",
    "SCHEMA_VERSION", "parse_spec", "find_datasets", "score_match", "reanalyze",
    "compute_qc", "aggregate", "project", "__version__",
]
