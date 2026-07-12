"""Core data types for Omics Oracle.

Every type here is a plain dataclass with a ``to_dict()`` that yields JSON-serializable
output, so the whole engine can cross a process/HTTP boundary into another app without
callers touching internals. Keep this module dependency-light (stdlib only).
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Optional

# Bump when the serialized shape of OracleReport changes in a breaking way.
SCHEMA_VERSION = "0.1.0"

# Omics modalities we can reason about. Only "rna-seq" is wired end-to-end in v0.1.
OMICS_TYPES = ("rna-seq", "proteomics", "atac-seq", "scrna-seq", "metabolomics", "unknown")

Severity = str  # "blocking" | "major" | "minor"


def _clean(obj: Any) -> Any:
    """Recursively convert dataclasses/containers into JSON-safe primitives."""
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    return obj


@dataclass
class ExperimentSpec:
    """Structured description of a *planned* experiment, extracted from the prompt."""
    raw_prompt: str
    organism: Optional[str] = None          # canonical, e.g. "Homo sapiens"
    tissue: Optional[str] = None            # e.g. "colorectal tumor", "plasma EV"
    contrast_case: Optional[str] = None     # e.g. "tumor"
    contrast_control: Optional[str] = None  # e.g. "adjacent normal"
    omics_type: str = "unknown"
    platform: Optional[str] = None          # e.g. "Illumina"
    n_per_group: Optional[int] = None
    parser: str = "heuristic"               # "llm" | "heuristic"
    keywords: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class Mismatch:
    field: str
    expected: Any
    found: Any
    severity: Severity
    note: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class DatasetMatch:
    accession: str
    source: str                              # e.g. "GEO"
    title: str = ""
    organism: Optional[str] = None
    omics_type: str = "unknown"
    platform: Optional[str] = None
    n_samples: Optional[int] = None
    summary: str = ""
    url: str = ""
    similarity: float = 0.0                  # 0..1, explainable (see finder.py)
    mismatches: list[Mismatch] = field(default_factory=list)
    judge: Optional[dict] = None             # LLM relevance verdict (see judge.py), or None
    # match tier (how this study relates to the plan): "direct" | "proxy" | "pathway" | "none"
    match_tier: str = ""
    match_tier_reason: str = ""

    @property
    def blocking(self) -> bool:
        return any(m.severity == "blocking" for m in self.mismatches)

    @property
    def judged_irrelevant(self) -> bool:
        return bool(self.judge) and self.judge.get("verdict") == "irrelevant"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["blocking"] = self.blocking
        return d


@dataclass
class CountMatrix:
    """Genes x samples expression matrix plus per-sample metadata.

    ``data_type`` records what the values actually are, because it decides which DE method is
    valid: only ``raw`` counts may go to pydeseq2. Per-study supplementary tables are frequently
    normalized (TPM/FPKM/log), so this must be detected and honored, never assumed.
    """
    accession: str
    genes: list[str]
    samples: list[str]
    counts: Any                              # pandas.DataFrame (genes x samples)
    sample_meta: Any = None                  # pandas.DataFrame (samples x fields) or None
    provenance: str = ""
    data_type: str = "raw"                   # "raw" | "normalized" | "unknown"
    source_file: str = ""                    # the supplementary filename it came from
    source_db: str = "GEO"                   # which repository the counts came from
    tier: str = "primary"                    # "primary" (GEO per-study) | "clinical" (Maayan/ARCHS4)

    def to_dict(self) -> dict:
        return {
            "accession": self.accession,
            "n_genes": len(self.genes),
            "n_samples": len(self.samples),
            "data_type": self.data_type,
            "source_file": self.source_file,
            "source_db": self.source_db,
            "tier": self.tier,
            "provenance": self.provenance,
        }


@dataclass
class QCSummary:
    n_samples: int
    n_genes_detected: int
    median_library_size: float
    median_zero_fraction: float
    flags: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class DEResult:
    accession: str
    contrast: str                            # human-readable description of the contrast used
    n_case: int
    n_control: int
    table: Any = None                        # pandas.DataFrame: gene, log2fc, pvalue, padj
    qc: Optional[QCSummary] = None
    status: str = "ok"                       # "ok" | "skipped" | "error"
    note: str = ""
    source_db: str = "GEO"                   # repository the counts came from
    tier: str = "primary"                    # "primary" | "clinical"

    def n_de(self, padj: float = 0.05, lfc: float = 1.0) -> int:
        if self.table is None or len(self.table) == 0:
            return 0
        t = self.table
        return int(((t["padj"] < padj) & (t["log2fc"].abs() >= lfc)).sum())

    def to_dict(self) -> dict:
        return {
            "accession": self.accession,
            "contrast": self.contrast,
            "n_case": self.n_case,
            "n_control": self.n_control,
            "status": self.status,
            "note": self.note,
            "source_db": self.source_db,
            "tier": self.tier,
            "n_de_padj05_lfc1": self.n_de(),
            "qc": self.qc.to_dict() if self.qc else None,
        }


@dataclass
class GeneConcordance:
    gene: str
    n_datasets: int          # datasets where the gene was testable
    n_up: int
    n_down: int
    direction_agreement: float   # max(n_up, n_down) / n_datasets, 0..1
    mean_log2fc: float
    per_dataset: dict = field(default_factory=dict)   # accession -> log2fc

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class Concordance:
    n_datasets: int
    genes: list[GeneConcordance] = field(default_factory=list)
    note: str = ""

    def to_dict(self) -> dict:
        return {
            "n_datasets": self.n_datasets,
            "note": self.note,
            "top_concordant": [g.to_dict() for g in self.genes[:50]],
        }


@dataclass
class Projection:
    expected_n_de: Optional[tuple[int, int]] = None      # (low, high) range
    effect_size_prior: Optional[str] = None
    qc_pitfalls: list[str] = field(default_factory=list)
    power_read: str = ""
    feasibility: str = ""                                 # "encouraging" | "cautionary" | "thin-evidence"
    caveats: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        d = asdict(self)
        if self.expected_n_de is not None:
            d["expected_n_de"] = list(self.expected_n_de)
        return d


@dataclass
class OracleReport:
    spec: ExperimentSpec
    matches: list[DatasetMatch]
    de_results: list[DEResult]
    concordance: Optional[Concordance]
    projection: Projection
    narrative: str = ""
    comparison: Optional[dict] = None        # cross-dataset Venn/gene-search data (see compare.py)
    enrichment: Optional[dict] = None        # per-dataset + shared pathways (see enrich.py)
    heatmap: Optional[dict] = None           # gene x dataset log2FC matrix (see heatmap.py)
    schema_version: str = SCHEMA_VERSION
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "schema_version": self.schema_version,
            "spec": self.spec.to_dict(),
            "matches": [m.to_dict() for m in self.matches],
            "de_results": [d.to_dict() for d in self.de_results],
            "concordance": self.concordance.to_dict() if self.concordance else None,
            "comparison": self.comparison,
            "enrichment": self.enrichment,
            "heatmap": self.heatmap,
            "projection": self.projection.to_dict(),
            "narrative": self.narrative,
            "warnings": self.warnings,
        }

    def to_markdown(self) -> str:
        from .report import render_markdown
        return render_markdown(self)
