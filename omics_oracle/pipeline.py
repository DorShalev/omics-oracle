"""The orchestrator — the single stable entry point another app imports.

    from omics_oracle import project_experiment
    report = project_experiment("RNA-seq of colorectal tumor vs normal, human, ~20/arm")

Every stage is defensive: a failure in fetch/reanalysis degrades to metadata-only and is
recorded in ``report.warnings`` rather than aborting the run.
"""
from __future__ import annotations

from typing import Optional, Sequence

from .llm import LLM
from .spec import parse_spec
from .finder import find_datasets
from .fetch import fetch_counts
from .reanalyze import reanalyze
from .aggregate import aggregate
from .project import project
from .report import build_report
from .sources.base import DatasetSource
from .sources.geo import GEOSource
from .sources.archs4 import ARCHS4Source
from .types import OracleReport, DEResult


def project_experiment(prompt: str,
                       max_datasets: int = 5,
                       reanalyze_top: int = 3,
                       h5_path: Optional[str] = None,
                       llm_model: Optional[str] = None,
                       judge: bool = True,
                       judge_top: int = 5,
                       finder_sources: Optional[Sequence[DatasetSource]] = None,
                       counts_sources: Optional[Sequence[DatasetSource]] = None) -> OracleReport:
    """Run the full prompt -> projection pipeline and return a serializable OracleReport."""
    warnings: list[str] = []
    llm = LLM(model=llm_model) if llm_model else LLM()
    if llm.available:
        warnings.append(f"LLM backend: {llm.backend_name} (spec parsing + narrative phrasing only; "
                        "numbers come from data, not the model).")
    else:
        warnings.append("No LLM backend available (no Claude CLI, no API key) — using the "
                        "heuristic spec parser and a templated narrative (numbers are unaffected).")

    # 1. spec
    spec = parse_spec(prompt, llm=llm)

    # 2. find + rank
    finders = list(finder_sources) if finder_sources is not None else [GEOSource()]
    finder_names = {getattr(s, "name", "") for s in finders}
    if spec.omics_type in ("proteomics", "metabolomics") and finder_names <= {"GEO"}:
        warnings.append(f"GEO is not the authoritative index for {spec.omics_type}; matches may "
                        "be off-assay. PRIDE (proteomics) / MetaboLights are the right sources "
                        "(not wired in v0.1) — treat these results as weak.")
    try:
        matches = find_datasets(spec, sources=finders, limit=max_datasets)
    except Exception as e:
        matches = []
        warnings.append(f"Dataset finder failed ({type(e).__name__}: {e}); no matches.")

    # 2b. LLM relevance judge — semantic anti-fabrication check on the top matches.
    if judge and llm.available and matches:
        from .judge import judge_matches
        try:
            verdicts = judge_matches(spec, matches, llm=llm, limit=judge_top)
            for m in matches:
                v = verdicts.get(m.accession)
                if v is not None:
                    m.judge = v.to_dict()
        except Exception as e:
            warnings.append(f"Relevance judge failed ({type(e).__name__}); matches shown unjudged.")
    elif judge and not llm.available:
        warnings.append("Relevance judge skipped (no LLM backend) — matches are ranked by the "
                        "deterministic mismatch model only, not semantically verified.")

    # 3-4. fetch per-study counts + reanalyze the top non-blocking matches.
    # Default counts backend: each matched study's OWN supplementary expression table from GEO
    # (small per-study download). A local ARCHS4 H5 is an explicit opt-in override via --h5.
    if counts_sources is not None:
        counts_srcs = list(counts_sources)
    elif h5_path:
        counts_srcs = [ARCHS4Source(h5_path=h5_path)]
    else:
        from .sources.geo_suppl import GeoSupplSource
        counts_srcs = [GeoSupplSource()]

    de_results: list[DEResult] = []
    # reanalyze only matches that clear BOTH the structural filter and the judge
    usable = [m for m in matches if not m.blocking and not m.judged_irrelevant]
    for m in usable[:reanalyze_top]:
        cm = fetch_counts(m, counts_srcs)
        # surface each source's honest refusal reasons (no silent skips)
        for s in counts_srcs:
            for w in getattr(s, "last_warnings", None) or []:
                if w not in warnings:
                    warnings.append(w)
        if cm is None:
            de_results.append(DEResult(accession=m.accession, contrast="(none)", n_case=0,
                                       n_control=0, status="skipped",
                                       note="no usable per-study counts table for this dataset"))
            continue
        de_results.append(reanalyze(cm, spec, llm=llm))

    # 5. concordance
    concordance = aggregate(de_results) if any(d.status == "ok" for d in de_results) else None

    # 6. projection
    projection = project(spec, matches, de_results, concordance)

    # 7. report
    return build_report(spec, matches, de_results, concordance, projection,
                        llm=llm, warnings=warnings)


def _supports(src: DatasetSource) -> bool:
    try:
        return src.supports_counts()
    except Exception:
        return False
