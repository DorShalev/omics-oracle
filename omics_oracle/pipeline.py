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
                       reanalyze_pool: int = 15,
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
    # Pull a larger candidate POOL than we display, so reanalysis can keep walking down the
    # ranking until it finds studies that actually deposit raw counts (many deposit only TPM).
    pool_limit = max(max_datasets, reanalyze_pool)
    try:
        pool = find_datasets(spec, sources=finders, limit=pool_limit)
    except Exception as e:
        pool = []
        warnings.append(f"Dataset finder failed ({type(e).__name__}: {e}); no matches.")
    matches = pool[:max_datasets]   # what the report displays + what the judge scores

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
    # Auto-expand: walk DOWN the ranked pool attempting reanalysis until we collect
    # `reanalyze_top` studies that actually yield DE (many deposit only normalized tables),
    # capped so we never hammer the network. Blocking/judged-irrelevant matches are skipped.
    ok = 0
    attempts = 0
    max_attempts = min(len(pool), max(reanalyze_top * 4, 8))
    for m in pool:
        if ok >= reanalyze_top or attempts >= max_attempts:
            break
        if m.blocking or m.judged_irrelevant:
            continue
        attempts += 1
        cm = fetch_counts(m, counts_srcs)
        for s in counts_srcs:                       # surface honest refusal reasons
            for w in getattr(s, "last_warnings", None) or []:
                if w not in warnings:
                    warnings.append(w)
        if cm is None:
            de_results.append(DEResult(accession=m.accession, contrast="(none)", n_case=0,
                                       n_control=0, status="skipped",
                                       note="no usable per-study counts table for this dataset"))
            continue
        res = reanalyze(cm, spec, llm=llm)
        de_results.append(res)
        if res.status == "ok":
            ok += 1
    # make sure every reanalyzed study is shown in the report, even if below the display cut
    shown = {m.accession for m in matches}
    matches = matches + [m for m in pool if m.accession in {d.accession for d in de_results}
                         and m.accession not in shown]

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
