"""LLM relevance judge — a second, adversarial check that a matched dataset genuinely fits
the research question, rather than just sharing keywords.

The deterministic ranker (finder.py) catches *structural* mismatches (organism, assay, missing
tissue/contrast tokens). The judge catches *semantic* ones the tokens miss: a study that
mentions "colorectal" only in its background, a cell line masquerading as tissue, a different
contrast axis (drug response vs tumor-vs-normal), etc. It reads the dataset's REAL metadata
(title + summary from the source) and returns a grounded verdict with a reason.

Crucially, the judge never invents evidence: it is told to answer "uncertain" when the
metadata is insufficient, and it can only *lower* trust in a match, never fabricate one.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional

from .llm import LLM
from .types import DatasetMatch, ExperimentSpec

VERDICTS = ("relevant", "partial", "irrelevant", "uncertain")


@dataclass
class Verdict:
    accession: str
    verdict: str          # one of VERDICTS
    confidence: float     # 0..1, self-reported by the judge
    reason: str
    judged_by: str        # backend name, or "none"

    def to_dict(self) -> dict:
        return asdict(self)


_SYS = (
    "You are a skeptical bioinformatics reviewer verifying whether a PUBLIC dataset genuinely "
    "matches a PLANNED experiment. You are given the plan and the dataset's real metadata. "
    "Decide if the dataset could serve as prior evidence for the plan.\n"
    "Return ONLY a JSON object: {\"verdict\": one of [relevant, partial, irrelevant, uncertain], "
    "\"confidence\": 0..1, \"reason\": short string}.\n"
    "Rules: judge the ACTUAL biological system, contrast, and assay — not mere keyword overlap. "
    "A dataset that only mentions the disease in passing, uses a different contrast axis, or is a "
    "cell line when tissue was requested is 'partial' or 'irrelevant'. If the metadata is too "
    "thin to tell, answer 'uncertain'. Never assume facts not present in the metadata."
)


def judge_match(spec: ExperimentSpec, match: DatasetMatch, llm: Optional[LLM] = None) -> Verdict:
    llm = llm if llm is not None else LLM()
    if not llm.available:
        return Verdict(match.accession, "uncertain", 0.0,
                       "no LLM backend available to verify semantic relevance", "none")
    user = _prompt(spec, match)
    obj = llm.complete_json(_SYS, user)
    if not obj:
        return Verdict(match.accession, "uncertain", 0.0,
                       "judge did not return a parseable verdict", llm.backend_name)
    verdict = str(obj.get("verdict", "uncertain")).lower()
    if verdict not in VERDICTS:
        verdict = "uncertain"
    try:
        conf = float(obj.get("confidence", 0.0))
    except (TypeError, ValueError):
        conf = 0.0
    return Verdict(match.accession, verdict, max(0.0, min(1.0, conf)),
                   str(obj.get("reason", ""))[:400], llm.backend_name)


def judge_matches(spec: ExperimentSpec, matches: list[DatasetMatch],
                  llm: Optional[LLM] = None, limit: int = 5) -> dict[str, Verdict]:
    """Judge up to `limit` top matches; returns {accession: Verdict}."""
    llm = llm if llm is not None else LLM()
    out: dict[str, Verdict] = {}
    for m in matches[:limit]:
        out[m.accession] = judge_match(spec, m, llm)
    return out


def _prompt(spec: ExperimentSpec, match: DatasetMatch) -> str:
    return (
        "PLANNED EXPERIMENT:\n"
        f"  organism: {spec.organism}\n"
        f"  tissue/cell-type: {spec.tissue}\n"
        f"  contrast: {spec.contrast_case} vs {spec.contrast_control}\n"
        f"  omics: {spec.omics_type}\n"
        f"  free text: {spec.raw_prompt}\n\n"
        "CANDIDATE DATASET (real metadata):\n"
        f"  accession: {match.accession} ({match.source})\n"
        f"  organism: {match.organism}\n"
        f"  assay(from source): {match.omics_type}\n"
        f"  n_samples: {match.n_samples}\n"
        f"  title: {match.title}\n"
        f"  summary: {match.summary[:1500]}\n"
    )
