"""Match-tier classification — how each matched study relates to the user's plan.

  Tier 1 "direct"  : same organism + tissue + contrast (a real disease/design match).
  Tier 2 "proxy"   : related biology — partial metadata match (other cancer, a model system,
                     same tissue different contrast, etc.).
  Tier 3 "pathway" : the VIRTUAL LAYER — the study may not match the disease, but its DE genes
                     enrich the SAME pathways as the user's signal. Insight beyond the disease.
  "none"           : neither a metadata nor a pathway match.

Tier 1/2 come from the deterministic mismatch flags + the LLM relevance judge; Tier 3 comes
from the enrichment result (studies that contribute to shared/surfaced pathways). Pathway can
also annotate a proxy study to highlight its cross-disease value.
"""
from __future__ import annotations

from typing import Optional

from .types import DatasetMatch, ExperimentSpec

TIER_ORDER = {"direct": 0, "proxy": 1, "pathway": 2, "none": 3}


def _pathway_accessions(enrichment: Optional[dict]) -> set[str]:
    accs: set[str] = set()
    if not enrichment:
        return accs
    for p in enrichment.get("shared_pathways", []):
        accs |= set(p.get("sources_significant", []))
    for p in enrichment.get("surfaced", []):
        accs |= set(p.get("significant_in", []))
        accs |= {s.get("source") for s in p.get("subthreshold_in", [])}
    accs.discard("USER")
    return accs


def classify_tiers(spec: ExperimentSpec, matches: list[DatasetMatch],
                   enrichment: Optional[dict] = None) -> list[DatasetMatch]:
    """Assign match_tier / match_tier_reason to each match in place; return the list."""
    pathway_accs = _pathway_accessions(enrichment)

    for m in matches:
        verdict = (m.judge or {}).get("verdict")
        has_major = any(x.severity in ("blocking", "major") for x in m.mismatches)
        in_pathway = m.accession in pathway_accs

        if not m.blocking and verdict == "relevant" and not has_major:
            tier, reason = "direct", "same organism, tissue and contrast (judge: relevant)"
        elif not m.blocking and verdict != "irrelevant":
            reason = "related biology / partial metadata match"
            if verdict == "partial":
                reason = "partial metadata match (judge: partial)"
            tier = "proxy"
        elif in_pathway:
            tier, reason = "pathway", ("different disease/assay, but its DE genes enrich the same "
                                       "pathways as your signal (virtual layer)")
        else:
            tier = "none"
            reason = ("blocked on organism/assay" if m.blocking
                      else "not a metadata or pathway match")

        # a metadata proxy that is ALSO a pathway match: keep tier=proxy but flag the extra value
        if in_pathway and tier == "proxy":
            reason += "; also shares your enriched pathways"

        m.match_tier, m.match_tier_reason = tier, reason
    return matches
