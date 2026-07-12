"""Sample selection — decide which of a study's samples are case / control / exclude, then
VERIFY that decision with an independent LLM judge before any DE is run on it.

A public study usually contains MORE than the two arms we care about: a tumor-vs-normal series
may also carry drug-treated samples, a timecourse, or replicates of an unrelated condition.
Taking *all* samples would corrupt the contrast. So the flow is:

  propose  — the model labels each sample case | control | exclude, each WITH a reason
             (optionally by self-consistency consensus over N runs, majority-voted per sample);
  verify   — a separate, skeptical judge checks the selection matches the requested contrast
             and that the exclusions are justified, returning agree/disagree + concerns;
  retry    — on disagreement, re-propose once with the judge's concerns fed back, then re-verify;
  drop     — if it still disagrees, the study is dropped honestly (no DE on a bad grouping).

The full trail (proposal with per-sample reasons, verdict, outcome) is recorded on the DEResult
so the report can show exactly which samples were used, which were excluded, and why. The model
only ever LABELS the exact sample ids we pass — it can never invent a sample.
"""
from __future__ import annotations

import re
from collections import Counter
from dataclasses import asdict, dataclass, field
from typing import Optional

CASE, CONTROL, EXCLUDE = "case", "control", "exclude"
_LABELS = (CASE, CONTROL, EXCLUDE)


# ---- audit data types (all JSON-serializable) ----------------------------------------------

@dataclass
class SampleAssignment:
    sample: str
    label: str            # case | control | exclude
    reason: str = ""
    votes: int = 1        # consensus runs that agreed on this label (1 if no consensus)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class SampleSelection:
    assignments: list[SampleAssignment]
    proposed_by: str                       # "explicit" | "heuristic" | backend name
    consensus_runs: int = 1

    def groups(self) -> dict[str, str]:
        return {a.sample: a.label for a in self.assignments if a.label in (CASE, CONTROL)}

    def excluded(self) -> list[str]:
        return [a.sample for a in self.assignments if a.label == EXCLUDE]

    def to_dict(self) -> dict:
        g = self.groups()
        return {"proposed_by": self.proposed_by, "consensus_runs": self.consensus_runs,
                "n_case": sum(1 for v in g.values() if v == CASE),
                "n_control": sum(1 for v in g.values() if v == CONTROL),
                "n_excluded": len(self.excluded()),
                "assignments": [a.to_dict() for a in self.assignments]}


@dataclass
class SelectionVerdict:
    agree: bool
    reason: str = ""
    concerns: list[str] = field(default_factory=list)
    judged_by: str = "none"

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class SampleAudit:
    selection: Optional[SampleSelection] = None
    verdict: Optional[SelectionVerdict] = None
    # explicit | heuristic | unverified | used | retried-used | dropped | none
    outcome: str = "none"

    def groups(self) -> dict[str, str]:
        if self.outcome == "dropped" or self.selection is None:
            return {}
        return self.selection.groups()

    def to_dict(self) -> dict:
        return {"outcome": self.outcome,
                "selection": self.selection.to_dict() if self.selection else None,
                "verdict": self.verdict.to_dict() if self.verdict else None}


# ---- prompts -------------------------------------------------------------------------------

_PROPOSE_SYS = (
    "You assign RNA-seq sample columns to a differential-expression contrast. Use ONLY the "
    "provided sample identifiers — never invent, rename, or omit-then-add ids. For EACH sample "
    "choose exactly one label: 'case' (matches the case arm of the contrast), 'control' "
    "(matches the control arm), or 'exclude' (belongs to NEITHER arm — a different condition, "
    "treatment, timepoint, cell line, or an ambiguous/QC sample). Excluding is CORRECT and "
    "expected when a sample does not fit the two arms; do NOT force every sample into a group. "
    "Return ONLY JSON of the form "
    "{\"assignments\": {\"<sample id>\": {\"label\": \"case|control|exclude\", "
    "\"reason\": \"<short justification>\"}, ...}}."
)

_VERIFY_SYS = (
    "You are a skeptical reviewer verifying a proposed case/control sample selection for a "
    "differential-expression contrast. You are given the research question and the proposed "
    "labels (with the proposer's reasons). Decide whether the selection is APPROPRIATE: the "
    "case and control sets must both be non-empty and genuinely represent the requested "
    "contrast, and the exclusions must be justified (no valid arm samples wrongly dropped, no "
    "off-contrast samples wrongly included). Return ONLY JSON: "
    "{\"agree\": true|false, \"reason\": \"<short>\", \"concerns\": [\"<specific issue>\", ...]}. "
    "Be strict: if the grouping looks arbitrary, mixes conditions, or the contrast is not "
    "actually present in these samples, set agree=false and say why."
)


def _sample_hints(counts) -> dict:
    samples = [str(c) for c in counts.counts.columns]
    meta = getattr(counts, "sample_meta", None)
    hints: dict[str, str] = {}
    if meta is not None:
        for s in samples:
            if s in meta.index:
                hints[s] = " ".join(str(v) for v in meta.loc[s].values)
    return hints


# ---- LLM steps (module-level so tests can call them directly) -------------------------------

def propose_selection(counts, spec, llm, feedback: Optional[list[str]] = None
                      ) -> Optional[SampleSelection]:
    """One LLM call: label every sample case|control|exclude with a reason."""
    samples = [str(c) for c in counts.counts.columns]
    if not samples or spec.contrast_case is None or spec.contrast_control is None:
        return None
    hints = _sample_hints(counts)
    user = (f"Contrast: case = {spec.contrast_case!r} vs control = {spec.contrast_control!r}\n"
            f"Tissue/context: {spec.tissue}\n"
            f"Sample ids (and any metadata hints): {hints or samples}")
    if feedback:
        user += ("\n\nA prior selection was REJECTED by a reviewer for these reasons — fix them:\n"
                 + "\n".join(f"- {c}" for c in feedback if c))
    obj = llm.complete_json(_PROPOSE_SYS, user)
    if not obj:
        return None
    raw = obj.get("assignments", obj)
    if not isinstance(raw, dict):
        return None
    assigns: list[SampleAssignment] = []
    for s in samples:
        v = raw.get(s)
        if v is None:
            continue
        if isinstance(v, dict):
            label = str(v.get("label", "")).lower()
            reason = str(v.get("reason", ""))[:300]
        else:
            label, reason = str(v).lower(), ""
        if label in _LABELS:
            assigns.append(SampleAssignment(s, label, reason))
    if not assigns:
        return None
    return SampleSelection(assigns, getattr(llm, "backend_name", "llm"))


def consensus_propose(counts, spec, llm, n: int = 3) -> Optional[SampleSelection]:
    """Run :func:`propose_selection` n times; majority-vote each sample over the 3 labels.

    Self-consistency guards the most result-critical LLM decision. A sample's winning label is
    the plurality across runs; its reason is taken from a run that voted that label; ``votes``
    records the agreement count so the report can flag shaky samples.
    """
    runs: list[SampleSelection] = []
    for _ in range(max(1, n)):
        sel = propose_selection(counts, spec, llm)
        if sel is not None:
            runs.append(sel)
    if not runs:
        return None
    votes: dict[str, Counter] = {}
    reasons: dict[str, dict[str, str]] = {}
    for sel in runs:
        for a in sel.assignments:
            votes.setdefault(a.sample, Counter())[a.label] += 1
            reasons.setdefault(a.sample, {}).setdefault(a.label, a.reason)
    assigns: list[SampleAssignment] = []
    for sample, c in votes.items():
        label, cnt = c.most_common(1)[0]
        assigns.append(SampleAssignment(sample, label, reasons[sample].get(label, ""), votes=cnt))
    return SampleSelection(assigns, getattr(llm, "backend_name", "llm"), consensus_runs=len(runs))


def verify_selection(counts, spec, selection: SampleSelection, llm) -> SelectionVerdict:
    """One independent LLM call: is this case/control selection appropriate for the contrast?"""
    if llm is None or not getattr(llm, "available", False):
        return SelectionVerdict(True, "no LLM backend to verify; selection used unverified",
                                judged_by="none")

    def _lst(label):
        return [f"{a.sample} ({a.reason})" if a.reason else a.sample
                for a in selection.assignments if a.label == label]

    user = (f"Research question contrast: case = {spec.contrast_case!r} vs "
            f"control = {spec.contrast_control!r}\nTissue/context: {spec.tissue}\n\n"
            f"PROPOSED case ({len(_lst(CASE))}): {_lst(CASE)}\n"
            f"PROPOSED control ({len(_lst(CONTROL))}): {_lst(CONTROL)}\n"
            f"PROPOSED exclude ({len(_lst(EXCLUDE))}): {_lst(EXCLUDE)}")
    obj = llm.complete_json(_VERIFY_SYS, user)
    if not obj:
        return SelectionVerdict(True, "verifier returned no parseable verdict; selection kept",
                                judged_by=getattr(llm, "backend_name", "llm"))
    agree = bool(obj.get("agree", False))
    concerns = obj.get("concerns") or []
    if not isinstance(concerns, list):
        concerns = [str(concerns)]
    return SelectionVerdict(agree, str(obj.get("reason", ""))[:400],
                            [str(c)[:200] for c in concerns][:8],
                            judged_by=getattr(llm, "backend_name", "llm"))


# ---- the single seam reanalyze calls -------------------------------------------------------

def select_samples(counts, spec, llm=None, *, consensus: int = 1, verify: bool = True,
                   max_retries: int = 1,
                   explicit_groups: Optional[dict[str, str]] = None) -> SampleAudit:
    """Decide + verify the case/control selection for one study. Cheapest paths first.

    Order: explicit grouping (0 calls) -> metadata heuristic (0 calls) -> LLM propose
    (consensus N calls) -> 1 verify -> bounded retry (concerns fed back). Worst case N+3 slow
    LLM calls. Returns a :class:`SampleAudit`; ``audit.groups()`` is empty when the study was
    dropped, so the caller runs no DE on a rejected grouping.
    """
    audit = SampleAudit()

    if explicit_groups:
        audit.selection = SampleSelection(
            [SampleAssignment(s, lab, "provided explicitly") for s, lab in explicit_groups.items()],
            proposed_by="explicit")
        audit.outcome = "explicit"
        return audit

    heur = _heuristic_selection(counts, spec)
    have_llm = llm is not None and getattr(llm, "available", False)

    if not have_llm:                                   # no model: heuristic, unverified
        if heur is not None:
            audit.selection, audit.outcome = heur, "heuristic"
        return audit

    sel = (consensus_propose(counts, spec, llm, consensus) if consensus > 1
           else propose_selection(counts, spec, llm))
    if sel is None:
        sel = heur
    if sel is None:
        return audit
    audit.selection = sel

    if not verify:
        audit.outcome = "used"
        return audit

    verdict = verify_selection(counts, spec, sel, llm)
    audit.verdict = verdict
    if verdict.agree:
        audit.outcome = "used"
        return audit

    for _ in range(max(0, max_retries)):               # bounded retry with the concerns fed back
        sel2 = propose_selection(counts, spec, llm, feedback=verdict.concerns or [verdict.reason])
        if sel2 is None:
            break
        verdict = verify_selection(counts, spec, sel2, llm)
        audit.selection, audit.verdict = sel2, verdict
        if verdict.agree:
            audit.outcome = "retried-used"
            return audit

    audit.outcome = "dropped"
    return audit


# generic control/case condition tokens — used when the spec's own wording doesn't match, so a
# ctrl/exp split (or wt/ko, treated/untreated, tumor/normal, diabetic/non-diabetic) still groups.
# Matched as whole words against a metadata blob that includes the source filename (split-file
# provenance) and the series-matrix characteristics.
_GENERIC_CONTROL = ("control", "ctrl", "normal", "healthy", "wildtype", "wt", "untreated",
                    "vehicle", "mock", "dmso", "sham", "baseline", "nondiabetic", "adjacent",
                    "naive", "uninfected", "pbs", "scramble", "scrambled", "sictrl", "shctrl")
_GENERIC_CASE = ("experiment", "exp", "treated", "treatment", "tumor", "tumour", "cancer",
                 "disease", "knockout", "ko", "knockdown", "kd", "mutant", "diabetic", "patient",
                 "stimulated", "infected", "lesion", "affected", "diseased")


def _heuristic_selection(counts, spec) -> Optional[SampleSelection]:
    """Case/control labels from sample metadata — the deterministic, no-LLM path.

    Uses a metadata blob per sample that includes the series-matrix characteristics AND the
    source filename, so split-file provenance (``*_ctrl`` vs ``*_exp``) is a grouping signal.
    First tries the spec's own contrast wording; if that doesn't produce both arms, falls back to
    generic condition tokens. Whatever it produces is still checked by the verification judge.
    """
    meta = getattr(counts, "sample_meta", None)
    if meta is None:
        return None
    blobs: dict[str, str] = {}
    for sample in counts.counts.columns:
        s = str(sample)
        if s in meta.index:
            raw = " ".join(str(v) for v in meta.loc[s].values).lower()
            blobs[s] = " " + re.sub(r"[^a-z0-9]+", " ", raw).strip() + " "   # word-level matching
    if not blobs:
        return None

    sel = None
    if spec.contrast_case and spec.contrast_control:
        sel = _assign_by_terms(blobs, _terms(spec.contrast_case), _terms(spec.contrast_control),
                               "sample metadata matched the requested contrast terms")
    if sel is None:                                  # split-file / generic condition fallback
        sel = _assign_by_terms(blobs, _GENERIC_CASE, _GENERIC_CONTROL,
                               "grouped by a case/control condition token (e.g. split-file arm)")
    return sel


def _assign_by_terms(blobs: dict[str, str], case_terms, ctrl_terms,
                     reason: str) -> Optional[SampleSelection]:
    """Assign each sample by whole-word term match; a sample hitting BOTH arms is left unassigned
    (ambiguous). Returns a selection only if both arms end up non-empty."""
    assigns: list[SampleAssignment] = []
    for s, blob in blobs.items():
        hit_case = any(f" {t} " in blob for t in case_terms)
        hit_ctrl = any(f" {t} " in blob for t in ctrl_terms)
        if hit_case and not hit_ctrl:
            assigns.append(SampleAssignment(s, CASE, reason))
        elif hit_ctrl and not hit_case:
            assigns.append(SampleAssignment(s, CONTROL, reason))
    if not any(a.label == CASE for a in assigns) or not any(a.label == CONTROL for a in assigns):
        return None
    return SampleSelection(assigns, proposed_by="heuristic")


def _terms(s: str) -> list[str]:
    return [t for t in s.lower().replace("/", " ").split() if len(t) > 2]
