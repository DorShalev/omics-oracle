"""Sample-selection propose/verify/retry seam — offline, LLM stubbed at complete_json."""
import numpy as np
import pandas as pd

from omics_oracle.types import CountMatrix, ExperimentSpec
from omics_oracle.selection import (select_samples, propose_selection, verify_selection,
                                     consensus_propose, SampleSelection, SampleAssignment)


class FakeLLM:
    """Dispatches on the system prompt: verify prompts say 'verifying', propose prompts don't."""
    available = True
    backend_name = "fake"

    def __init__(self, propose=None, verify=None):
        self._p, self._v = propose, verify

    def complete_json(self, system, user):
        return self._v if "verif" in system.lower() else self._p


def _counts(cols):
    df = pd.DataFrame(np.arange(len(cols) * 4).reshape(4, len(cols)) + 1.0,
                      index=[f"G{i}" for i in range(4)], columns=cols)
    return CountMatrix(accession="GSE1", genes=list(df.index), samples=list(df.columns), counts=df)


def _spec():
    return ExperimentSpec(raw_prompt="x", organism="Homo sapiens", tissue="colon",
                          contrast_case="tumor", contrast_control="normal", omics_type="rna-seq")


def test_explicit_groups_bypass_llm():
    audit = select_samples(_counts(["a", "b"]), _spec(), llm=FakeLLM(),
                           explicit_groups={"a": "case", "b": "control"})
    assert audit.outcome == "explicit"
    assert audit.groups() == {"a": "case", "b": "control"}


def test_propose_supports_exclude():
    prop = {"assignments": {
        "t1": {"label": "case", "reason": "tumor"},
        "n1": {"label": "control", "reason": "normal"},
        "drug1": {"label": "exclude", "reason": "drug-treated, off-contrast"}}}
    sel = propose_selection(_counts(["t1", "n1", "drug1"]), _spec(), FakeLLM(propose=prop))
    assert sel.groups() == {"t1": "case", "n1": "control"}
    assert sel.excluded() == ["drug1"]


def test_verify_agree_uses_selection():
    prop = {"assignments": {"t1": {"label": "case"}, "n1": {"label": "control"}}}
    llm = FakeLLM(propose=prop, verify={"agree": True, "reason": "clean tumor/normal split"})
    audit = select_samples(_counts(["t1", "n1"]), _spec(), llm=llm)
    assert audit.outcome == "used" and audit.verdict.agree
    assert audit.groups() == {"t1": "case", "n1": "control"}


def test_verify_disagree_then_dropped():
    prop = {"assignments": {"x1": {"label": "case"}, "x2": {"label": "control"}}}
    llm = FakeLLM(propose=prop,
                  verify={"agree": False, "reason": "arbitrary split",
                          "concerns": ["case/control not distinguishable from these ids"]})
    audit = select_samples(_counts(["x1", "x2"]), _spec(), llm=llm, max_retries=1)
    assert audit.outcome == "dropped"
    assert audit.groups() == {}                 # rejected grouping is NOT used for DE
    assert audit.verdict.agree is False


def test_no_llm_falls_back_to_heuristic_unverified():
    # metadata lets the heuristic split without any model
    df = pd.DataFrame(np.arange(8).reshape(4, 2) + 1.0,
                      index=[f"G{i}" for i in range(4)], columns=["s_tumor", "s_normal"])
    meta = pd.DataFrame({"desc": ["tumor sample", "normal sample"]},
                        index=["s_tumor", "s_normal"])
    cm = CountMatrix(accession="GSE1", genes=list(df.index), samples=list(df.columns),
                     counts=df, sample_meta=meta)

    class Dead:
        available = False
        backend_name = "none"

    audit = select_samples(cm, _spec(), llm=Dead())
    assert audit.outcome == "heuristic"
    assert audit.groups() == {"s_tumor": "case", "s_normal": "control"}


def test_consensus_majority_vote_over_three_labels():
    runs = [
        SampleSelection([SampleAssignment("a", "case"), SampleAssignment("b", "control"),
                         SampleAssignment("c", "exclude")], "fake"),
        SampleSelection([SampleAssignment("a", "case"), SampleAssignment("b", "control"),
                         SampleAssignment("c", "exclude")], "fake"),
        SampleSelection([SampleAssignment("a", "control"), SampleAssignment("b", "control"),
                         SampleAssignment("c", "case")], "fake"),
    ]
    calls = iter(runs)

    class SeqLLM:
        available = True
        backend_name = "fake"

        def complete_json(self, system, user):
            sel = next(calls)
            return {"assignments": {x.sample: {"label": x.label} for x in sel.assignments}}

    sel = consensus_propose(_counts(["a", "b", "c"]), _spec(), SeqLLM(), n=3)
    g = sel.groups()
    assert g["a"] == "case" and g["b"] == "control"   # a: 2/3 case, b: 3/3 control
    assert sel.excluded() == ["c"]                     # c: 2/3 exclude
