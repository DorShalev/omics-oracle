"""Judge logic with a stub LLM (no network)."""
from omics_oracle.types import ExperimentSpec, DatasetMatch
from omics_oracle.judge import judge_match, judge_matches
from omics_oracle.llm import LLM


class _StubLLM(LLM):
    def __init__(self, obj):
        self._obj = obj
    @property
    def available(self): return True
    @property
    def backend_name(self): return "stub"
    def complete_json(self, system, user): return self._obj


def _spec():
    return ExperimentSpec(raw_prompt="x", organism="Homo sapiens", tissue="colorectal",
                          contrast_case="tumor", contrast_control="normal", omics_type="rna-seq")


def test_judge_relevant():
    m = DatasetMatch(accession="GSE1", source="GEO", title="colorectal tumor vs normal")
    v = judge_match(_spec(), m, llm=_StubLLM({"verdict": "relevant", "confidence": 0.9, "reason": "matches"}))
    assert v.verdict == "relevant" and v.confidence == 0.9 and v.judged_by == "stub"


def test_judge_irrelevant_and_bounds():
    m = DatasetMatch(accession="GSE2", source="GEO", title="HEK293 drug screen")
    v = judge_match(_spec(), m, llm=_StubLLM({"verdict": "irrelevant", "confidence": 5, "reason": "cell line"}))
    assert v.verdict == "irrelevant"
    assert v.confidence == 1.0  # clamped to [0,1]


def test_judge_unknown_verdict_defaults_uncertain():
    m = DatasetMatch(accession="GSE3", source="GEO")
    v = judge_match(_spec(), m, llm=_StubLLM({"verdict": "banana", "confidence": 0.5, "reason": "?"}))
    assert v.verdict == "uncertain"


def test_judge_no_llm_is_uncertain_not_fabricated():
    class _Off(LLM):
        @property
        def available(self): return False
    m = DatasetMatch(accession="GSE4", source="GEO")
    v = judge_match(_spec(), m, llm=_Off())
    assert v.verdict == "uncertain" and v.judged_by == "none"


def test_judge_matches_limit():
    ms = [DatasetMatch(accession=f"GSE{i}", source="GEO") for i in range(10)]
    verdicts = judge_matches(_spec(), ms, llm=_StubLLM({"verdict": "relevant", "confidence": 1, "reason": "x"}), limit=3)
    assert len(verdicts) == 3
