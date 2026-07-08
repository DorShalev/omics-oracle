from omics_oracle.spec import parse_spec
from omics_oracle.llm import LLM


class _NoLLM(LLM):
    @property
    def available(self):  # force heuristic path deterministically
        return False


def test_heuristic_extracts_core_fields():
    spec = parse_spec(
        "RNA-seq of colorectal tumor vs adjacent normal in human patients, Illumina, 20 per group",
        llm=_NoLLM(),
    )
    assert spec.parser == "heuristic"
    assert spec.organism == "Homo sapiens"
    assert spec.omics_type == "rna-seq"
    assert spec.platform == "Illumina"
    assert spec.n_per_group == 20
    assert spec.contrast_case and spec.contrast_control
    assert "colorectal" in (spec.tissue or "")


def test_heuristic_proteomics_and_mouse():
    spec = parse_spec("LC-MS proteomics of mouse liver, treated vs control, n=6", llm=_NoLLM())
    assert spec.organism == "Mus musculus"
    assert spec.omics_type == "proteomics"
    assert spec.n_per_group == 6


def test_unknown_omics_is_labeled():
    spec = parse_spec("some vague plan about cells", llm=_NoLLM())
    assert spec.omics_type == "unknown"
    assert spec.parser == "heuristic"
