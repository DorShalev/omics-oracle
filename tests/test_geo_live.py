"""Opt-in live integration test against NCBI GEO. Run with OMICS_ORACLE_LIVE=1."""
import os
import pytest

from omics_oracle.spec import parse_spec
from omics_oracle.llm import LLM
from omics_oracle.finder import find_datasets
from omics_oracle.sources.geo import GEOSource

LIVE = os.environ.get("OMICS_ORACLE_LIVE") == "1"


class _NoLLM(LLM):
    @property
    def available(self): return False


@pytest.mark.skipif(not LIVE, reason="set OMICS_ORACLE_LIVE=1 to run network tests")
def test_geo_finds_colorectal_rnaseq():
    spec = parse_spec("RNA-seq of human colorectal cancer tumor vs normal", llm=_NoLLM())
    matches = find_datasets(spec, sources=[GEOSource()], limit=5)
    assert matches, "expected at least one GEO match"
    top = matches[0]
    assert top.accession.startswith("GSE")
    assert top.organism  # metadata populated
    # ranking must expose similarity and mismatch structure
    assert 0.0 <= top.similarity <= 1.0
