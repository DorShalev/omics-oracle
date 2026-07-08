"""Keep the test suite deterministic and offline w.r.t. the LLM.

Tests exercise the *logic* (parsing rules, ranking, DE, concordance, report), not the model.
Forcing the LLM backend off means project_experiment() uses the heuristic parser and the
templated narrative, so no test shells out to the Claude CLI or hits an API. The opt-in live
GEO test (OMICS_ORACLE_LIVE=1) still runs its network round-trip.
"""
import os
import pytest


@pytest.fixture(autouse=True)
def _no_llm(monkeypatch):
    monkeypatch.setenv("OMICS_ORACLE_LLM", "none")
    yield
