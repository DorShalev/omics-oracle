"""Offline tests for ARCHS4 remote source index logic (series membership, SC filtering).

No network: we inject in-memory index arrays and bypass the cache-file check.
"""
import numpy as np

from omics_oracle.sources.archs4_remote import ARCHS4RemoteSource


def _src_with_index(series, geo, scprob=None):
    s = ARCHS4RemoteSource(url="http://unused", index_dir="/nonexistent")
    s._series = np.array(series, dtype=object)
    s._geo = np.array(geo, dtype=object)
    s._genes = np.array(["A", "B"], dtype=object)
    s._scprob = None if scprob is None else np.array(scprob)
    # make _load_index a no-op since we injected arrays
    s._load_index = lambda: None
    return s


def test_series_membership_simple():
    s = _src_with_index(["GSE1", "GSE2", "GSE1"], ["GSM1", "GSM2", "GSM3"])
    idx = s._sample_rows("GSE1")
    assert list(idx) == [0, 2]


def test_series_membership_multi_series_joined():
    # a sample can belong to several series, joined by comma/semicolon
    s = _src_with_index(["GSE1,GSE9", "GSE2", "GSE3;GSE1"], ["GSM1", "GSM2", "GSM3"])
    idx = s._sample_rows("GSE1")
    assert list(idx) == [0, 2]


def test_single_cell_samples_filtered():
    s = _src_with_index(["GSE1", "GSE1", "GSE1"], ["GSM1", "GSM2", "GSM3"],
                        scprob=[0.01, 0.99, 0.2])
    idx = s._sample_rows("GSE1")
    assert list(idx) == [0, 2]   # the p=0.99 single-cell sample is dropped


def test_no_match_returns_empty():
    s = _src_with_index(["GSE1", "GSE2"], ["GSM1", "GSM2"])
    assert s._sample_rows("GSE_absent").size == 0


def test_supports_counts_false_without_index():
    s = ARCHS4RemoteSource(url="http://unused", index_dir="/definitely/not/here")
    assert s.supports_counts() is False
