"""series_matrix parsing + column matching, and split-file / generic heuristic grouping."""
import numpy as np
import pandas as pd

from omics_oracle.series_matrix import _parse, build_sample_meta
from omics_oracle.selection import _heuristic_selection
from omics_oracle.types import CountMatrix, ExperimentSpec


def test_parse_series_matrix_header():
    lines = [
        '!Sample_title\t"1-002"\t"1-008"\t"n-1"',
        '!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"',
        '!Sample_source_name_ch1\t"Tumor"\t"Tumor"\t"Normal"',
        '!Sample_characteristics_ch1\t"tissue: Tumor"\t"tissue: Tumor"\t"tissue: Normal"',
        '!Sample_characteristics_ch1\t"treatment: ctrl"\t"treatment: exp"\t"treatment: none"',
    ]
    df = _parse(lines)
    assert list(df.index) == ["GSM1", "GSM2", "GSM3"]
    assert df.loc["GSM1", "title"] == "1-002"
    assert df.loc["GSM1", "characteristics"] == "tissue: Tumor; treatment: ctrl"
    assert df.loc["GSM3", "source"] == "Normal"


def test_match_columns_by_title_and_gsm():
    meta = _parse([
        '!Sample_title\t"1-002"\t"1-008"',
        '!Sample_geo_accession\t"GSM1"\t"GSM2"',
        '!Sample_source_name_ch1\t"Tumor"\t"Normal"',
        '!Sample_characteristics_ch1\t"tissue: Tumor"\t"tissue: Normal"',
    ])
    m = build_sample_meta(["1-002", "GSM2"], meta)          # one by title, one by GSM id
    assert m.loc["1-002", "gsm"] == "GSM1"
    assert m.loc["1-002", "characteristics"] == "tissue: Tumor"
    assert m.loc["GSM2", "gsm"] == "GSM2"
    assert m.loc["GSM2", "source_name"] == "Normal"


def _cm(cols, meta_df):
    df = pd.DataFrame(np.arange(len(cols) * 4).reshape(4, len(cols)) + 1.0,
                      index=[f"G{i}" for i in range(4)], columns=cols)
    return CountMatrix(accession="GSE1", genes=list(df.index), samples=cols, counts=df,
                       sample_meta=meta_df)


def test_heuristic_uses_series_matrix_characteristics():
    # real condition labels from the series matrix let the spec-term path group correctly
    meta = pd.DataFrame({"characteristics": ["tissue: tumor", "tissue: tumor",
                                             "tissue: normal", "tissue: normal"]},
                        index=["a", "b", "c", "d"])
    spec = ExperimentSpec(raw_prompt="x", contrast_case="tumor", contrast_control="normal",
                          omics_type="rna-seq")
    sel = _heuristic_selection(_cm(["a", "b", "c", "d"], meta), spec)
    assert sel.groups() == {"a": "case", "b": "case", "c": "control", "d": "control"}


def test_heuristic_split_file_provenance_generic_fallback():
    # spec says tumor/normal but the study split files are ctrl/exp -> generic fallback groups them
    meta = pd.DataFrame({"source_file": ["s_ctrl.txt", "s_ctrl.txt", "s_exp.txt", "s_exp.txt"]},
                        index=["s1", "s2", "s3", "s4"])
    spec = ExperimentSpec(raw_prompt="x", contrast_case="tumor", contrast_control="normal",
                          omics_type="rna-seq")
    sel = _heuristic_selection(_cm(["s1", "s2", "s3", "s4"], meta), spec)
    assert sel.groups() == {"s1": "control", "s2": "control", "s3": "case", "s4": "case"}


def test_heuristic_ambiguous_sample_left_unassigned():
    # a sample whose blob hits BOTH arms is not forced into either
    meta = pd.DataFrame({"characteristics": ["tumor", "normal", "tumor normal"]},
                        index=["a", "b", "c"])
    spec = ExperimentSpec(raw_prompt="x", contrast_case="tumor", contrast_control="normal",
                          omics_type="rna-seq")
    sel = _heuristic_selection(_cm(["a", "b", "c"], meta), spec)
    assert sel.groups() == {"a": "case", "b": "control"}     # c omitted (ambiguous)
