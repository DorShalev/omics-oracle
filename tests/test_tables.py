"""Table parsing + data-type detection — the honesty-critical unit tests (no network)."""
import gzip
import io

import numpy as np
import pandas as pd

from omics_oracle.tables import (read_expression_table, classify_data_type,
                                 is_probably_expression_file, filename_score)


def _counts_frame(seed=0, integer=True, scale=1.0):
    rng = np.random.default_rng(seed)
    genes = [f"ENSG{i:08d}" for i in range(50)]
    cols = [f"S{i}" for i in range(6)]
    m = rng.integers(0, 5000, size=(50, 6)).astype(float) * scale
    if not integer:
        m = m + rng.random((50, 6))  # fractional (RSEM-like)
    return pd.DataFrame(m if not integer else m.astype(int), index=genes, columns=cols)


# ---- classifier ------------------------------------------------------------------------

def test_integer_counts_are_raw():
    assert classify_data_type(_counts_frame(integer=True), "GSE1_counts.tsv") == "raw"


def test_rsem_fractional_counts_are_raw():
    df = _counts_frame(integer=False)
    assert classify_data_type(df, "GSE1_raw_counts.txt") == "raw"


def test_tpm_colsums_constant_are_normalized():
    # every column normalized to the same total (~1e6) => TPM/CPM
    df = _counts_frame(integer=False)
    df = df / df.sum(axis=0) * 1e6
    assert classify_data_type(df, "GSE1_expr.txt") == "normalized"


def test_log_values_with_negatives_are_normalized():
    df = _counts_frame() .apply(lambda c: np.log2(c + 1)) - 5  # introduces negatives
    assert classify_data_type(df, "GSE1_expr.txt") == "normalized"


def test_small_float_values_are_normalized():
    df = _counts_frame(integer=False) / 5000.0 * 20  # small floats, max ~20
    assert classify_data_type(df, "GSE1_expr.txt") == "normalized"


def test_filename_tpm_overrides_integer_values():
    df = _counts_frame(integer=True)
    assert classify_data_type(df, "GSE1_TPM.txt") == "normalized"


def test_filename_raw_but_normalized_values_is_unknown():
    df = _counts_frame(integer=False) / 5000.0 * 20  # looks normalized
    # filename claims raw counts, values disagree => refuse (unknown), never fake
    assert classify_data_type(df, "GSE1_raw_counts.txt") == "unknown"


# ---- filename heuristics ----------------------------------------------------------------

def test_expression_file_recognition():
    assert is_probably_expression_file("GSE1_raw_counts.tsv.gz")
    assert not is_probably_expression_file("GSE1_differential_expression.txt.gz")
    assert not is_probably_expression_file("GSE1_clinical_details.txt.gz")
    assert not is_probably_expression_file("GSE1.bam")


def test_raw_scores_above_normalized():
    assert filename_score("GSE1_raw_counts.tsv.gz") > filename_score("GSE1_TPM.txt.gz")
    assert filename_score("GSE1_clinical.txt.gz") < 0


# ---- parser round-trips (each real-world format) ----------------------------------------

def _df():
    return pd.DataFrame({"S1": [10, 200, 3000], "S2": [12, 180, 2500]},
                        index=["ENSG1", "ENSG2", "ENSG3"])


def test_parse_tsv_gz():
    df = _df()
    raw = gzip.compress(df.to_csv(sep="\t").encode())
    out = read_expression_table(raw, "x_counts.tsv.gz")
    assert list(out.columns) == ["S1", "S2"] and out.shape == (3, 2)


def test_parse_csv():
    df = _df()
    out = read_expression_table(df.to_csv().encode(), "x_counts.csv")
    assert out.loc["ENSG2", "S1"] == 200


def test_parse_quoted_space_delimited():
    # mimic GSE313035: quoted labels, space separated, header lacks index-col name
    body = '"S1" "S2"\n"A1BG" 12.9 21.0\n"A2M" 373.9 175.5\n'
    out = read_expression_table(body.encode(), "x_TPM.txt")
    assert "A1BG" in out.index and out.shape[1] == 2


def test_parse_rejects_nontable():
    try:
        read_expression_table(b"not,a,table\njust text here\n", "x.txt")
    except ValueError:
        return
    # a garbage file should not masquerade as an expression matrix
    assert False, "expected ValueError on non-table input"
