"""Robust parsing + typing of heterogeneous per-study expression tables.

GEO supplementary tables vary wildly: tab/comma/whitespace delimiters, quoted fields, a gene-id
column that may lack a header, Ensembl-with-version or symbol IDs, .gz or .xlsx, and values that
may be raw counts (possibly fractional, e.g. RSEM) or normalized (TPM/FPKM/log). This module turns
raw bytes into a clean genes x samples numeric DataFrame and classifies what the values are — the
latter decides whether pydeseq2 is valid downstream.
"""
from __future__ import annotations

import gzip
import io
import re
from typing import Optional

import numpy as np
import pandas as pd

# filename tokens (checked in priority order; "normalized" beats "counts")
_NORM_TOKENS = ("tpm", "fpkm", "rpkm", "cpm", "normali", "vst", "rlog", "zscore", "z-score",
                "log2", "logcpm", "scaled", "rpm")
_RAW_TOKENS = ("raw_count", "rawcount", "raw.count", "readcount", "read_count", "featurecount",
               "feature_count", "expected_count", "est_count", "genecount", "gene_count", "counts",
               "count")
_SKIP_TOKENS = ("differential", "_deg", "deg_", "de_results", "diffexp", "results", "clinical",
                "metadata", "samplesheet", "sample_sheet", "annotation", "readme", "filelist",
                "fpkm_tracking", "gtf", "peaks", "mirna", "isoform")

# Delimited bulk expression tables we can parse. `.tab` is common for FPKM/count matrices and
# was previously ignored — real bulk studies (e.g. per-sample FPKM matrices) were being missed.
_TABLE_EXT = (".csv", ".txt", ".tsv", ".tab", ".csv.gz", ".txt.gz", ".tsv.gz", ".tab.gz", ".xlsx")

# Single-cell deposition signatures — NOT bulk expression tables. A 10x study deposits the
# barcodes/features/matrix.mtx triplet, an .h5 feature-barcode matrix, or a Seurat .rds/.h5ad/.loom;
# these are cell x gene and must not be fed to bulk DE. We detect and set them aside honestly.
_SINGLE_CELL_TOKENS = ("barcodes.tsv", "features.tsv", "genes.tsv", "matrix.mtx",
                       "_bc_matrix", "filtered_feature", "raw_feature", "molecule_info")
_SINGLE_CELL_EXT = (".mtx", ".mtx.gz", ".h5", ".h5ad", ".loom", ".rds")


def is_single_cell_file(filename: str) -> bool:
    """True if the filename is a single-cell deposition (10x matrix, .h5/.rds/.h5ad/.loom)."""
    low = filename.lower()
    return low.endswith(_SINGLE_CELL_EXT) or any(t in low for t in _SINGLE_CELL_TOKENS)


def is_probably_expression_file(filename: str) -> bool:
    low = filename.lower()
    if is_single_cell_file(low):
        return False
    if not low.endswith(_TABLE_EXT):
        return False
    return not any(t in low for t in _SKIP_TOKENS)


def filename_score(filename: str) -> int:
    """Higher = more likely the count/expression matrix we want (raw preferred over normalized)."""
    low = filename.lower()
    if any(t in low for t in _SKIP_TOKENS):
        return -100
    score = 0
    if any(t in low for t in ("raw_count", "rawcount", "featurecount", "expected_count",
                              "readcount", "read_count")):
        score += 10
    elif "count" in low and not any(t in low for t in _NORM_TOKENS):
        score += 7
    if any(t in low for t in _NORM_TOKENS):
        score += 3            # usable, but normalized — lower priority than raw
    if any(t in low for t in ("matrix", "expression", "genes", "gene_")):
        score += 1
    return score


def decompress(raw: bytes, filename: str) -> bytes:
    if filename.lower().endswith(".gz") or raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw)
    return raw


def read_expression_table(raw: bytes, filename: str, min_cols: int = 2) -> pd.DataFrame:
    """Return a genes x samples numeric DataFrame, or raise ValueError if it can't be parsed.

    ``min_cols`` is the minimum number of sample columns required; use 1 for a per-sample file
    (one gene->count column) and the default 2 for a study-level matrix.
    """
    low = filename.lower()
    if low.endswith(".xlsx"):
        df = pd.read_excel(io.BytesIO(raw), index_col=0)
        return _finalize(df)

    data = decompress(raw, filename)
    text = data.decode("utf-8", errors="replace")
    # skip obvious comment/preamble lines beginning with '#'
    lines = [ln for ln in text.splitlines() if ln.strip()]
    lines = [ln for ln in lines if not ln.startswith("#")]
    text = "\n".join(lines)

    best = None
    for kw in ({"sep": "\t"}, {"sep": ","}, {"sep": r"\s+", "engine": "python"}):
        try:
            df = pd.read_csv(io.StringIO(text), index_col=0, quotechar='"', **kw)
        except Exception:
            continue
        try:
            fin = _finalize(df)
        except ValueError:
            continue
        numeric_frac = np.isfinite(fin.to_numpy(dtype=float)).mean() if fin.size else 0.0
        cand = (fin.shape[1], round(float(numeric_frac), 3), fin)
        if best is None or cand[:2] > best[:2]:
            best = cand
    if best is None or best[0] < min_cols or best[1] < 0.8:
        raise ValueError(f"could not parse a numeric genes x samples table from {filename}")
    return best[2]


def _finalize(df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.shape[1] < 1:
        raise ValueError("no columns")
    # strip quotes/whitespace from index + column labels
    df = df.copy()
    df.index = [_clean_label(i) for i in df.index]
    df.columns = [_clean_label(c) for c in df.columns]
    # keep only numeric sample columns (drop annotation columns like gene_name/length/biotype)
    num = df.apply(pd.to_numeric, errors="coerce")
    keep = [c for c in num.columns if num[c].notna().mean() >= 0.8]
    if len(keep) < 1:
        raise ValueError("no numeric sample columns")
    num = num[keep]
    # index must look like gene ids (mostly non-numeric labels)
    idx = pd.Index([str(i) for i in num.index])
    nonnumeric = np.mean([not _looks_numeric(i) for i in idx]) if len(idx) else 0.0
    if nonnumeric < 0.5:
        raise ValueError("index does not look like gene identifiers")
    num = num[~idx.duplicated(keep="first")]
    num.index = idx[~idx.duplicated(keep="first")]
    return num.dropna(how="all")


def classify_data_type(df: pd.DataFrame, filename: str) -> str:
    """Return 'raw' | 'normalized' | 'unknown'.

    Only 'raw' (integer OR fractional count-scale, e.g. RSEM/salmon expected_counts) is valid for
    pydeseq2 downstream. We combine the filename signal with a value signal and, when they
    disagree, return 'unknown' (refuse) rather than guess — honesty over convenience.
    """
    fname = _filename_signal(filename)
    value = _value_signal(df)
    if fname == "normalized":
        return "normalized"          # depositors label normalized tables deliberately; trust it
    if fname == "raw":
        return "raw" if value != "normalized" else "unknown"   # conflict => refuse
    return value                     # no filename signal => go by the values


def _filename_signal(filename: str) -> Optional[str]:
    low = filename.lower()
    if any(t in low for t in _NORM_TOKENS):
        return "normalized"
    if any(t in low for t in ("raw_count", "rawcount", "raw.count", "featurecount",
                              "expected_count", "est_count", "readcount", "read_count")):
        return "raw"
    return None


def _value_signal(df: pd.DataFrame) -> str:
    vals = df.to_numpy(dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return "unknown"
    if (vals < 0).any():
        return "normalized"                              # negatives => log/z-scored
    colsums = df.sum(axis=0).to_numpy(dtype=float)
    colsums = colsums[np.isfinite(colsums) & (colsums > 0)]
    cv = float(colsums.std() / colsums.mean()) if colsums.size >= 2 and colsums.mean() else 1.0
    vmax = float(vals.max())
    near_int = float(np.mean(np.abs(vals - np.round(vals)) < 1e-3))
    # TPM/CPM: every library normalized to ~the same total => column sums cluster tightly
    if cv < 0.02 and vmax > 100:
        return "normalized"
    # raw counts: integer-valued, big dynamic range, library sizes vary
    if near_int > 0.9 and vmax > 100:
        return "raw"
    # expected_counts (RSEM/salmon): fractional but count-scale (large max, varying library sizes)
    if vmax > 1000 and cv > 0.05:
        return "raw"
    # small fractional non-negative values => almost certainly normalized (TPM/log)
    if vmax <= 100:
        return "normalized"
    return "unknown"


# ---- helpers ---------------------------------------------------------------------------

def _clean_label(x) -> str:
    return str(x).strip().strip('"').strip("'").strip()


def _looks_numeric(s: str) -> bool:
    return bool(re.fullmatch(r"[-+]?\d*\.?\d+([eE][-+]?\d+)?", str(s).strip()))
