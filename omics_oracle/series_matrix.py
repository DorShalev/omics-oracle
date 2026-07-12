"""Real per-sample metadata from a GEO ``series_matrix.txt.gz``.

Count-table columns are often cryptic (``1-002``); the series matrix carries the truth — each
sample's ``!Sample_title`` / ``!Sample_source_name`` / ``!Sample_characteristics`` (tissue,
disease, treatment, genotype…). We fetch and parse just the header (stopping before the giant
expression block), then match it back onto the count columns so the sample selector groups on
real condition labels instead of guessing from opaque ids.
"""
from __future__ import annotations

import gzip
import re
from typing import Optional

import pandas as pd

_MATRIX_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def _matrix_url(accession: str) -> str:
    num = accession[3:]
    prefix = num[:-3] if len(num) > 3 else ""
    return (f"{_MATRIX_BASE}/GSE{prefix}nnn/{accession}/matrix/"
            f"{accession}_series_matrix.txt.gz")


def fetch_sample_metadata(accession: str, session, timeout: float = 60.0,
                          max_lines: int = 4000) -> Optional[pd.DataFrame]:
    """Return a DataFrame indexed by GSM with columns [title, source, characteristics], or None.

    Streams the gzip and stops at ``!series_matrix_table_begin`` so only the (small) metadata
    header is read, never the full expression block.
    """
    url = _matrix_url(accession)
    try:
        with session.get(url, stream=True, timeout=timeout) as r:
            r.raise_for_status()
            r.raw.decode_content = False              # the file is gzip on disk; decode it ourselves
            sample_lines: list[str] = []
            with gzip.GzipFile(fileobj=r.raw) as gz:
                for i, raw_line in enumerate(gz):
                    if i > max_lines:
                        break
                    line = raw_line.decode("utf-8", "replace")
                    if line.startswith("!series_matrix_table_begin"):
                        break
                    if line.startswith("!Sample_"):
                        sample_lines.append(line.rstrip("\n"))
    except Exception:
        return None
    return _parse(sample_lines)


def _parse(sample_lines: list[str]) -> Optional[pd.DataFrame]:
    fields: dict[str, list[str]] = {}
    char_rows: list[list[str]] = []
    for line in sample_lines:
        parts = line.split("\t")
        key = parts[0]
        vals = [p.strip().strip('"') for p in parts[1:]]
        if key == "!Sample_characteristics_ch1":
            char_rows.append(vals)
        else:
            fields.setdefault(key, vals)              # first occurrence wins
    gsms = fields.get("!Sample_geo_accession", [])
    n = len(gsms)
    if n == 0:
        return None
    titles = _pad(fields.get("!Sample_title", []), n)
    sources = _pad(fields.get("!Sample_source_name_ch1", []), n)
    # each characteristics line is one attribute ("tissue: Tumor") across all samples; join per sample
    chars = ["" for _ in range(n)]
    for row in char_rows:
        for i, v in enumerate(_pad(row, n)):
            if v:
                chars[i] = f"{chars[i]}; {v}".lstrip("; ")
    df = pd.DataFrame({"title": titles, "source": sources, "characteristics": chars},
                      index=[g.strip() for g in gsms])
    df.index.name = "gsm"
    return df


def build_sample_meta(columns, meta: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Attach series-matrix metadata to the count columns, indexed by the count column name.

    Columns are matched to samples by GSM id or by (normalized) title; unmatched columns keep an
    empty blob. Returns None if nothing matched (caller keeps its filename-based meta).
    """
    cols = [str(c) for c in columns]
    by_gsm = {g.upper(): g for g in meta.index}
    by_title = {}
    for g, t in meta["title"].items():
        nt = _norm(t)
        if nt:
            by_title.setdefault(nt, g)

    matched: dict[str, Optional[str]] = {}
    for c in cols:
        gsm = None
        m = re.search(r"GSM\d+", c, re.IGNORECASE)
        if m and m.group(0).upper() in by_gsm:
            gsm = by_gsm[m.group(0).upper()]
        if gsm is None:
            nc = _norm(c)
            if nc in by_title:
                gsm = by_title[nc]
            elif nc:
                for nt, g in by_title.items():
                    if nc == nt or (len(nc) >= 4 and (nc in nt or nt in nc)):
                        gsm = g
                        break
        matched[c] = gsm

    if not any(matched.values()):
        return None
    rows = []
    for c in cols:
        g = matched[c]
        if g is not None:
            r = meta.loc[g]
            rows.append({"gsm": g, "title": r["title"], "source_name": r["source"],
                         "characteristics": r["characteristics"]})
        else:
            rows.append({"gsm": "", "title": "", "source_name": "", "characteristics": ""})
    return pd.DataFrame(rows, index=cols)


def _pad(vals: list[str], n: int) -> list[str]:
    return (list(vals) + [""] * n)[:n]


def _norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())
