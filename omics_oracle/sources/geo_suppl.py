"""Per-study counts from GEO supplementary files — the default counts backend.

On request, for a matched GEO series, this lists the study's own supplementary files, picks the
expression table (raw counts preferred over normalized), downloads only that file (a few MB),
parses it, detects whether the values are raw counts or normalized, and returns a CountMatrix.
No bulk database is ever downloaded.

Honesty: it never fabricates. If no usable table exists it returns None with a recorded reason;
if the table is normalized it still returns the matrix (tagged data_type="normalized") so the
reanalysis layer can report it and refuse DE rather than run pydeseq2 on invalid input.
"""
from __future__ import annotations

import os
import re
from typing import Optional

import pandas as pd
import requests
from requests.adapters import HTTPAdapter

try:
    from urllib3.util.retry import Retry
except Exception:  # pragma: no cover
    Retry = None

from ..types import CountMatrix, DatasetMatch, ExperimentSpec
from ..tables import (read_expression_table, classify_data_type, is_probably_expression_file,
                      filename_score)

_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
_HREF = re.compile(r'href="([^"?]+)"')
# tokens that distinguish per-condition split files of the same matrix
_SPLIT_TOKENS = re.compile(
    r"(ctrl|control|exp|experiment|treat(?:ed|ment)?|case|tumou?r|normal|wt|ko|kd|"
    r"mut|wildtype|disease|healthy)", re.IGNORECASE)


class GeoSupplSource:
    name = "GEO-suppl"

    def __init__(self, cache_dir: str = "data/geo_suppl",
                 max_bytes: Optional[int] = None, timeout: float = 60.0):
        self.cache_dir = cache_dir
        self.max_bytes = max_bytes or int(os.environ.get("GEO_SUPPL_MAX_BYTES", str(50 * 1024 ** 2)))
        self.timeout = timeout
        self.last_warnings: list[str] = []
        self.session = requests.Session()
        if Retry is not None:
            retry = Retry(total=3, backoff_factor=0.5,
                          status_forcelist=(429, 500, 502, 503, 504),
                          allowed_methods=("GET", "HEAD"))
            self.session.mount("https://", HTTPAdapter(max_retries=retry))

    # -- interface ------------------------------------------------------------------------

    def search(self, spec: ExperimentSpec, limit: int = 10) -> list[DatasetMatch]:
        return []  # GEO E-utilities is the finder; this only supplies counts

    def supports_counts(self) -> bool:
        return True  # works on-demand per accession; no index/download prerequisite

    def fetch_counts(self, match: DatasetMatch) -> Optional[CountMatrix]:
        self.last_warnings = []
        acc = match.accession
        try:
            files = self._list_suppl(acc)
        except Exception as e:
            self._warn(f"{acc}: could not list supplementary files ({type(e).__name__})")
            return None
        candidates = [(f, u) for f, u in files if is_probably_expression_file(f)]
        if not candidates:
            self._warn(f"{acc}: no supplementary expression table found "
                       f"(files: {', '.join(f for f, _ in files[:6]) or 'none'})")
            return None

        group = self._best_group(candidates)
        frames, used = [], []
        for fname, url in group:
            raw = self._download(acc, fname, url)
            if raw is None:
                continue
            try:
                df = read_expression_table(raw, fname)
            except ValueError as e:
                self._warn(f"{acc}: could not parse {fname} ({e})")
                continue
            df.columns = [f"{c}" for c in df.columns]
            frames.append((fname, df))
            used.append(fname)
        if not frames:
            self._warn(f"{acc}: found candidate table(s) but none parsed")
            return None

        counts, meta, dtype = self._assemble(frames)
        if counts is None:
            self._warn(f"{acc}: candidate tables did not share a gene index; cannot merge")
            return None
        if dtype == "normalized":
            self._warn(f"{acc}: supplementary table is normalized ({used[0]}); "
                       "DE requires raw counts and will be skipped (not faked)")
        elif dtype == "unknown":
            self._warn(f"{acc}: could not confirm the table is raw counts ({used[0]}); DE skipped")

        return CountMatrix(
            accession=acc, genes=list(counts.index), samples=list(counts.columns),
            counts=counts, sample_meta=meta, data_type=dtype, source_file=", ".join(used),
            source_db="GEO", tier="primary",
            provenance=f"GEO supplementary file(s): {', '.join(used)}",
        )

    # -- discovery ------------------------------------------------------------------------

    def _suppl_url(self, accession: str) -> str:
        num = accession[3:]
        prefix = num[:-3] if len(num) > 3 else ""
        return f"{_BASE}/GSE{prefix}nnn/{accession}/suppl/"

    def _list_suppl(self, accession: str) -> list[tuple[str, str]]:
        url = self._suppl_url(accession)
        r = self.session.get(url, timeout=self.timeout)
        r.raise_for_status()
        out: list[tuple[str, str]] = []
        for href in _HREF.findall(r.text):
            name = href.split("/")[-1]
            if not name or name in ("..",) or href.startswith("/") and name == accession:
                continue
            if name.startswith(accession) or name.lower().endswith(
                    (".gz", ".csv", ".txt", ".tsv", ".xlsx")):
                out.append((name, url + name))
        # de-dup preserving order
        seen, uniq = set(), []
        for n, u in out:
            if n not in seen:
                seen.add(n); uniq.append((n, u))
        return uniq

    def _best_group(self, candidates: list[tuple[str, str]]) -> list[tuple[str, str]]:
        """Pick the best single table, or a set of per-condition split files to merge."""
        # signature = filename with split/replicate tokens removed => split files collapse together
        def sig(fname: str) -> str:
            s = _SPLIT_TOKENS.sub("", fname.lower())
            s = re.sub(r"[_\-.]?(rep)?\d+", "", s)
            return s
        groups: dict[str, list[tuple[str, str]]] = {}
        for f, u in candidates:
            groups.setdefault(sig(f), []).append((f, u))
        # score each group by its best filename_score; prefer larger groups on ties (splits)
        best = max(groups.values(),
                   key=lambda g: (max(filename_score(f) for f, _ in g), len(g)))
        best.sort(key=lambda fu: filename_score(fu[0]), reverse=True)
        # a genuine split set (>1 file) is merged; otherwise just the single best file
        return best if len(best) > 1 else best[:1]

    # -- download (capped, cached) --------------------------------------------------------

    def _download(self, accession: str, fname: str, url: str) -> Optional[bytes]:
        cdir = os.path.join(self.cache_dir, accession)
        cpath = os.path.join(cdir, fname)
        if os.path.exists(cpath) and os.path.getsize(cpath) > 0:
            with open(cpath, "rb") as fh:
                return fh.read()
        try:
            head = self.session.head(url, timeout=self.timeout, allow_redirects=True)
            clen = int(head.headers.get("content-length", "0") or 0)
            if clen and clen > self.max_bytes:
                self._warn(f"{accession}: {fname} is {clen/1e6:.0f} MB > cap "
                           f"{self.max_bytes/1e6:.0f} MB; skipped")
                return None
        except Exception:
            pass
        try:
            buf = bytearray()
            with self.session.get(url, stream=True, timeout=self.timeout) as r:
                r.raise_for_status()
                for chunk in r.iter_content(chunk_size=1 << 20):
                    buf.extend(chunk)
                    if len(buf) > self.max_bytes:
                        self._warn(f"{accession}: {fname} exceeded {self.max_bytes/1e6:.0f} MB "
                                   "while streaming; skipped")
                        return None
        except Exception as e:
            self._warn(f"{accession}: download failed for {fname} ({type(e).__name__})")
            return None
        os.makedirs(cdir, exist_ok=True)
        with open(cpath, "wb") as fh:
            fh.write(buf)
        return bytes(buf)

    # -- assembly -------------------------------------------------------------------------

    def _assemble(self, frames: list[tuple[str, pd.DataFrame]]):
        if len(frames) == 1:
            fname, df = frames[0]
            meta = self._sample_meta(df.columns, fname)
            return df, meta, classify_data_type(df, fname)
        # merge per-condition split files on the shared gene index
        merged = None
        origin: dict[str, str] = {}
        for fname, df in frames:
            df = df.loc[:, ~df.columns.duplicated(keep="first")]
            for c in df.columns:
                origin[str(c)] = fname
            merged = df if merged is None else merged.join(df, how="inner")
        if merged is None or merged.shape[1] < 2 or merged.shape[0] == 0:
            return None, None, "unknown"
        meta = pd.DataFrame(
            {"column_name": [str(c) for c in merged.columns],
             "source_file": [origin.get(str(c), "") for c in merged.columns]},
            index=[str(c) for c in merged.columns])
        # type from the first (highest-scored) file's name
        return merged, meta, classify_data_type(merged, frames[0][0])

    def _sample_meta(self, columns, fname: str) -> pd.DataFrame:
        cols = [str(c) for c in columns]
        return pd.DataFrame({"column_name": cols, "source_file": [fname] * len(cols)}, index=cols)

    def _warn(self, msg: str) -> None:
        self.last_warnings.append(msg)
