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
                      is_single_cell_file, filename_score)

_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
_HREF = re.compile(r'href="([^"?]+)"')
# tokens that distinguish per-condition split files of the same matrix
_SPLIT_TOKENS = re.compile(
    r"(ctrl|control|exp|experiment|treat(?:ed|ment)?|case|tumou?r|normal|wt|ko|kd|"
    r"mut|wildtype|disease|healthy)", re.IGNORECASE)
# assay-type tokens stripped from per-sample filenames so a sample's count + FPKM files unify
_ASSAY_TOKENS = re.compile(
    r"[_\-.]?(raw[_\-.]?counts?|read[_\-.]?counts?|feature[_\-.]?counts?|expected[_\-.]?counts?|"
    r"gene[_\-.]?counts?|counts?|fpkm|rpkm|tpm|cpm|rpm|normali\w*|vst|rlog)", re.IGNORECASE)


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
            # Distinguish "single-cell (set aside)" from "genuinely no data". Recent GEO top-hits
            # are dominated by scRNA-seq, deposited as 10x .mtx/.h5 or Seurat .rds — cell x gene,
            # not a bulk table. A lone GSE..._RAW.tar usually IS the data; peek inside to see which.
            tar = next((fu for fu in files if fu[0].lower().endswith("_raw.tar")), None)
            if tar is not None:
                kind = self._classify_tar(acc, tar[0], tar[1])
                if kind == "single-cell":
                    match.modality = "single-cell"
                    self._warn(f"{acc}: single-cell RNA-seq (10x matrix / .h5 / .rds in "
                               "RAW.tar) — set aside; bulk DE not applicable (single-cell "
                               "analysis is planned)")
                    return None
                cm = self._assemble_from_tar(acc, tar[0], tar[1], match)
                if cm is not None:
                    return cm
            if any(is_single_cell_file(f) for f, _ in files):
                match.modality = "single-cell"
                self._warn(f"{acc}: single-cell RNA-seq deposition (.mtx/.h5/.rds) — set aside; "
                           "bulk DE not applicable (single-cell analysis is planned)")
                return None
            self._warn(f"{acc}: no supplementary expression table found "
                       f"(files: {', '.join(f for f, _ in files[:6]) or 'none'})")
            return None
        match.modality = "bulk"

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
            self._warn(f"{acc}: supplementary table is normalized ({used[0]}); DE will use a "
                       "Welch t-test on log values (not DESeq2, which needs raw counts)")
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

    def _download(self, accession: str, fname: str, url: str, cap: Optional[int] = None,
                  silent: bool = False, prefix_only: bool = False) -> Optional[bytes]:
        """Fetch a supplementary file (capped, cached). ``cap`` overrides the default size cap
        (RAW.tars need a larger one); ``prefix_only`` returns just the first ``cap`` bytes without
        caching (used to peek a tar's headers); ``silent`` suppresses size/skip warnings."""
        limit = cap or self.max_bytes
        cdir = os.path.join(self.cache_dir, accession)
        cpath = os.path.join(cdir, fname)
        if not prefix_only and os.path.exists(cpath) and os.path.getsize(cpath) > 0:
            with open(cpath, "rb") as fh:
                return fh.read()
        if not prefix_only:
            try:
                head = self.session.head(url, timeout=self.timeout, allow_redirects=True)
                clen = int(head.headers.get("content-length", "0") or 0)
                if clen and clen > limit:
                    if not silent:
                        self._warn(f"{accession}: {fname} is {clen/1e6:.0f} MB > cap "
                                   f"{limit/1e6:.0f} MB; skipped")
                    return None
            except Exception:
                pass
        try:
            buf = bytearray()
            with self.session.get(url, stream=True, timeout=self.timeout) as r:
                r.raise_for_status()
                for chunk in r.iter_content(chunk_size=1 << 20):
                    buf.extend(chunk)
                    if prefix_only and len(buf) >= limit:
                        return bytes(buf[:limit])          # enough to read tar headers
                    if not prefix_only and len(buf) > limit:
                        if not silent:
                            self._warn(f"{accession}: {fname} exceeded {limit/1e6:.0f} MB "
                                       "while streaming; skipped")
                        return None
        except Exception as e:
            if not silent:
                self._warn(f"{accession}: download failed for {fname} ({type(e).__name__})")
            return None
        if prefix_only:
            return bytes(buf)
        os.makedirs(cdir, exist_ok=True)
        with open(cpath, "wb") as fh:
            fh.write(buf)
        return bytes(buf)

    # -- RAW.tar (per-sample bulk count files bundled into one archive) --------------------

    def _classify_tar(self, acc: str, name: str, url: str) -> str:
        """Peek at a RAW.tar's first members to tell single-cell from bulk without full download.

        Tar stores each member as a 512-byte header (name at offset 0) followed by its data; a
        small prefix is enough to read the first member name(s), which reveal the modality —
        `..._matrix.mtx.gz`/`..._bc_matrix.h5` (single-cell) vs `..._counts.txt.gz` (bulk).
        """
        import tarfile, io as _io
        prefix = self._download(acc, name, url, cap=2 * 1024 ** 2, silent=True, prefix_only=True)
        if not prefix:
            return "unknown"
        names: list[str] = []
        try:
            with tarfile.open(fileobj=_io.BytesIO(prefix), mode="r:") as tf:
                for m in tf:                          # stops when the prefix runs out
                    if m.isfile():
                        names.append(m.name.split("/")[-1])
                    if len(names) >= 6:
                        break
        except Exception:
            pass
        if not names:                                 # couldn't parse headers from the prefix
            return "unknown"
        if any(is_single_cell_file(n) for n in names):
            return "single-cell"
        return "bulk"

    def _assemble_from_tar(self, acc: str, name: str, url: str,
                           match: DatasetMatch) -> Optional[CountMatrix]:
        """Download a bulk RAW.tar and stitch its per-sample count files into one matrix."""
        import tarfile, io as _io
        raw = self._download(acc, name, url, cap=int(os.environ.get(
            "GEO_TAR_MAX_BYTES", str(400 * 1024 ** 2))))
        if raw is None:
            return None
        # each sample may appear as several files (e.g. read_count.xlsx AND FPKM.xlsx); keep the
        # highest-scored file per sample id (raw counts beat normalized) so we don't double-count.
        best: dict[str, tuple[int, pd.Series]] = {}
        try:
            with tarfile.open(fileobj=_io.BytesIO(raw), mode="r:") as tf:
                for m in tf:
                    base = m.name.split("/")[-1]
                    if not m.isfile() or is_single_cell_file(base) or not is_probably_expression_file(base):
                        continue
                    fh = tf.extractfile(m)
                    if fh is None:
                        continue
                    try:
                        s = self._read_one_sample(fh.read(), base)
                    except Exception:
                        continue                       # one unreadable member must not abort the tar
                    if s is None or len(s) <= 100:
                        continue
                    sid, score = self._sample_id(base), filename_score(base)
                    if sid not in best or score > best[sid][0]:
                        best[sid] = (score, s)
        except Exception as e:
            self._warn(f"{acc}: could not read {name} ({type(e).__name__})")
            return None
        series = {sid: v[1] for sid, v in best.items()}
        if len(series) < 2:
            self._warn(f"{acc}: {name} had no per-sample count tables to assemble")
            return None
        counts = pd.DataFrame(series).dropna(how="all")
        counts = counts[counts.sum(axis=1) > 0]
        if counts.shape[0] < 100 or counts.shape[1] < 2:
            return None
        match.modality = "bulk"
        dtype = classify_data_type(counts, name)
        self._warn(f"{acc}: assembled a {counts.shape[0]}×{counts.shape[1]} matrix from "
                   f"{counts.shape[1]} per-sample files inside {name} ({dtype})")
        meta = self._sample_meta(counts.columns, name)
        return CountMatrix(accession=acc, genes=list(counts.index), samples=list(counts.columns),
                           counts=counts, sample_meta=meta, data_type=dtype, source_file=name,
                           source_db="GEO", tier="primary",
                           provenance=f"assembled from per-sample files in {name}")

    def _read_one_sample(self, raw: bytes, fname: str) -> Optional[pd.Series]:
        """Parse one per-sample count file to a gene->value Series (take the count column)."""
        try:
            df = read_expression_table(raw, fname, min_cols=1)
        except ValueError:
            return None
        if df.shape[1] == 0:
            return None
        # featureCounts-style files carry Chr/Start/End/Strand/Length before the count column;
        # the sample count is the LAST numeric column. Single-column files just use that column.
        col = df.columns[-1]
        s = pd.to_numeric(df[col], errors="coerce").dropna()
        return s if len(s) else None

    @staticmethod
    def _sample_id(fname: str) -> str:
        """Sample id from a per-sample filename, with the assay-type token stripped so the same
        sample's count file and FPKM file collapse to ONE id (GSM..._Veh1_read_count +
        GSM..._Veh1_FPKM -> GSM..._Veh1)."""
        base = re.sub(r"\.(txt|tsv|csv|tab|xlsx)(\.gz)?$", "", fname, flags=re.IGNORECASE)
        base = _ASSAY_TOKENS.sub("", base)
        base = re.sub(r"[_\-.]{2,}", "_", base).strip("_-.")
        if base:
            return base
        m = re.match(r"(GSM\d+)", fname)
        return m.group(1) if m else fname

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
