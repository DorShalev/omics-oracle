"""ARCHS4 counts WITHOUT downloading the 58 GB H5 — stream only the samples you need.

ARCHS4's `/data/expression` is chunked as (2000 genes x 1 sample), so a single sample is a
handful of small HDF5 chunks. This source reads just the columns of a matched GEO series over
HTTP range requests, using a one-time cached per-sample index (series_id / geo_accession /
gene symbols). Typical cost per dataset: tens of MB, not 58 GB.

The full-file local backend (`ARCHS4Source`) remains available for offline / heavy use; both
implement the same DatasetSource interface and return identical CountMatrix objects.
"""
from __future__ import annotations

import os
from typing import Optional

import numpy as np
import pandas as pd

from ..types import CountMatrix, DatasetMatch, ExperimentSpec

DEFAULT_URL = os.environ.get(
    "ARCHS4_URL", "https://s3.k8s.maayanlab.cloud/archs4/files/human_gene_v2.latest.h5")
DEFAULT_INDEX_DIR = os.environ.get("ARCHS4_INDEX_DIR", "/home/omics-oracle/data/archs4/index")

# Exclude likely single-cell samples from bulk DE (ARCHS4 provides a probability per sample).
SC_PROB_CUTOFF = 0.5
MAX_SAMPLES = int(os.environ.get("ARCHS4_MAX_SAMPLES", "200"))  # safety cap per dataset fetch


class ARCHS4RemoteSource:
    name = "ARCHS4"

    def __init__(self, url: Optional[str] = None, index_dir: Optional[str] = None,
                 block_size: int = 2 ** 22):
        self.url = url or DEFAULT_URL
        self.index_dir = index_dir or DEFAULT_INDEX_DIR
        self.block_size = block_size
        self._series = None
        self._geo = None
        self._genes = None
        self._scprob = None
        self._h5 = None

    # -- index (cached) -------------------------------------------------------------------

    def _index_ready(self) -> bool:
        return all(os.path.exists(os.path.join(self.index_dir, f))
                   for f in ("series_id.npy", "geo_accession.npy", "genes.npy"))

    def _load_index(self):
        if self._series is None:
            self._series = np.load(os.path.join(self.index_dir, "series_id.npy"), allow_pickle=True)
            self._geo = np.load(os.path.join(self.index_dir, "geo_accession.npy"), allow_pickle=True)
            self._genes = np.load(os.path.join(self.index_dir, "genes.npy"), allow_pickle=True)
            scp = os.path.join(self.index_dir, "singlecellprobability.npy")
            self._scprob = np.load(scp, allow_pickle=True) if os.path.exists(scp) else None

    def _sample_rows(self, accession: str) -> np.ndarray:
        """Row indices of samples belonging to a GEO series (series_id may list several)."""
        self._load_index()
        acc = accession.strip()

        def member(s):
            s = s if isinstance(s, str) else str(s)
            # series_id entries can be joined by comma/semicolon/tab
            return acc in [p.strip() for p in s.replace(";", ",").replace("\t", ",").split(",")]

        mask = np.array([member(s) for s in self._series], dtype=bool)
        idx = np.nonzero(mask)[0]
        if self._scprob is not None and len(idx):
            # drop single-cell samples for bulk RNA-seq DE
            idx = idx[np.array([self._scprob[i] <= SC_PROB_CUTOFF for i in idx])]
        return idx

    # -- H5 handle ------------------------------------------------------------------------

    def _open(self):
        if self._h5 is None:
            import h5py
            import fsspec
            f = fsspec.open(self.url, block_size=self.block_size).open()
            self._h5 = h5py.File(f, "r")
        return self._h5

    # -- interface ------------------------------------------------------------------------

    def search(self, spec: ExperimentSpec, limit: int = 10) -> list[DatasetMatch]:
        # GEO is the finder; this source only supplies counts.
        return []

    def supports_counts(self) -> bool:
        return self._index_ready()

    def fetch_counts(self, match: DatasetMatch) -> Optional[CountMatrix]:
        if not self._index_ready():
            return None
        try:
            idx = self._sample_rows(match.accession)
        except Exception:
            return None
        if idx.size == 0:
            return None
        if idx.size > MAX_SAMPLES:
            idx = idx[:MAX_SAMPLES]
        idx = np.sort(idx)
        try:
            h5 = self._open()
            # expression is (genes, samples); pull just our sample columns
            expr = h5["data/expression"][:, idx]            # genes x n_idx, uint32
            gsm = [self._decode(self._geo[i]) for i in idx]
            genes = [self._decode(g) for g in self._genes]
            meta = self._sample_meta(h5, idx, gsm)
        except Exception:
            return None
        df = pd.DataFrame(np.asarray(expr), index=genes, columns=gsm)
        df = df[~df.index.duplicated(keep="first")]
        return CountMatrix(
            accession=match.accession, genes=list(df.index), samples=gsm, counts=df,
            sample_meta=meta,
            provenance=f"ARCHS4 v2 remote stream ({len(gsm)} samples, {os.path.basename(self.url)})",
        )

    def _sample_meta(self, h5, idx, gsm) -> pd.DataFrame:
        fields = ["characteristics_ch1", "source_name_ch1", "title", "instrument_model"]
        data = {}
        for fld in fields:
            try:
                col = h5[f"meta/samples/{fld}"][:][idx]
                data[fld] = [self._decode(v) for v in col]
            except Exception:
                continue
        return pd.DataFrame(data, index=gsm)

    @staticmethod
    def _decode(x):
        return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)
