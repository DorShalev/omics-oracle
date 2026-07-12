"""Bulk-coverage: .tab tables, single-cell set-aside, and RAW.tar per-sample assembly."""
import gzip
import io
import tarfile

import numpy as np
import pandas as pd

from omics_oracle.tables import is_probably_expression_file, is_single_cell_file
from omics_oracle.sources.geo_suppl import GeoSupplSource
from omics_oracle.types import DatasetMatch


def test_tab_files_are_expression_tables():
    assert is_probably_expression_file("GSE1_liver_FPKM.tab.gz")
    assert is_probably_expression_file("GSE1_counts.tab")


def test_single_cell_signatures_detected_and_excluded():
    for f in ["GSM1_matrix.mtx.gz", "GSM1_barcodes.tsv.gz", "GSM1_features.tsv.gz",
              "GSE1_filtered_feature_bc_matrix.h5", "GSE1_seurat.rds", "GSE1.h5ad"]:
        assert is_single_cell_file(f), f
        assert not is_probably_expression_file(f), f   # never fed to bulk DE


def _tar_bytes(members: dict[str, bytes]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:") as tf:
        for name, data in members.items():
            info = tarfile.TarInfo(name=name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
    return buf.getvalue()


def _per_sample_counts(gsm, label, rng, kind="read_count", depth=1.0):
    # realistic raw counts: overdispersed with a per-sample library depth so column sums vary
    # (uniform, equal-sum columns would look like TPM to the raw/normalized classifier)
    genes = [f"G{i}" for i in range(200)]
    vals = (rng.gamma(0.4, 300, size=200) * depth).round().astype(int)
    body = "gene\tcount\n" + "\n".join(f"{g}\t{v}" for g, v in zip(genes, vals))
    return f"GSM{gsm}_{label}_{kind}.txt", gzip.compress(body.encode())


def test_bulk_raw_tar_assembled(monkeypatch):
    rng = np.random.default_rng(0)
    members = {}
    for gsm, label, depth in [(1, "Veh1", 1.0), (2, "Veh2", 1.6), (3, "KO1", 0.6), (4, "KO2", 2.2)]:
        name, data = _per_sample_counts(gsm, label, rng, depth=depth)
        members[name] = data
    tar = _tar_bytes(members)
    s = GeoSupplSource()
    monkeypatch.setattr(s, "_list_suppl",
                        lambda acc: [("GSE9_RAW.tar", "http://x/raw"), ("filelist.txt", "http://x/fl")])
    monkeypatch.setattr(s, "_download",
                        lambda acc, fn, url, **kw: tar if fn.endswith("_RAW.tar") else None)
    m = DatasetMatch(accession="GSE9", source="GEO")
    cm = s.fetch_counts(m)
    assert cm is not None and m.modality == "bulk"
    assert cm.counts.shape == (200, 4)          # 4 per-sample files -> 4 columns
    assert cm.data_type == "raw"                 # integer counts detected as raw
    assert any("assembled" in w for w in s.last_warnings)


def test_single_cell_raw_tar_set_aside(monkeypatch):
    # a 10x triplet inside the tar -> classified single-cell, not analyzed
    members = {"GSM1_barcodes.tsv.gz": gzip.compress(b"AAA\nCCC\n"),
               "GSM1_features.tsv.gz": gzip.compress(b"G1\tG1\nG2\tG2\n"),
               "GSM1_matrix.mtx.gz": gzip.compress(b"%%MatrixMarket\n2 2 2\n1 1 5\n2 2 7\n")}
    tar = _tar_bytes(members)
    s = GeoSupplSource()
    monkeypatch.setattr(s, "_list_suppl", lambda acc: [("GSE10_RAW.tar", "http://x/raw")])
    monkeypatch.setattr(s, "_download", lambda acc, fn, url, **kw: tar)
    m = DatasetMatch(accession="GSE10", source="GEO")
    cm = s.fetch_counts(m)
    assert cm is None and m.modality == "single-cell"
    assert any("single-cell" in w for w in s.last_warnings)


def test_per_sample_dedup_prefers_raw_over_fpkm(monkeypatch):
    # same sample deposited twice (read_count + FPKM); assembler keeps ONE column, the raw one
    rng = np.random.default_rng(1)
    members = {}
    for gsm, label in [(1, "A"), (2, "B")]:
        n_cnt, d_cnt = _per_sample_counts(gsm, label, rng, kind="read_count")
        n_fpkm, d_fpkm = _per_sample_counts(gsm, label, rng, kind="FPKM")
        members[n_cnt] = d_cnt
        members[n_fpkm] = d_fpkm
    tar = _tar_bytes(members)
    s = GeoSupplSource()
    monkeypatch.setattr(s, "_list_suppl", lambda acc: [("GSE11_RAW.tar", "http://x/raw")])
    monkeypatch.setattr(s, "_download", lambda acc, fn, url, **kw: tar)
    cm = s.fetch_counts(DatasetMatch(accession="GSE11", source="GEO"))
    assert cm is not None
    assert cm.counts.shape[1] == 2               # 2 samples, not 4 (FPKM+count collapsed)
    assert set(cm.samples) == {"GSM1_A", "GSM2_B"}
