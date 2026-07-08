"""ARCHS4 H5 download helper — resumable, parallel, size-disclosing.

The full ARCHS4 counts matrix is a single large H5 (~58 GB human gene-level), not many files.
You usually do NOT need it: the remote streaming source (`sources/archs4_remote.py`) reads only
the samples you query. Download the full file only for offline or heavy repeated use.

This helper streams the file in parallel byte-range segments (the server supports range GETs),
resumes a partial download, and refuses to start if disk is short.
"""
from __future__ import annotations

import os
import threading
from typing import Optional

import requests

# species -> canonical remote object (via archs4py config, with a static fallback)
_FALLBACK_URLS = {
    "human": "https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.latest.h5",
    "mouse": "https://s3.dev.maayanlab.cloud/archs4/files/mouse_gene_v2.latest.h5",
}


def build_remote_index(index_dir: str = "data/archs4/index", url: Optional[str] = None,
                       on_log=print) -> str:
    """Cache the ARCHS4 per-sample index (series_id / geo_accession / gene symbols) locally.

    One-time, slow (reads ~1.1M vlen-string records over HTTP). Enables the remote counts
    source to map GEO series -> sample columns without the 58 GB download.
    """
    import time
    import h5py
    import fsspec
    import numpy as np

    url = url or DEFAULT_URL
    os.makedirs(index_dir, exist_ok=True)
    f = fsspec.open(url, block_size=2 ** 22).open()
    h5 = h5py.File(f, "r")
    on_log(f"ARCHS4 version {h5['meta/info/version'][()]} — caching sample index to {index_dir}")

    def dec(a):
        return np.array([x.decode() if isinstance(x, (bytes, bytearray)) else str(x) for x in a])

    for name in ("series_id", "geo_accession", "singlecellprobability"):
        t = time.time()
        arr = h5[f"meta/samples/{name}"][:]
        if arr.dtype == object:
            arr = dec(arr)
        np.save(os.path.join(index_dir, f"{name}.npy"), arr)
        on_log(f"  cached {name} ({len(arr)}) in {time.time()-t:.1f}s")
    np.save(os.path.join(index_dir, "genes.npy"), dec(h5["meta/genes/symbol"][:]))
    on_log("  cached gene symbols")
    return index_dir


DEFAULT_URL = os.environ.get(
    "ARCHS4_URL", "https://s3.k8s.maayanlab.cloud/archs4/files/human_gene_v2.latest.h5")


def resolve_url(species: str) -> str:
    species = species.lower()
    try:
        import archs4py
        conf = archs4py.utils.get_config()
        return conf["GENE_COUNTS"][species.upper()]["latest"]["primary"]
    except Exception:
        if species not in _FALLBACK_URLS:
            raise ValueError(f"unknown species '{species}' (use 'human' or 'mouse')")
        return _FALLBACK_URLS[species]


def remote_size(url: str, timeout: float = 30.0) -> Optional[int]:
    r = requests.head(url, allow_redirects=True, timeout=timeout)
    n = r.headers.get("content-length")
    return int(n) if n else None


def _free_bytes(path: str) -> int:
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize


def download_archs4_h5(species: str = "human", dest_dir: str = "data/archs4",
                       version: str = "latest", connections: int = 8,
                       assume_yes: bool = False, on_log=print) -> str:
    """Download the ARCHS4 gene-count H5. Returns the local path.

    Safeguards: discloses the size and requires confirmation (unless assume_yes), checks free
    disk (needs size + 5 GB headroom), and resumes an existing partial file.
    """
    url = resolve_url(species) if version == "latest" else resolve_url(species)
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, os.path.basename(url))

    size = remote_size(url)
    size_gb = (size / 1e9) if size else None
    on_log(f"ARCHS4 {species} gene counts")
    on_log(f"  url : {url}")
    on_log(f"  size: {size_gb:.1f} GB" if size_gb else "  size: unknown")
    on_log(f"  dest: {dest}")

    have = os.path.getsize(dest) if os.path.exists(dest) else 0
    if size and have == size:
        on_log("  already fully downloaded.")
        return dest
    if have:
        on_log(f"  resuming from {have / 1e9:.1f} GB already on disk.")

    if size:
        free = _free_bytes(dest_dir)
        need = size - have + 5 * 1024 ** 3
        if free < need:
            raise RuntimeError(f"insufficient disk: need ~{need/1e9:.0f} GB free, have "
                               f"{free/1e9:.0f} GB in {dest_dir}")

    if not assume_yes:
        # This helper is normally called from the CLI, which handles interactive confirmation.
        raise RuntimeError("refusing to download without confirmation; pass assume_yes=True "
                           "(CLI: add --yes) after reviewing the size above")

    if not size or connections <= 1 or not _supports_range(url):
        _download_stream(url, dest, have, size, on_log)
    else:
        _download_parallel_resumable(url, dest, size, connections, on_log)
    final = os.path.getsize(dest)
    if size and final != size:
        raise RuntimeError(f"size mismatch after download: got {final}, expected {size}")
    on_log(f"done: {dest} ({final/1e9:.1f} GB)")
    return dest


def _supports_range(url: str) -> bool:
    try:
        r = requests.get(url, headers={"Range": "bytes=0-0"}, stream=True, timeout=30,
                         allow_redirects=True)
        return r.status_code == 206
    except Exception:
        return False


def _download_stream(url: str, dest: str, have: int, size: Optional[int], on_log) -> None:
    headers = {"Range": f"bytes={have}-"} if have else {}
    with requests.get(url, headers=headers, stream=True, timeout=60, allow_redirects=True) as r:
        r.raise_for_status()
        mode = "ab" if have else "wb"
        done = have
        with open(dest, mode) as fh:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                fh.write(chunk)
                done += len(chunk)
                if size and done % (512 * 1024 * 1024) < 1024 * 1024:
                    on_log(f"  {done/1e9:.1f}/{size/1e9:.1f} GB ({100*done/size:.0f}%)")


def _download_parallel_resumable(url: str, dest: str, size: int, connections: int,
                                 on_log) -> None:
    """Download `size` bytes in N segments, each to its own resumable .partK file.

    A segment resumes from whatever it already has (Range from its current length), so an
    interrupted run continues instead of restarting. Segments are concatenated into `dest`
    once all are complete.
    """
    seg = size // connections
    bounds = [(i * seg, (size - 1 if i == connections - 1 else (i + 1) * seg - 1))
              for i in range(connections)]
    parts = [f"{dest}.part{i}" for i in range(connections)]
    errors: list[BaseException] = []
    lock = threading.Lock()
    total_have = sum(os.path.getsize(p) for p in parts if os.path.exists(p))
    on_log(f"parallel: {connections} segments"
           + (f", resuming ({total_have/1e9:.1f} GB already present)" if total_have else ""))

    def worker(i, lo, hi):
        part = parts[i]
        target = hi - lo + 1
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have >= target:
            return
        try:
            with requests.get(url, headers={"Range": f"bytes={lo+have}-{hi}"}, stream=True,
                              timeout=120, allow_redirects=True) as r:
                r.raise_for_status()
                with open(part, "ab") as fh:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        fh.write(chunk)
        except BaseException as e:  # noqa: BLE001
            with lock:
                errors.append(e)

    threads = [threading.Thread(target=worker, args=(i, lo, hi), daemon=True)
               for i, (lo, hi) in enumerate(bounds)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    if errors:
        raise RuntimeError(f"parallel download failed ({len(errors)} segment error(s)): "
                           f"{errors[0]}; re-run to resume the incomplete segments")
    # verify each segment, then concatenate
    for i, (lo, hi) in enumerate(bounds):
        got = os.path.getsize(parts[i])
        if got != hi - lo + 1:
            raise RuntimeError(f"segment {i} incomplete ({got}/{hi-lo+1}); re-run to resume")
    with open(dest, "wb") as out:
        for p in parts:
            with open(p, "rb") as pf:
                while True:
                    buf = pf.read(8 * 1024 * 1024)
                    if not buf:
                        break
                    out.write(buf)
    for p in parts:
        os.remove(p)
