"""One-time: cache the ARCHS4 per-sample index so remote counts work without the 58 GB file.

Thin wrapper over omics_oracle.data.build_remote_index (kept for `python scripts/...` use).
Configure with ARCHS4_URL and ARCHS4_INDEX_DIR; defaults target human gene-level counts.
"""
import os
import sys

from omics_oracle.data import build_remote_index


def main():
    index_dir = os.environ.get("ARCHS4_INDEX_DIR", "/home/omics-oracle/data/archs4/index")
    build_remote_index(index_dir=index_dir, url=os.environ.get("ARCHS4_URL"))
    print("DONE")


if __name__ == "__main__":
    sys.exit(main())
