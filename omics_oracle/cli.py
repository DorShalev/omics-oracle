"""CLI:
  Project an experiment:   python -m omics_oracle "<prompt>" [--json] [--h5 PATH] ...
  Enable remote counts:    python -m omics_oracle --build-index [--species human]
  Download full H5:        python -m omics_oracle --download-h5 human --yes
"""
from __future__ import annotations

import argparse
import json
import sys


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="omics_oracle",
        description="Sanity-check a planned omics experiment against public data.")
    ap.add_argument("prompt", nargs="?", help="natural-language description of the planned experiment")
    ap.add_argument("--max-datasets", type=int, default=5, help="max matched datasets (default 5)")
    ap.add_argument("--reanalyze-top", type=int, default=3, help="how many top matches to reanalyze")
    ap.add_argument("--h5", default=None, help="path to a local ARCHS4 H5 for reprocessed counts")
    ap.add_argument("--model", default=None, help="LiteLLM/OpenRouter model id (else Claude CLI/heuristic)")
    ap.add_argument("--no-judge", action="store_true", help="disable the LLM relevance judge")
    ap.add_argument("--json", action="store_true", help="emit JSON instead of Markdown")
    # maintenance actions
    ap.add_argument("--build-index", action="store_true",
                    help="cache the ARCHS4 sample index to enable remote (no-58GB) counts")
    ap.add_argument("--download-h5", metavar="SPECIES", default=None,
                    help="download the full ARCHS4 gene-count H5 (human|mouse)")
    ap.add_argument("--h5-dir", default="data/archs4", help="dest dir for --download-h5")
    ap.add_argument("--connections", type=int, default=8, help="parallel segments for --download-h5")
    ap.add_argument("--yes", action="store_true", help="skip the download size confirmation")
    args = ap.parse_args(argv)

    if args.build_index:
        return _build_index(args)
    if args.download_h5:
        return _download(args)
    if not args.prompt:
        ap.error("a prompt is required (or use --build-index / --download-h5)")

    from .pipeline import project_experiment
    report = project_experiment(
        args.prompt, max_datasets=args.max_datasets, reanalyze_top=args.reanalyze_top,
        h5_path=args.h5, llm_model=args.model, judge=not args.no_judge,
    )
    if args.json:
        json.dump(report.to_dict(), sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        sys.stdout.write(report.to_markdown() + "\n")
    return 0


def _build_index(args) -> int:
    import os
    from .data import build_remote_index
    index_dir = os.environ.get("ARCHS4_INDEX_DIR", os.path.join(args.h5_dir, "index"))
    print("Building ARCHS4 remote index (one-time; reads ~1.1M sample records over HTTP)…")
    build_remote_index(index_dir=index_dir)
    print(f"index ready at {index_dir} — remote counts are now enabled.")
    return 0


def _download(args) -> int:
    from .data import download_archs4_h5, resolve_url, remote_size
    species = args.download_h5.lower()
    if not args.yes:
        url = resolve_url(species)
        size = remote_size(url)
        gb = f"{size/1e9:.1f} GB" if size else "unknown size"
        print(f"About to download the full ARCHS4 {species} H5 ({gb}) to {args.h5_dir}.")
        print("This is ONE large file. For most use you don't need it — the remote index streams")
        print("only the samples you query. Re-run with --yes to proceed.")
        return 2
    try:
        download_archs4_h5(species, dest_dir=args.h5_dir, connections=args.connections,
                           assume_yes=True)
    except Exception as e:
        print(f"download failed: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
