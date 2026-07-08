"""Minimal CLI:  python -m omics_oracle "<experiment prompt>" [--json] [--h5 PATH] ..."""
from __future__ import annotations

import argparse
import json
import sys

from .pipeline import project_experiment


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="omics_oracle",
        description="Sanity-check a planned omics experiment against public data.")
    ap.add_argument("prompt", help="natural-language description of the planned experiment")
    ap.add_argument("--max-datasets", type=int, default=5, help="max matched datasets (default 5)")
    ap.add_argument("--reanalyze-top", type=int, default=3, help="how many top matches to reanalyze")
    ap.add_argument("--h5", default=None, help="path to a local ARCHS4 H5 for reprocessed counts")
    ap.add_argument("--model", default=None, help="LiteLLM/OpenRouter model id (else heuristic)")
    ap.add_argument("--json", action="store_true", help="emit JSON instead of Markdown")
    args = ap.parse_args(argv)

    report = project_experiment(
        args.prompt, max_datasets=args.max_datasets, reanalyze_top=args.reanalyze_top,
        h5_path=args.h5, llm_model=args.model,
    )
    if args.json:
        json.dump(report.to_dict(), sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        sys.stdout.write(report.to_markdown() + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
