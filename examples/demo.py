"""Minimal end-to-end demo of the stable API. Run: python examples/demo.py

Hits NCBI GEO live for real matches. Reanalysis is metadata-only unless ARCHS4_H5 is set.
"""
from omics_oracle import project_experiment

PROMPT = "RNA-seq of human colorectal cancer tumor vs adjacent normal tissue, Illumina, 20 per group"

if __name__ == "__main__":
    report = project_experiment(PROMPT, max_datasets=5)
    print(report.to_markdown())

    # The same report as a serializable dict — this is what another app would consume:
    d = report.to_dict()
    print("\n--- machine-readable summary ---")
    print("schema:", d["schema_version"], "| feasibility:", d["projection"]["feasibility"])
    for m in d["matches"]:
        print(f'  {m["accession"]}  sim={m["similarity"]}  blocking={m["blocking"]}')
