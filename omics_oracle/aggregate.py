"""Cross-experiment aggregator: per-gene concordance across datasets.

Concordance = how consistently a gene moves in the same direction across the datasets where
it was testable. We report direction agreement (fraction of datasets sharing the majority
sign), mean effect size, and the per-dataset context so a caller can always trace a claim.

Effect sizes are only comparable because ARCHS4 counts are uniformly reprocessed and each
dataset's DE uses the same pipeline; we still restrict to genes tested in >=2 datasets.
"""
from __future__ import annotations

from typing import Iterable

import numpy as np

from .types import Concordance, DEResult, GeneConcordance


def aggregate(de_results: Iterable[DEResult], min_datasets: int = 2,
              padj: float = 0.1, lfc: float = 0.5, top: int = 200) -> Concordance:
    ok = [d for d in de_results if d.status == "ok" and d.table is not None and len(d.table)]
    if len(ok) < min_datasets:
        return Concordance(n_datasets=len(ok),
                           note=f"only {len(ok)} dataset(s) yielded DE tables; need >= {min_datasets} "
                                "for cross-experiment concordance")

    # collect per-gene log2fc across datasets, keeping only clearly-moving genes per dataset
    per_gene: dict[str, dict[str, float]] = {}
    for d in ok:
        t = d.table
        sig = t[(t["padj"].fillna(1.0) < padj) & (t["log2fc"].abs() >= lfc)]
        for gene, l2fc in zip(sig["gene"], sig["log2fc"]):
            per_gene.setdefault(str(gene), {})[d.accession] = float(l2fc)

    genes: list[GeneConcordance] = []
    for gene, byds in per_gene.items():
        if len(byds) < min_datasets:
            continue
        vals = np.array(list(byds.values()), dtype=float)
        n_up = int((vals > 0).sum())
        n_down = int((vals < 0).sum())
        agreement = max(n_up, n_down) / len(vals)
        genes.append(GeneConcordance(
            gene=gene, n_datasets=len(byds), n_up=n_up, n_down=n_down,
            direction_agreement=round(agreement, 3),
            mean_log2fc=round(float(vals.mean()), 3),
            per_dataset={k: round(v, 3) for k, v in byds.items()},
        ))

    # rank: more supporting datasets, then stronger agreement, then larger effect
    genes.sort(key=lambda g: (g.n_datasets, g.direction_agreement, abs(g.mean_log2fc)), reverse=True)

    note = (f"{len(genes)} gene(s) differentially expressed in >= {min_datasets} datasets "
            f"(padj<{padj}, |log2fc|>={lfc}).")
    return Concordance(n_datasets=len(ok), genes=genes[:top], note=note)
