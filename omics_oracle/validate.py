"""Deterministic narrative validation — anti-hallucination guard.

The LLM phrases the narrative; the numbers and genes come from data. This checks that every
gene-symbol-looking token the narrative mentions actually appears in the analyzed data (the
harmonized symbols from comparison/concordance/heatmap). Anything else is surfaced as a
"verify" warning so an invented gene can't slip through unnoticed. Conservative by design: it
flags (never silently rewrites), and a domain stoplist keeps assay/pathway acronyms out.
"""
from __future__ import annotations

import re
from typing import Iterable

# uppercase tokens that look gene-ish but are assay/pathway/domain terms, not genes to validate
_STOP = {
    "RNA", "DNA", "MRNA", "CDNA", "DE", "DEG", "DEGS", "GEO", "SRA", "NCBI", "TPM", "FPKM",
    "RPKM", "CPM", "VST", "QC", "LLM", "PAI", "KEGG", "GO", "ID", "IDS", "FC", "LOG2FC",
    "PADJ", "FDR", "PCA", "ARCHS4", "HGNC", "GSM", "GSE", "GDS", "GPL", "US", "USA", "EU",
    "MSS", "MSI", "PFS", "OS", "KO", "KD", "WT", "CRC", "PDAC", "FLOX", "EV", "EVS", "AUC",
    "MAPK", "PI3K", "AKT", "MTOR", "NF", "TCA", "ATP", "II", "III", "IV", "3D", "OOC",
    "N", "AND", "OR", "THE", "A", "AN", "VS", "HCT", "BXPC3", "SW480", "KBM7",
}
_GENE_TOKEN = re.compile(r"\b[A-Z][A-Z0-9]{1,8}\b")
_ACCESSION = re.compile(r"^(GSE|GDS|GSM|GPL|SRP|SRR|PRJNA)\d+$")


def collect_allowed_genes(*sources: Iterable) -> set[str]:
    """Union of legitimate gene symbols from comparison/concordance/heatmap structures."""
    allowed: set[str] = set()
    for src in sources:
        if not src:
            continue
        if isinstance(src, dict):
            for g in src.get("genes", []) or []:
                if isinstance(g, dict) and g.get("gene"):
                    allowed.add(str(g["gene"]).upper())
        else:  # iterable of objects with .gene
            for g in src:
                gene = getattr(g, "gene", None)
                if gene:
                    allowed.add(str(gene).upper())
    return allowed


def validate_narrative(text: str, allowed_genes: set[str]) -> list[str]:
    """Return gene-symbol-looking tokens in the narrative not present in the analyzed data."""
    if not text:
        return []
    flagged: list[str] = []
    seen: set[str] = set()
    for tok in _GENE_TOKEN.findall(text):
        up = tok.upper()
        if up in seen or up in _STOP or _ACCESSION.match(up):
            continue
        # require a gene-symbol shape: contains a digit, or is a plausible 3-8 letter symbol
        looks_gene = any(c.isdigit() for c in up) or (3 <= len(up) <= 8)
        if not looks_gene:
            continue
        if up not in allowed_genes:
            seen.add(up)
            flagged.append(tok)
    return flagged
