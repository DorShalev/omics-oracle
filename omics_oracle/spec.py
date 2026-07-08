"""Experiment-spec parser: prompt -> ExperimentSpec.

Primary path is the LLM (structured JSON, validated into the dataclass). Fallback is a
deterministic keyword heuristic so the engine runs with no key. The chosen path is recorded
in ``spec.parser`` so honesty is preserved downstream.
"""
from __future__ import annotations

import re
from typing import Optional

from .llm import LLM
from .types import ExperimentSpec, OMICS_TYPES

# ---- canonical vocabularies for the heuristic path -------------------------------------

_ORGANISMS = {
    "human": "Homo sapiens", "homo sapiens": "Homo sapiens", "patient": "Homo sapiens",
    "mouse": "Mus musculus", "murine": "Mus musculus", "mus musculus": "Mus musculus",
    "rat": "Rattus norvegicus", "yeast": "Saccharomyces cerevisiae",
    "zebrafish": "Danio rerio", "drosophila": "Drosophila melanogaster", "fly": "Drosophila melanogaster",
}

_OMICS = {
    "rna-seq": "rna-seq", "rnaseq": "rna-seq", "rna seq": "rna-seq", "transcriptom": "rna-seq",
    "bulk rna": "rna-seq", "mrna-seq": "rna-seq",
    "scrna": "scrna-seq", "single-cell rna": "scrna-seq", "single cell rna": "scrna-seq",
    "atac": "atac-seq",
    "proteom": "proteomics", "mass spec": "proteomics", "lc-ms": "proteomics", "dia-nn": "proteomics",
    "metabolom": "metabolomics", "lipidom": "metabolomics",
}

_PLATFORMS = ["illumina", "novaseq", "nextseq", "hiseq", "nanopore", "pacbio", "10x", "bgi", "mgi"]

# tissue/condition tokens worth carrying into search (kept broad; the finder does overlap)
_STOP = set("""a an the of for and or with without vs versus from to in on at by is are was
    experiment experiments study studies sample samples data dataset patients patient control
    controls case cases healthy planned rna seq sequencing profiling expression per group each
    about approximately roughly around n= reads""".split())


def parse_spec(prompt: str, llm: Optional[LLM] = None) -> ExperimentSpec:
    llm = llm if llm is not None else LLM()
    if llm.available:
        spec = _parse_with_llm(prompt, llm)
        if spec is not None:
            return spec
    return _parse_heuristic(prompt)


# ---- LLM path --------------------------------------------------------------------------

_SYS = (
    "You extract a structured description of a PLANNED omics experiment from the user's text. "
    "Return ONLY a JSON object with keys: organism (canonical Latin binomial or null), tissue "
    "(free text or null), contrast_case (the perturbed/disease group or null), contrast_control "
    "(the baseline group or null), omics_type (one of "
    + ", ".join(OMICS_TYPES) + "), platform (or null), n_per_group (integer or null), keywords "
    "(array of 3-8 salient search terms). Do not invent details that are not implied by the text."
)


def _parse_with_llm(prompt: str, llm: LLM) -> Optional[ExperimentSpec]:
    obj = llm.complete_json(_SYS, prompt)
    if not obj:
        return None
    omics = str(obj.get("omics_type") or "unknown").lower()
    if omics not in OMICS_TYPES:
        omics = _norm_omics(omics) or "unknown"
    n = obj.get("n_per_group")
    try:
        n = int(n) if n is not None else None
    except (TypeError, ValueError):
        n = None
    kws = obj.get("keywords") or []
    if not isinstance(kws, list):
        kws = []
    return ExperimentSpec(
        raw_prompt=prompt,
        organism=obj.get("organism") or None,
        tissue=obj.get("tissue") or None,
        contrast_case=obj.get("contrast_case") or None,
        contrast_control=obj.get("contrast_control") or None,
        omics_type=omics,
        platform=obj.get("platform") or None,
        n_per_group=n,
        parser="llm",
        keywords=[str(k) for k in kws][:8],
    )


# ---- heuristic path --------------------------------------------------------------------

def _norm_omics(text: str) -> Optional[str]:
    t = text.lower()
    for needle, canon in _OMICS.items():
        if needle in t:
            return canon
    return None


def _parse_heuristic(prompt: str) -> ExperimentSpec:
    low = prompt.lower()

    organism = None
    for needle, canon in _ORGANISMS.items():
        if re.search(rf"\b{re.escape(needle)}\b", low):
            organism = canon
            break

    omics = _norm_omics(low) or "unknown"

    platform = next((p for p in _PLATFORMS if p in low), None)
    if platform:
        platform = platform.capitalize()

    n_per_group = _extract_n(low)

    case, control = _extract_contrast(prompt)
    tissue = _extract_tissue(low)
    keywords = _keywords(low)

    return ExperimentSpec(
        raw_prompt=prompt, organism=organism, tissue=tissue,
        contrast_case=case, contrast_control=control, omics_type=omics,
        platform=platform, n_per_group=n_per_group, parser="heuristic",
        keywords=keywords,
    )


def _extract_n(low: str) -> Optional[int]:
    # "20 per arm/group/condition", "n=20", "20 patients ... 20 controls" -> take the smaller-ish signal
    m = re.search(r"(\d{1,4})\s*(?:per|/)\s*(?:arm|group|condition|cohort)", low)
    if m:
        return int(m.group(1))
    m = re.search(r"\bn\s*=\s*(\d{1,4})", low)
    if m:
        return int(m.group(1))
    m = re.search(r"(\d{1,4})\s+(?:patients?|subjects?|samples?|donors?|replicates?)", low)
    if m:
        return int(m.group(1))
    return None


def _extract_contrast(prompt: str) -> tuple[Optional[str], Optional[str]]:
    # "<case> vs <control>" / "<case> versus <control>" / "<case> compared to <control>"
    m = re.search(r"([A-Za-z0-9'\-/ ]{2,40}?)\s+(?:vs\.?|versus|compared to|against)\s+([A-Za-z0-9'\-/ ]{2,40})",
                  prompt, re.IGNORECASE)
    if not m:
        return None, None
    case = _trim_phrase(m.group(1))
    control = _trim_phrase(m.group(2))
    return case or None, control or None


# Words to strip from a contrast phrase — assay/organism/filler, but NOT group labels like
# "control"/"tumor"/"treated" (those are exactly what we want to keep).
_CONTRAST_STRIP = set("""rna-seq rnaseq rna seq sequencing scrna atac proteomics metabolomics
    lc-ms mass spec expression profiling transcriptomics bulk mrna of the a an in from for
    with human mouse murine rat yeast zebrafish drosophila fly patients patient samples sample
    cells cell tissue tissues""".split())


def _trim_phrase(s: str) -> str:
    toks = [t for t in re.split(r"\s+", s.strip()) if t]
    toks = [t for t in toks if t.lower() not in _CONTRAST_STRIP]
    # keep the trailing 1-3 tokens that describe the group
    return " ".join(toks[-3:]).strip(" ,.;:")


def _extract_tissue(low: str) -> Optional[str]:
    tissues = ["plasma", "serum", "exosome", "ev", "extracellular vesicle", "tumor", "tumour",
               "biopsy", "blood", "liver", "brain", "colorectal", "colon", "breast", "lung",
               "pancrea", "kidney", "muscle", "skin", "bone marrow", "csf", "urine"]
    hits = [t for t in tissues if t in low]
    return " ".join(dict.fromkeys(hits)) or None if hits else None


def _keywords(low: str) -> list[str]:
    words = re.findall(r"[a-z0-9][a-z0-9'\-]{2,}", low)
    seen: list[str] = []
    for w in words:
        if w in _STOP or w.isdigit():
            continue
        if w not in seen:
            seen.append(w)
    return seen[:8]
