"""Report builder: assemble the OracleReport and render honest, cited Markdown.

The narrative is optionally phrased by the LLM but is constrained to the numbers already
computed; when no LLM is configured a deterministic template is used. Either way the numbers
are identical because they come from the data.
"""
from __future__ import annotations

from typing import Optional

from .llm import LLM
from .types import (Concordance, DatasetMatch, DEResult, ExperimentSpec, OracleReport,
                    Projection)


def build_report(spec: ExperimentSpec, matches: list[DatasetMatch], de_results: list[DEResult],
                 concordance: Optional[Concordance], projection: Projection,
                 llm: Optional[LLM] = None, warnings: Optional[list[str]] = None,
                 comparison: Optional[dict] = None,
                 enrichment: Optional[dict] = None,
                 heatmap: Optional[dict] = None) -> OracleReport:
    report = OracleReport(spec=spec, matches=matches, de_results=de_results,
                          concordance=concordance, projection=projection,
                          comparison=comparison, enrichment=enrichment, heatmap=heatmap,
                          warnings=warnings or [])
    report.narrative = _narrative(report, llm)
    return report


def _narrative(report: OracleReport, llm: Optional[LLM]) -> str:
    llm = llm if llm is not None else LLM()
    if llm.available:
        text = llm.complete(_SYS, _facts_block(report))
        if text:
            return text.strip()
    return _template_narrative(report)


_SYS = (
    "You are a cautious bioinformatics reviewer. Write a short (4-8 sentence) plain-language "
    "summary of an experiment sanity-check, using ONLY the facts provided. Do not introduce "
    "any numbers, gene names, or datasets not present in the facts. Lead with whether prior "
    "public data supports the plan, name the closest datasets by accession, state the main "
    "caveats/mismatches, and end with the feasibility read. Never overstate confidence."
)


def _facts_block(report: OracleReport) -> str:
    d = report.to_dict()
    lines = [f"SPEC: {d['spec']}"]
    lines.append("MATCHES:")
    for m in report.matches[:5]:
        mm = "; ".join(f"{x.field}:{x.severity}" for x in m.mismatches) or "none"
        lines.append(f"  {m.accession} sim={m.similarity} organism={m.organism} "
                     f"omics={m.omics_type} n={m.n_samples} mismatches=[{mm}] title={m.title[:80]}")
    lines.append(f"REANALYZED: {[dd.to_dict() for dd in report.de_results]}")
    if report.concordance:
        lines.append(f"CONCORDANCE: {report.concordance.note}; "
                     f"top={[g.gene for g in report.concordance.genes[:8]]}")
    lines.append(f"PROJECTION: {report.projection.to_dict()}")
    return "\n".join(lines)


def _template_narrative(report: OracleReport) -> str:
    usable = [m for m in report.matches if not m.blocking]
    parts: list[str] = []
    if usable:
        top = usable[0]
        parts.append(f"The closest public dataset is {top.accession} "
                     f"(similarity {top.similarity}), {top.organism or 'organism n/a'}, "
                     f"{top.n_samples if top.n_samples is not None else '?'} samples.")
    else:
        parts.append("No candidate dataset cleared the blocking-mismatch filter, so prior "
                     "public evidence for this exact plan is weak.")
    reanalyzed = [d for d in report.de_results if d.status == "ok"]
    if reanalyzed:
        rng = report.projection.expected_n_de
        parts.append(f"Reanalyzing {len(reanalyzed)} prior study(ies) gives an expected "
                     f"DE-gene count roughly in {list(rng) if rng else 'n/a'} "
                     f"(padj<0.05, |log2FC|>=1).")
    else:
        parts.append("Reprocessed counts were not available for reanalysis, so DE-count and "
                     "effect-size priors are not estimated (metadata-only pass).")
    if report.concordance and report.concordance.genes:
        genes = ", ".join(g.gene for g in report.concordance.genes[:5])
        parts.append(f"Genes most consistently differential across studies: {genes}.")
    parts.append(f"Overall feasibility read: {report.projection.feasibility}. "
                 f"{report.projection.power_read}")
    return " ".join(parts)


def render_markdown(report: OracleReport) -> str:
    s = report.spec
    P = report.projection
    L: list[str] = []
    L.append("# Omics Oracle — experiment sanity check")
    L.append("")
    L.append(f"> _{s.raw_prompt}_")
    L.append("")
    L.append(f"**Parsed spec** (`{s.parser}` parser): organism=**{s.organism or '—'}**, "
             f"tissue=**{s.tissue or '—'}**, omics=**{s.omics_type}**, "
             f"contrast=**{(s.contrast_case or '?')} vs {(s.contrast_control or '?')}**, "
             f"platform=**{s.platform or '—'}**, N/group=**{s.n_per_group or '—'}**.")
    L.append("")

    if report.warnings:
        L.append("**Warnings:** " + "; ".join(report.warnings))
        L.append("")

    L.append("## Narrative")
    L.append(report.narrative or "_(none)_")
    L.append("")

    L.append("## Matched datasets (ranked)")
    if not report.matches:
        L.append("_No candidate datasets returned by the finder._")
    else:
        judged = any(m.judge for m in report.matches)
        tiered = any(m.match_tier for m in report.matches)
        head = "| Rank | Accession | Sim |"
        sep = "|---|---|---|"
        if tiered:
            head += " Tier |"; sep += "---|"
        head += " Organism | Omics | N | Blocking? |"
        sep += "---|---|---|---|"
        if judged:
            head += " Judge |"; sep += "---|"
        head += " Mismatches |"; sep += "---|"
        L.append(head); L.append(sep)
        for i, m in enumerate(report.matches, 1):
            mm = "; ".join(f"{x.field}({x.severity})" for x in m.mismatches) or "—"
            row = f"| {i} | [{m.accession}]({m.url}) | {m.similarity} |"
            if tiered:
                row += f" {m.match_tier or '—'} |"
            row += (f" {m.organism or '—'} | {m.omics_type} | "
                    f"{m.n_samples if m.n_samples is not None else '—'} | "
                    f"{'⛔' if m.blocking else '—'} |")
            if judged:
                jv = (m.judge or {}).get("verdict", "—")
                row += f" {jv} |"
            row += f" {mm} |"
            L.append(row)
        L.append("")
        # judge reasons for anything not cleanly relevant — the anti-fabrication trail
        jnotes = [(m.accession, m.judge) for m in report.matches
                  if m.judge and m.judge.get("verdict") in ("partial", "irrelevant", "uncertain")]
        if jnotes:
            L.append("**Relevance verification (LLM judge):**")
            for acc, j in jnotes:
                L.append(f"- {acc}: **{j.get('verdict')}** (conf {j.get('confidence')}) — {j.get('reason')}")
            L.append("")
        # spell out the most important mismatch notes (the 'HEK293 not colorectal' caveats)
        notes = [x.note for m in report.matches for x in m.mismatches
                 if x.severity in ("blocking", "major") and x.note]
        if notes:
            L.append("**Mismatch caveats:**")
            for n in dict.fromkeys(notes):
                L.append(f"- {n}")
            L.append("")

    L.append("## Per-dataset reanalysis")
    if not report.de_results:
        L.append("_No datasets were reanalyzed._")
    else:
        for d in report.de_results:
            src = f"{d.tier}: {d.source_db}"
            if getattr(d, "method", ""):
                src += f"; {d.method}"
            label = f" _[{src}]_" if getattr(d, "source_db", None) else ""
            L.append(f"- **{d.accession}**{label} — {d.status.upper()}: {d.note}")
            if d.qc:
                L.append(f"  - QC: {d.qc.n_samples} samples, median lib size "
                         f"{d.qc.median_library_size:,.0f}, median zero-fraction {d.qc.median_zero_fraction}."
                         + (" Flags: " + "; ".join(d.qc.flags) if d.qc.flags else ""))
            sel = getattr(d, "sample_selection", None)
            if sel and sel.get("selection"):
                s = sel["selection"]
                line = (f"  - Sample selection ({sel.get('outcome')}): "
                        f"case={s.get('n_case')}, control={s.get('n_control')}, "
                        f"excluded={s.get('n_excluded')}")
                v = sel.get("verdict")
                if v:
                    line += f" — judge {'✅ agreed' if v.get('agree') else '⚠️ disagreed'}"
                    if not v.get("agree") and v.get("concerns"):
                        line += f" ({'; '.join(v['concerns'][:2])})"
                L.append(line)
            if d.status == "ok":
                L.append(f"  - DE: {d.n_de()} genes at padj<0.05 & |log2FC|>=1 "
                         f"(case={d.n_case}, control={d.n_control}).")
    L.append("")

    if report.concordance:
        L.append("## Cross-experiment concordance")
        L.append(report.concordance.note)
        if report.concordance.genes:
            L.append("")
            L.append("| Gene | Datasets | Up | Down | Agreement | Mean log2FC |")
            L.append("|---|---|---|---|---|---|")
            for g in report.concordance.genes[:20]:
                L.append(f"| {g.gene} | {g.n_datasets} | {g.n_up} | {g.n_down} | "
                         f"{g.direction_agreement} | {g.mean_log2fc} |")
        L.append("")

    if report.comparison:
        c = report.comparison
        L.append("## Cross-dataset comparison (Venn + gene search)")
        L.append("Reanalyzed studies: " + ", ".join(
            f"{d['accession']} [{d['tier']}] ({d['n_de']} DE)" for d in c["datasets"]))
        if c.get("pairwise"):
            L.append("")
            L.append("**Pairwise DE-gene overlap:**")
            for p in sorted(c["pairwise"], key=lambda x: x["shared"], reverse=True)[:10]:
                L.append(f"- {p['a']} ∩ {p['b']}: **{p['shared']}** shared "
                         f"(Jaccard {p['jaccard']}; {p['a_only']}/{p['b_only']} unique)")
        top = [g for g in c["genes"] if g["n_datasets"] >= 2][:15]
        if top:
            L.append("")
            L.append("**Genes DE across the most studies** (search by gene in the JSON `comparison.genes`):")
            L.append("| Gene | #datasets | up/down | mean log2FC |")
            L.append("|---|---|---|---|")
            for g in top:
                L.append(f"| {g['gene']} | {g['n_datasets']} | {g['up']}/{g['down']} | {g['mean_log2fc']} |")
        if c.get("user"):
            u = c["user"]
            L.append("")
            L.append(f"**Your gene list ({u['n_genes']} genes) overlap:** " + "; ".join(
                f"{v['accession']} {v['shared']} shared (J={v['jaccard']})" for v in u["vs_datasets"]))
        L.append("")

    if report.heatmap and report.heatmap.get("genes"):
        h = report.heatmap
        cols = " × ".join(f"{d['accession']}[{d['match_tier'] or '?'}]" for d in h["datasets"])
        L.append("## Concordance heatmap (gene × dataset log2FC)")
        L.append(f"Matrix: {len(h['genes'])} genes ({h['scope']}) × {len(h['datasets'])} datasets "
                 f"({cols}). Full data in JSON `heatmap` (UI toggles Top-N / Up-Down over `heatmap.genes`).")
        L.append("")

    if report.enrichment and report.enrichment.get("sources"):
        e = report.enrichment
        L.append(f"## Pathway enrichment ({e['library']})")
        if e.get("shared_pathways"):
            L.append("**Shared pathways** (enriched across studies):")
            for p in e["shared_pathways"][:10]:
                L.append(f"- **{p['term']}** — significant in {p['n_significant']} "
                         f"({', '.join(p['sources_significant'])})")
        if e.get("surfaced"):
            L.append("")
            L.append("**Worth a look — strong in some studies, sub-threshold (but genes moving) in others:**")
            for p in e["surfaced"][:8]:
                subs = ", ".join(f"{s['source']}({s['overlap_n']} genes)" for s in p["subthreshold_in"])
                L.append(f"- **{p['term']}** — significant in {', '.join(p['significant_in'])}; "
                         f"sub-threshold in {subs}")
        if e.get("note"):
            L.append("")
            L.append(f"_{e['note']}_")
        L.append("")

    L.append("## Projection")
    if P.expected_n_de:
        L.append(f"- **Expected #DE genes:** ~{P.expected_n_de[0]}–{P.expected_n_de[1]} "
                 f"(from reanalyzed priors).")
    if P.effect_size_prior:
        L.append(f"- **Effect-size prior:** {P.effect_size_prior}.")
    L.append(f"- **Power / feasibility:** {P.power_read}")
    L.append(f"- **Overall feasibility:** **{P.feasibility}**.")
    if P.qc_pitfalls:
        L.append("- **Expected QC pitfalls:**")
        for p in P.qc_pitfalls:
            L.append(f"  - {p}")
    if P.caveats:
        L.append("- **Caveats:**")
        for c in P.caveats:
            L.append(f"  - {c}")
    L.append("")
    L.append(f"<sub>schema {report.schema_version} · real data only · "
             f"projections are cited to the datasets above and are not a substitute for "
             f"formal power analysis.</sub>")
    return "\n".join(L)
