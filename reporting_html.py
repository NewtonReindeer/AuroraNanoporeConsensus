"""HTML report rendering for HCV and HBV analyses."""
import os
import re
from datetime import datetime
from html import escape as html_escape

from .config import (
    HBV_MEDICATIONS,
    HBV_MIN_COVERAGE,
    HBV_MIN_REPORTABLE_COVERAGE,
    HBV_MIN_VARIANT_FREQUENCY,
    HBV_S_ESCAPE_ASSOCIATIONS,
    HCV_HTML_SCRIPT_VERSION,
    HCV_MIN_COVERAGE,
    HCV_MIN_VARIANT_FREQUENCY,
    HCV_NS5A_REPRESENTATIVE_MEDICATIONS,
)
from .filesystem import (
    ensure_dir,
    find_hbv_report_logo,
    image_file_to_data_uri,
    run_samples_folder,
    safe_report_filename,
)
from .hbv import (
    hbv_combination_interpretations,
    hbv_detected_calls,
    hbv_status_for_medication,
    hbv_swedish_evidence_label,
    hbv_swedish_interpretation,
)


def hcv_detected_calls(row):
    """Flatten reportable HCV NS5A calls already produced by the existing caller."""
    calls = []
    for detail in row.get("hcv_debug_details", []) or []:
        label = str(detail.get("position", "") or "")
        reference_aa = str(detail.get("reference_aa", "") or "")
        depth = int(detail.get("depth", 0) or 0) if str(detail.get("depth", "")).strip() else 0
        reported = str(detail.get("reported_calls", "") or "").strip()
        if not label or not reported or reported.lower().startswith(("auto-frame failed", "hcv resistance analysis failed")):
            continue
        for part in reported.split("|"):
            match = re.search(r"([A-Z*])\s*:\s*([0-9.]+)%", part.strip())
            if not match:
                continue
            alt_aa = match.group(1)
            frequency = float(match.group(2))
            if alt_aa == reference_aa:
                continue
            calls.append({
                "position": label,
                "mutation": f"{label}{alt_aa}",
                "reference_aa": reference_aa,
                "alt_aa": alt_aa,
                "frequency": frequency,
                "depth": depth,
                "count": int(round(depth * frequency / 100.0)) if depth else 0,
                "population": "Majoritetspopulation" if frequency >= 50.0 else "Minoritetspopulation",
            })
    return calls

def hcv_html_summary(row, calls):
    quality = str(row.get("hcv_quality", "") or "")
    interpretation = str(row.get("hcv_interpretation", "") or "")
    if quality == "Not Included":
        return "Subtypen ingår inte i den validerade NS5A-panelen i denna analysversion."
    if quality != "✓":
        return "Resultatet är inte fullständigt bedömbart på grund av otillräcklig paneltäckning eller analysfel."
    if not calls:
        return "Inga NS5A-resistensassocierade substitutioner påvisades över rapporteringströskeln."
    findings = ", ".join(f"{item['mutation']} ({item['frequency']:.1f} %)" for item in calls)
    return f"Påvisade HCV NS5A-RAS: {findings}. Samlad bedömning: {interpretation}."

def hcv_status_class(row, calls):
    quality = str(row.get("hcv_quality", "") or "")
    if quality not in ("✓",):
        return "status-unknown"
    if any(item["frequency"] >= 50.0 for item in calls):
        return "status-high"
    if calls:
        return "status-medium"
    return "status-good"

def hcv_panel_coverage_html(row):
    rows = []
    for detail in row.get("hcv_debug_details", []) or []:
        label = str(detail.get("position", "") or "")
        if not label or label == "panel":
            continue
        depth = int(detail.get("depth", 0) or 0) if str(detail.get("depth", "")).strip() else 0
        reported = str(detail.get("reported_calls", "") or "").strip() or "Inget panelrelevant fynd"
        qc = "Godkänd" if depth >= HCV_MIN_COVERAGE else "Otillräcklig"
        rows.append(
            "<tr>"
            f"<td>{html_escape(label)}</td>"
            f"<td>{html_escape(str(detail.get('reference_aa', '') or ''))}</td>"
            f"<td>{html_escape(str(detail.get('reference_codon', '') or ''))}</td>"
            f"<td>{depth}</td>"
            f"<td>{html_escape(reported)}</td>"
            f"<td>{qc}</td>"
            "</tr>"
        )
    return "".join(rows) or "<tr><td colspan='6' class='empty'>Paneldata saknas.</td></tr>"

def hcv_medication_rows_html(row, calls):
    """GLUE-inspired medication view without inventing drug-specific scores."""
    if str(row.get("hcv_quality", "") or "") != "✓":
        assessment = "Ej bedömbar på grund av otillräcklig paneltäckning"
        findings = "—"
        css = "drug-unknown"
    elif calls:
        assessment = "NS5A-RAS påvisad – läkemedelsspecifik effekt ej beräknad"
        findings = "; ".join(f"{x['mutation']} ({x['frequency']:.1f} %)" for x in calls)
        css = "drug-ras"
    else:
        assessment = "Ingen panel-RAS påvisad"
        findings = "—"
        css = "drug-none"
    return "".join(
        "<tr>"
        f"<td><strong>{html_escape(drug)}</strong></td>"
        f"<td><span class='drug-status {css}'>{html_escape(assessment)}</span></td>"
        f"<td>{html_escape(findings)}</td>"
        "</tr>"
        for drug in HCV_NS5A_REPRESENTATIVE_MEDICATIONS
    )

def write_hcv_html_report(row, report_path, run_tag):
    calls = hcv_detected_calls(row)
    panel_subtype = str(row.get('hcv_panel_subtype', '') or '').strip().upper()
    blast_subtype = str(row.get('blast_subtype', '') or row.get('blast_genotype_display', '') or '').strip()
    if panel_subtype:
        ns5a_analysis_status = f'{panel_subtype}-panel'
    elif blast_subtype:
        ns5a_analysis_status = f'Ej inkluderad för subtyp {blast_subtype}'
    else:
        ns5a_analysis_status = 'Ej analyserad'
    logo_uri = image_file_to_data_uri(find_hbv_report_logo(row, report_path))
    logo_html = f"<img class='report-logo' src='{logo_uri}' alt='Logotyp'>" if logo_uri else ""
    summary = hcv_html_summary(row, calls)
    mutation_rows = []
    for item in calls:
        mutation_rows.append(
            "<tr>"
            f"<td><strong>{html_escape(item['mutation'])}</strong></td>"
            f"<td>{item['frequency']:.1f} %</td>"
            f"<td>{html_escape(item['population'])}</td>"
            f"<td>{item['count']} / {item['depth']}</td>"
            f"<td>{html_escape(str(row.get('hcv_panel_subtype', '') or 'Ej analyserad'))}</td>"
            "</tr>"
        )
    if not mutation_rows:
        mutation_rows.append("<tr><td colspan='5' class='empty'>Inga rapporterbara NS5A-RAS påvisades.</td></tr>")

    html = f"""<!DOCTYPE html>
<html lang='sv'>
<head>
<meta charset='UTF-8'>
<meta name='viewport' content='width=device-width, initial-scale=1'>
<title>HCV-resistensrapport – {html_escape(str(row.get('sample','')))}</title>
<style>
:root {{ --blue:#dbeeff; --blue-dark:#185b8d; --line:#c9d3dc; --text:#1e2a32; --muted:#5d6b75; }}
* {{ box-sizing:border-box; }}
body {{ font-family:Arial,Helvetica,sans-serif; color:var(--text); margin:0; background:#eef2f5; }}
.report {{ max-width:1050px; margin:24px auto; background:#fff; padding:34px 42px; box-shadow:0 3px 18px rgba(0,0,0,.12); }}
h1 {{ margin:0 0 5px; font-size:28px; color:#123e5f; }}
.header {{ display:flex; justify-content:space-between; align-items:flex-start; gap:24px; margin-bottom:20px; }}
.title-block {{ flex:1; }} .report-logo {{ max-width:220px; max-height:105px; object-fit:contain; }}
.subtitle {{ color:var(--muted); }} .summary {{ background:#f4f8fb; border-left:6px solid var(--blue-dark); padding:12px 15px; margin:14px 0 18px; }}
h2 {{ font-size:18px; margin:28px 0 10px; border-bottom:2px solid var(--blue-dark); padding-bottom:6px; }}
table {{ width:100%; border-collapse:collapse; font-size:14px; }} th {{ background:var(--blue); text-align:left; }}
th,td {{ border:1px solid var(--line); padding:9px 10px; vertical-align:top; }} .meta td:first-child {{ width:28%; background:#f4f8fb; font-weight:bold; }}
.status {{ padding:14px 16px; border-radius:6px; font-weight:bold; margin:16px 0; }}
.status-good {{ background:#dff3e3; border-left:7px solid #3e9d54; }} .status-medium {{ background:#fff4c7; border-left:7px solid #d0a600; }}
.status-high {{ background:#ffe1d6; border-left:7px solid #d86032; }} .status-unknown {{ background:#e3e6e8; border-left:7px solid #6f777d; }}
.drug-status {{ font-weight:bold; }} .drug-none {{ color:#2c7b3e; }} .drug-ras {{ color:#a17900; }} .drug-unknown {{ color:#6f777d; }}
.empty {{ color:var(--muted); font-style:italic; }} .notice {{ background:#f5f5f5; border-left:5px solid #777; padding:12px 14px; margin-top:26px; font-size:12px; }}
.small {{ font-size:12px; color:var(--muted); word-break:break-all; }}
.report-footer {{
    margin-top: 24px;
    padding-top: 10px;
    border-top: 1px solid var(--line);
    color: var(--muted);
    font-size: 10px;
    line-height: 1.35;
    text-align: center;
}}
@media print {{ body {{ background:#fff; }} .report {{ box-shadow:none; margin:0; max-width:none; }} .report-logo {{ max-height:85px; }} }}
</style>
</head>
<body><main class='report'>
<div class='header'><div class='title-block'>
<h1>Hepatit C-virus – NS5A-resistensrapport</h1>
<div class='subtitle'>Läsnivåbaserad RAS-analys · Endast för forskningsbruk</div>
</div>{logo_html}</div>
<table class='meta'>
<tr><td>Prov-ID</td><td>{html_escape(str(row.get('sample','')))}</td></tr>
<tr><td>Genotyp/subtyp (BLAST)</td><td><strong>{html_escape(str(row.get('blast_genotype_display','') or 'Ej fastställd'))}</strong></td></tr>
<tr><td>BLAST toppträff</td><td><strong>{html_escape(str(row.get('blast_genotype_match_accession','') or row.get('blast_top_accession','') or 'Ej tillgänglig'))}</strong> · {html_escape(str(row.get('blast_genotype_match_organism','') or row.get('blast_top_organism','')))} · identitet {html_escape(str(row.get('blast_genotype_match_identity','') or row.get('blast_top_identity','')))} % · query coverage {html_escape(str(row.get('blast_genotype_match_query_coverage','') or row.get('blast_top_query_coverage','')))} %</td></tr>
<tr><td>NS5A-resistensanalys</td><td>{html_escape(ns5a_analysis_status)}</td></tr>
<tr><td>Bästa mappningsreferens</td><td>{html_escape(str(row.get('best_reference','')))}</td></tr>
<tr><td>Rapportdatum</td><td>{datetime.now().strftime('%Y-%m-%d %H:%M')}</td></tr>
<tr><td>Körning</td><td>{html_escape(run_tag)}</td></tr>
<tr><td>Rapporteringströskel</td><td>≥ {HCV_MIN_VARIANT_FREQUENCY * 100:.1f} %</td></tr>
<tr><td>Minsta kodontäckning</td><td>{HCV_MIN_COVERAGE} läsningar</td></tr>
<tr><td>Median nollskild täckning</td><td>{html_escape(str(row.get('hcv_median_nonzero_coverage','')))}×</td></tr>
<tr><td>Antal mappade läsningar</td><td>{html_escape(str(row.get('mapped_reads','')))}</td></tr>
</table>
<div class='summary'><strong>Sammanfattning:</strong> {html_escape(summary)}</div>
<div class='status {hcv_status_class(row, calls)}'>{html_escape(str(row.get('hcv_interpretation','')))}</div>

<h2>NS5A-läkemedelsöversikt</h2>
<table><thead><tr><th>Läkemedel</th><th>Panelbedömning</th><th>Påvisade RAS</th></tr></thead>
<tbody>{hcv_medication_rows_html(row, calls)}</tbody></table>
<p class='small'>Tabellen är en klassbaserad screeningvy. Den aktuella panelen beräknar inte GLUE-liknande läkemedelsspecifika resistenspoäng eller behandlingsrekommendationer.</p>

<h2>Påvisade NS5A-resistensassocierade substitutioner</h2>
<table><thead><tr><th>Mutation</th><th>Andel</th><th>Population</th><th>Mutantläsningar / djup</th><th>NS5A-panel</th></tr></thead>
<tbody>{''.join(mutation_rows)}</tbody></table>

<h2>Paneltäckning och positionskontroll</h2>
<table><thead><tr><th>Position</th><th>Referens-AA</th><th>Referenskodon</th><th>Kodontäckning</th><th>Rapporterat resultat</th><th>QC</th></tr></thead>
<tbody>{hcv_panel_coverage_html(row)}</tbody></table>

<h2>Teknisk kvalitet</h2>
<table class='meta'>
<tr><td>Kvalitetsstatus</td><td>{html_escape(str(row.get('hcv_quality','')))}</td></tr>
<tr><td>Övergripande QC</td><td>{html_escape(str(row.get('qc_verdict','')))}</td></tr>
<tr><td>QC-kommentar</td><td>{html_escape(str(row.get('qc_reason','')))}</td></tr>
<tr><td>Positioner med otillräcklig täckning</td><td>{html_escape(str(row.get('hcv_low_coverage','')))}</td></tr>
<tr><td>NS5A-koordinat</td><td>start {html_escape(str(row.get('hcv_ns5a_start_used','')))}</td></tr>
<tr><td>Genom-/måltäckning</td><td>{html_escape(str(row.get('target_or_genome_coverage_percent','')))} %</td></tr>
<tr><td>Medeldjup</td><td>{html_escape(str(row.get('mean_depth','')))}×</td></tr>
<tr><td>Mediandjup</td><td>{html_escape(str(row.get('median_depth','')))}×</td></tr>
<tr><td>Positioner ≥10×</td><td>{html_escape(str(row.get('positions_ge_10x_percent','')))} %</td></tr>
<tr><td>Positioner ≥30×</td><td>{html_escape(str(row.get('positions_ge_30x_percent','')))} %</td></tr>
<tr><td>Konsensuslängd</td><td>{html_escape(str(row.get('consensus_length','')))} nt</td></tr>
<tr><td>N-baser</td><td>{html_escape(str(row.get('n_count','')))} ({html_escape(str(row.get('n_percent','')))} %)</td></tr>
<tr><td>Medaka-modell</td><td>{html_escape(str(row.get('medaka_model','')))}</td></tr>
<tr><td>FASTQ-fil</td><td class='small'>{html_escape(str(row.get('input_fastq','')))}</td></tr>
<tr><td>Konsensusfil</td><td class='small'>{html_escape(str(row.get('consensus_file','')))}</td></tr>
<tr><td>BAM-fil</td><td class='small'>{html_escape(str(row.get('bam_file','')))}</td></tr>
<tr><td>Referensfil</td><td class='small'>{html_escape(str(row.get('reference_file','')))}</td></tr>
<tr><td>BLAST-databas</td><td class='small'>{html_escape(str(row.get('blast_database','')))}</td></tr>
<tr><td>BLAST-databasdatum</td><td>{html_escape(str(row.get('blast_database_date','')))}</td></tr>
<tr><td>BLAST-status</td><td>{html_escape(str(row.get('blast_status','')))} {html_escape(str(row.get('blast_message','')))}</td></tr>
<tr><td>Scriptversion</td><td>{html_escape(str(row.get('hcv_script_version','')))} / {HCV_HTML_SCRIPT_VERSION}</td></tr>
</table>
<div class='notice'><strong>Viktig begränsning:</strong> Rapporten visar endast de subtypspecifika NS5A-positioner och substitutioner som ingår i Aurora-panelen. Den ersätter inte klinisk behandlingsbedömning eller expertgranskning.</div>
<footer class='report-footer'>
Vid frågor, tekniska problem eller oväntade analysresultat, kontakta skriptets författare Joakim Bedner Stenbäck (joakim.stenback@vgregion.se).
</footer>

</main></body></html>"""
    ensure_dir(os.path.dirname(report_path), must_be_writable=True)
    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write(html)

def generate_hcv_html_reports(hcv_rows, output_folder, run_tag):
    html_dir = os.path.join(run_samples_folder(output_folder, run_tag), "HCV_HTML_reports")
    ensure_dir(html_dir, must_be_writable=True)
    index_rows = []
    for row in hcv_rows:
        sample = str(row.get("sample", "") or "prov")
        filename = f"{safe_report_filename(sample)}_{run_tag}_HCV.html"
        report_path = os.path.join(html_dir, filename)
        write_hcv_html_report(row, report_path, run_tag)
        row["hcv_html_report"] = report_path
        index_rows.append(
            "<tr>"
            f"<td><a href='{html_escape(filename)}'>{html_escape(sample)}</a></td>"
            f"<td>{html_escape(str(row.get('blast_genotype_display','') or row.get('hcv_subtype','')))}</td>"
            f"<td>{html_escape(str(row.get('hcv_interpretation','')))}</td>"
            f"<td>{html_escape(str(row.get('hcv_quality','')))}</td>"
            "</tr>"
        )
    index_path = os.path.join(html_dir, f"index_{run_tag}.html")
    index_html = f"""<!DOCTYPE html><html lang='sv'><head><meta charset='UTF-8'><title>HCV-rapporter {html_escape(run_tag)}</title>
<style>body{{font-family:Arial;margin:35px;color:#1e2a32}}table{{border-collapse:collapse;width:100%}}th{{background:#dbeeff}}th,td{{border:1px solid #c9d3dc;padding:9px;text-align:left}}a{{color:#185b8d}}</style></head>
<body><h1>HCV NS5A-resistensrapporter</h1><p>Körning: {html_escape(run_tag)}</p>
<table><thead><tr><th>Prov</th><th>Genotyp/subtyp (BLAST)</th><th>Bedömning</th><th>Kvalitet</th></tr></thead><tbody>{''.join(index_rows)}</tbody></table></body></html>"""
    with open(index_path, "w", encoding="utf-8") as handle:
        handle.write(index_html)
    print(f"HCV HTML reports saved to {html_dir}")
    return html_dir

def _hbv_compact_pattern_interpretation(pattern):
    """Return a concise professional-facing interpretation of an HBV rule result."""
    text = str(pattern.get("interpretation", "") or "").strip()
    level = str(pattern.get("level", "") or "")
    category = str(pattern.get("category", "") or "")

    # hbv_apply_population_context() prepends explanatory boilerplate for major/minor
    # populations. Remove that boilerplate in the HTML report while retaining the
    # mutation/rule-specific interpretation from the original rule.
    prefixes = [
        r"^Ett etablerat resistensassocierat mönster påvisades i majoritetspopulation\s*\(högsta variantandel [^)]+\)\.\s*Fyndet ger starkt genotypiskt stöd för nedsatt läkemedelskänslighet\.\s*",
        r"^[^.]+ påvisades, men den högsta rapporterade variantandelen är [^.]+\.\s*Ingen dominerande resistensassocierad population påvisades\.\s*Fyndet är en varningssignal och bör bedömas tillsammans med behandlingshistorik och HBV-DNA-utveckling\.\s*",
    ]
    for pattern_re in prefixes:
        text = re.sub(pattern_re, "", text, flags=re.I)

    if text:
        return text
    if level == "major":
        return f"Resistensmönster förenligt med {category.lower()}."
    if level == "minor":
        return "Resistensassocierat mönster i minoritetspopulation."
    if level == "caution":
        return "Ofullständigt eller kontextberoende resistensfynd."
    if level == "unknown":
        return "Bedömning begränsas av teknisk kvalitet."
    return "Inget etablerat resistensmönster påvisat."


def hbv_combination_rows_html(row):
    """Always-visible compact table of rule-based HBV resistance combinations."""
    css = {
        "major": "pattern-high",
        "minor": "pattern-minor",
        "caution": "pattern-caution",
        "none": "pattern-none",
        "unknown": "pattern-unknown",
    }
    rows = []
    for pattern in hbv_combination_interpretations(row):
        rows.append(
            "<tr>"
            f"<td><strong>{html_escape(pattern['category'])}</strong></td>"
            f"<td><span class='pattern-status {css.get(pattern.get('level'), 'pattern-none')}'>{html_escape(pattern['title'])}</span></td>"
            f"<td>{html_escape(pattern.get('mutations') or '—')}</td>"
            f"<td>{html_escape(_hbv_compact_pattern_interpretation(pattern))}</td>"
            "</tr>"
        )
    return "".join(rows)


def hbv_compact_interpretation_html(row):
    """Small overall resistance summary without the former large teaching cards."""
    patterns = hbv_combination_interpretations(row)
    major = [p for p in patterns if p.get("level") == "major"]
    minor = [p for p in patterns if p.get("level") == "minor"]
    caution = [p for p in patterns if p.get("level") == "caution"]
    unknown = [p for p in patterns if p.get("level") == "unknown"]

    if unknown:
        text = "Resistensbedömningen är inte fullständigt bedömbar på grund av otillräcklig teknisk kvalitet."
        cls = "compact-unknown"
    elif major:
        names = ", ".join(dict.fromkeys(p["category"] for p in major))
        text = f"Resistens påvisad: {names}."
        cls = "compact-high"
    elif minor:
        names = ", ".join(dict.fromkeys(p["category"] for p in minor))
        text = f"Resistensassocierat fynd i minoritetspopulation: {names}."
        cls = "compact-minor"
    elif caution:
        names = ", ".join(dict.fromkeys(p["category"] for p in caution))
        text = f"Möjligt nedsatt känslighet/ofullständigt resistensmönster: {names}."
        cls = "compact-caution"
    else:
        text = "Inget etablerat regelbaserat resistensmönster påvisat."
        cls = "compact-none"

    return f"<div class='compact-interpretation {cls}'><strong>Tolkning:</strong> {html_escape(text)}</div>"


def hbv_overall_summary(row, calls):
    rt_status = str(row.get("hbv_rt_assessment_status", "") or "").upper()
    if not rt_status:
        rt_status = "FULL" if row.get("hbv_quality") == "✓" else "NOT_EVALUABLE"
    rt_calls = [item for item in calls if item.get("gene") == "RT"]
    s_calls = [item for item in calls if item.get("gene") == "S"]
    if rt_status == "NOT_EVALUABLE":
        return "RT-resistensbedömning är inte möjlig på grund av otillräckligt lokalt kodondjup."
    parts = []
    if rt_calls:
        parts.append("Påvisade RT-panelmutationer: " + ", ".join(f"{item['mutation']} ({item['frequency']} %)" for item in rt_calls) + ".")
    elif rt_status == "LIMITED":
        parts.append("Inga resistensassocierade mutationer påvisades vid de bedömbara RT-positionerna.")
    else:
        parts.append("Inga resistensassocierade RT-panelmutationer påvisades över rapporteringströskeln.")
    if s_calls:
        parts.append("Påvisade ytantigen/escape-panelmutationer: " + ", ".join(f"{item['mutation']} ({item['frequency']} %)" for item in s_calls) + ".")
    if rt_status == "LIMITED":
        parts.append("RT-bedömningen är begränsad till positioner med minst " + str(HBV_MIN_REPORTABLE_COVERAGE) + " kompletta kodonläsningar.")
    return " ".join(parts)


def _hbv_escape_evidence(item):
    """Return mutation-specific curated escape labels and concise interpretation."""
    mutation = str(item.get("mutation", "") or "")
    evidence = HBV_S_ESCAPE_ASSOCIATIONS.get(mutation, {})
    labels = list(evidence.get("associations", ()) or ())
    if not labels:
        labels = ["Escape-associerad – specifik kategori saknas i referenstabellen"]
    interpretation = str(
        evidence.get("interpretation", "")
        or "Escape-associerad substitution; specifik fenotypisk effekt fastställs inte av denna analys."
    )
    return labels, interpretation


def _hbv_escape_associations(item):
    return _hbv_escape_evidence(item)[0]


def _hbv_escape_short_interpretation(item):
    return _hbv_escape_evidence(item)[1]

def hbv_panel_coverage_rows(row):
    rows = []
    for detail in row.get("hbv_mutation_details", []) or []:
        depth = int(detail.get("depth", 0) or 0)
        calls = detail.get("calls", []) or []
        nonsyn_calls = []
        for call in calls:
            alt = str(call.get("aa", "") or "")
            if alt and alt != detail.get("reference_aa"):
                nonsyn_calls.append(f"{alt} ({float(call.get('frequency', 0) or 0):.1f} %)")
        result = ", ".join(nonsyn_calls) or "Inget panelrelevant fynd"
        if depth >= HBV_MIN_COVERAGE:
            qc = "Godkänd"
        elif depth >= HBV_MIN_REPORTABLE_COVERAGE:
            qc = "Begränsad – rapporterbar"
        else:
            qc = "Otillräcklig"
        rows.append(
            "<tr>"
            f"<td>{html_escape(str(detail.get('position','')))}</td>"
            f"<td>{html_escape(str(detail.get('reference_aa','')))}</td>"
            f"<td>{depth}</td>"
            f"<td>{html_escape(result)}</td>"
            f"<td>{qc}</td>"
            "</tr>"
        )
    return "".join(rows) or "<tr><td colspan='5' class='empty'>Paneldata saknas.</td></tr>"


def hbv_rt_observed_variant_rows(row):
    """Show strong non-panel amino-acid observations without classifying them as resistance."""
    rows = []
    for detail in row.get("hbv_mutation_details", []) or []:
        if detail.get("gene") != "RT":
            continue
        depth = int(detail.get("depth", 0) or 0)
        if depth < HBV_MIN_REPORTABLE_COVERAGE:
            continue
        extra = detail.get("observed_nonpanel_variants", []) or []
        if not extra:
            continue
        profile = detail.get("aa_profile_above_threshold", []) or []
        profile_text = ", ".join(
            f"{x.get('aa','')} {float(x.get('frequency',0) or 0):.1f}% ({int(x.get('count',0) or 0)})"
            for x in profile
        )
        finding = "Blandpopulation" if len(profile) >= 2 else "Avvikande aminosyra"
        coverage = "Godkänd" if depth >= HBV_MIN_COVERAGE else "Begränsat lokalt läsdjup"
        rows.append(
            "<tr>"
            f"<td><strong>{html_escape(str(detail.get('position','')))}</strong></td>"
            f"<td>{html_escape(finding)}</td>"
            f"<td>{html_escape(profile_text or '—')}</td>"
            f"<td>{depth}</td>"
            f"<td>{html_escape(coverage)}</td>"
            "</tr>"
        )
    return "".join(rows) or "<tr><td colspan='5' class='empty'>Inga ytterligare aminosyravarianter över rapporteringströskeln påvisades vid RT-panelpositionerna.</td></tr>"


def hbv_recombination_window_rows_html(row):
    windows = row.get("hbv_recombination_windows", []) or []
    if not windows:
        status = str(row.get('hbv_recombination_status', '') or '')
        if status == 'Not indicated':
            return "<tr><td colspan='6' class='empty'>Fönsteranalys ej indicerad för denna högkonfidensgenotyp.</td></tr>"
        if status == 'Not evaluable':
            return "<tr><td colspan='6' class='empty'>Fönsteranalys ej bedömbar på grund av otillräcklig eller ojämnt fördelad sekvenstäckning.</td></tr>"
        return "<tr><td colspan='6' class='empty'>Ingen fönsterdata tillgänglig.</td></tr>"
    rows = []
    for item in windows:
        genotype = str(item.get("genotype", "") or "Ambiguous")
        rows.append(
            "<tr>"
            f"<td>{html_escape(str(item.get('start','')))}-{html_escape(str(item.get('end','')))}</td>"
            f"<td><strong>{html_escape(genotype)}</strong></td>"
            f"<td>{html_escape(str(item.get('identity','') or '—'))}</td>"
            f"<td>{html_escape(str(item.get('qcov','') or '—'))}</td>"
            f"<td>{html_escape(str(item.get('margin','') or '—'))}</td>"
            f"<td>{html_escape(str(item.get('accession','') or '—'))}</td>"
            "</tr>"
        )
    return "".join(rows)

def _hbv_assessment_scope(row):
    """Separate resistance evaluability from genotype evaluability."""
    def _f(key, default=0.0):
        try:
            return float(row.get(key, default) or default)
        except (TypeError, ValueError):
            return float(default)

    qc = str(row.get("qc_verdict", "") or "").upper()
    rt_status = str(row.get("hbv_rt_assessment_status", "") or "").upper()
    if not rt_status:
        rt_status = "FULL" if str(row.get("hbv_quality", "")) == "✓" else "NOT_EVALUABLE"
    resistance_ok = rt_status in {"FULL", "LIMITED"}
    resistance_complete = rt_status == "FULL"

    s1_ok = (_f("hbv_system1_unique_median_depth") >= 10 and
             _f("hbv_system1_unique_positions_ge_10x_percent") >= 70)
    s2_ok = (_f("hbv_system2_unique_median_depth") >= 10 and
             _f("hbv_system2_unique_positions_ge_10x_percent") >= 70)

    genotype = str(row.get("hbv_curated_genotype_display", "") or row.get("blast_genotype_display", "") or "").strip()
    confidence = str(row.get("hbv_curated_genotype_confidence", "") or row.get("blast_genotype_confidence", "") or "").strip().upper()
    geno_status = str(row.get("hbv_genotyping_status", "") or "").strip().upper()
    identity = _f("hbv_curated_reference_identity")
    query_cov = _f("hbv_curated_reference_query_coverage")
    informative_bases = max(0, int(_f("hbv_genotype_callable_bases", max(0, int(_f("consensus_length")) - int(_f("n_count"))))))
    callable_percent = _f(
        "hbv_genotype_callable_percent",
        (informative_bases / max(_f("consensus_length"), 1.0) * 100.0),
    )

    # Genotyping can remain reportable from a substantial partial genome when
    # one HBV system is well covered and the available sequence gives a clear
    # genotype match. Resistance and recombination keep stricter requirements.
    available_sequence_covered = query_cov >= max(35.0, callable_percent - 10.0)
    genotype_ok = (
        geno_status == "OK"
        and genotype
        and not genotype.lower().startswith("indeterminate")
        and genotype.lower() not in {"ej fastställd", "unknown"}
        and confidence in {"HIGH", "MODERATE", "MEDIUM", "LIMITED", "BEGRÄNSAD"}
        and identity >= 90.0
        and available_sequence_covered
        and informative_bases >= 1000
        and (s1_ok or s2_ok)
    )

    if resistance_complete and genotype_ok:
        text, css = "Bedömbar för resistens och genotypning.", "status-good"
    elif resistance_complete:
        text, css = "Bedömbar för resistens; genotyp ej säkert fastställd.", "status-medium"
    elif resistance_ok and genotype_ok:
        text, css = "Resistensbedömning möjlig för täckta positioner; genotypning bedömbar.", "status-medium"
    elif resistance_ok:
        text, css = "Resistensbedömning möjlig för täckta positioner; genotyp ej säkert fastställd.", "status-medium"
    elif genotype_ok:
        text, css = "Ej bedömbar för resistens; bedömbar endast för genotypning.", "status-medium"
    else:
        text, css = "Ej bedömbar för resistens eller säker genotypning.", "status-unknown"

    return {
        "resistance_ok": resistance_ok,
        "resistance_complete": resistance_complete,
        "rt_status": rt_status,
        "genotype_ok": genotype_ok,
        "text": text,
        "css": css,
        "system1_genotype_ok": s1_ok,
        "system2_genotype_ok": s2_ok,
    }


def _hbv_confidence_sv(value):
    text = str(value or "").strip()
    return {
        "HIGH": "Hög",
        "MODERATE": "Måttlig",
        "MEDIUM": "Måttlig",
        "LIMITED": "Begränsad",
        "LOW": "Låg",
    }.get(text.upper(), text or "—")


def _hbv_recombination_status_sv(value):
    text = str(value or "").strip()
    return {
        "NOT RUN": "Ej körd",
        "NOT INDICATED": "Ej indicerad",
        "NOT EVALUABLE": "Ej bedömbar",
        "INDETERMINATE": "Ej bedömbar",
        "NO SIGNIFICANT INTERGENOTYPIC MOSAIC DETECTED": "Ingen signifikant intergenotypisk mosaik påvisad",
        "POSSIBLE INTERGENOTYPIC RECOMBINATION": "Möjligt intergenotypiskt rekombinant mönster",
        "GENOTYPE-I-COMPATIBLE MOSAIC": "Mönster förenligt med genotyp I",
        "KNOWN RECOMBINANT LINEAGE (GENOTYPE I)": "Känd rekombinant linje (genotyp I)",
    }.get(text.upper(), text or "Ej körd")


def write_hbv_html_report(row, report_path, run_tag):
    calls = hbv_detected_calls(row)
    rt_calls = [x for x in calls if x["gene"] == "RT"]
    s_calls = [x for x in calls if x["gene"] == "S"]
    logo_uri = image_file_to_data_uri(find_hbv_report_logo(row, report_path))
    logo_html = f"<img class='report-logo' src='{logo_uri}' alt='Logotyp'>" if logo_uri else ""
    assessment = _hbv_assessment_scope(row)
    if assessment["resistance_ok"]:
        summary_text = hbv_overall_summary(row, calls)
    elif assessment["genotype_ok"]:
        summary_text = "Resistensbedömning ej möjlig med tillräcklig säkerhet. HBV-genotyp kan rapporteras från tillgänglig sekvensdata."
    else:
        summary_text = "Varken resistensbedömning eller säker HBV-genotypning är möjlig med tillgänglig sekvensdata."

    def status_class():
        if not assessment["resistance_ok"]:
            return "status-unknown"
        if any(float(x.get("frequency", 0) or 0) >= 50 for x in calls):
            return "status-high"
        if calls:
            return "status-medium"
        return "status-good"

    medication_rows = []
    swedish_drugs = {
        "Lamivudine": "Lamivudin",
        "Adefovir": "Adefovir",
        "Entecavir": "Entekavir",
        "Tenofovir": "Tenofovir",
    }
    combination_patterns = hbv_combination_interpretations(row)
    for drug in HBV_MEDICATIONS:
        status, mutations, level = hbv_status_for_medication(row, drug)
        drug_sv = swedish_drugs.get(drug, drug)
        related_patterns = [p for p in combination_patterns if str(p.get("category", "")).lower() == drug_sv.lower()]
        if any(p.get("level") == "major" for p in related_patterns):
            short_status = "Resistens"
            short_interp = _hbv_compact_pattern_interpretation(next(p for p in related_patterns if p.get("level") == "major"))
            level = "major"
        elif any(p.get("level") == "minor" for p in related_patterns):
            short_status = "Resistensassocierat fynd – minoritetspopulation"
            short_interp = _hbv_compact_pattern_interpretation(next(p for p in related_patterns if p.get("level") == "minor"))
            level = "minor"
        elif any(p.get("level") == "caution" for p in related_patterns):
            short_status = "Möjligt nedsatt känslighet / ofullständigt mönster"
            short_interp = _hbv_compact_pattern_interpretation(next(p for p in related_patterns if p.get("level") == "caution"))
            level = "caution"
        elif level == "major":
            short_status = "Resistensassocierad mutation"
            short_interp = "Se påvisade mutationer och regelbaserad bedömning."
        elif level == "minor":
            short_status = "Resistensassocierad mutation – minoritetspopulation"
            short_interp = "Fyndet förekommer under 50 % och ska bedömas i kliniskt sammanhang."
        elif level == "caution":
            short_status = "Försiktighetsfynd"
            short_interp = "Association finns, men fyndet fastställer inte ensamt klinisk resistens."
        else:
            if assessment.get("resistance_complete"):
                short_status = "Ingen resistens påvisad i panelen"
                short_interp = "Inga associerade panelmutationer påvisade."
            else:
                short_status = "Ingen associerad mutation påvisad vid bedömbara positioner"
                short_interp = "Bedömningen är begränsad av lokal RT-täckning; ej bedömbara positioner kan inte uteslutas."

        medication_rows.append(
            f"<tr><td><strong>{html_escape(drug_sv)}</strong></td>"
            f"<td><span class='drug-status drug-{level}'>{html_escape(short_status)}</span></td>"
            f"<td>{html_escape(mutations or '—')}</td>"
            f"<td>{html_escape(short_interp)}</td></tr>"
        )

    def rt_mutation_rows(items, empty_text):
        if not items:
            return f"<tr><td colspan='7' class='empty'>{html_escape(empty_text)}</td></tr>"
        rows_html = []
        for item in items:
            meds = ", ".join(
                f"{swedish_drugs.get(drug, drug)}: {hbv_swedish_evidence_label(evidence)}"
                for drug, evidence in item.get("medications", {}).items()
            ) or "—"
            flag = item.get("flag", "")
            flag_sv = {
                "PRIORITY": "Prioriterat fynd",
                "ESTABLISHED": "Etablerat fynd",
                "CAUTION": "Kontextberoende/försiktighet",
                "ESCAPE": "Escape-associerat fynd",
            }.get(flag, flag or "—")
            rows_html.append(
                "<tr>"
                f"<td><strong>{html_escape(str(item['mutation']))}</strong></td>"
                f"<td>{item['frequency']} %</td>"
                f"<td>{html_escape(item['population'])}</td>"
                f"<td>{item['count']} / {item['depth']}</td>"
                f"<td>{html_escape('Begränsat lokalt läsdjup' if item.get('coverage_status') == 'limited' else 'Godkänd')}</td>"
                f"<td>{html_escape(flag_sv)}</td>"
                f"<td>{html_escape(meds)}</td>"
                "</tr>"
            )
        return "".join(rows_html)

    def escape_mutation_rows(items, empty_text):
        if not items:
            return f"<tr><td colspan='7' class='empty'>{html_escape(empty_text)}</td></tr>"
        rows_html = []
        for item in items:
            associations = _hbv_escape_associations(item)
            association_html = " ".join(f"<span class='escape-tag'>{html_escape(label)}</span>" for label in associations)
            rows_html.append(
                "<tr>"
                f"<td><strong>{html_escape(str(item['mutation']))}</strong></td>"
                f"<td>{item['frequency']} %</td>"
                f"<td>{html_escape(item['population'])}</td>"
                f"<td>{item['count']} / {item['depth']}</td>"
                f"<td>{association_html}</td>"
                f"<td>{html_escape(str(item.get('category','') or '—'))}</td>"
                f"<td>{html_escape(_hbv_escape_short_interpretation(item))}</td>"
                "</tr>"
            )
        return "".join(rows_html)

    generated_date = datetime.now().strftime("%Y-%m-%d %H:%M")
    qc_verdict = str(row.get("qc_verdict", "") or "Ej bedömd").upper()
    qc_badge_class = {"HIGH": "qc-high", "MEDIUM": "qc-medium", "LOW": "qc-low"}.get(qc_verdict, "qc-neutral")

    # System-specific unique-region breadth is what the HBV assay QC uses to
    # demonstrate that both amplification systems are represented.
    s1_cov = f"{row.get('hbv_system1_unique_positions_ge_10x_percent', 0)} %"
    s2_cov = f"{row.get('hbv_system2_unique_positions_ge_10x_percent', 0)} %"

    if assessment["resistance_ok"]:
        resistance_sections_html = f"""
<h2>Resistens- och kombinationsbedömning</h2>
{hbv_compact_interpretation_html(row)}
<table><thead><tr><th>Läkemedel/mönster</th><th>Bedömning</th><th>Mutationer i provet</th><th>Kort tolkning</th></tr></thead>
<tbody>{hbv_combination_rows_html(row)}</tbody></table>
<div class='small' style='margin-top:6px'><strong>Obs:</strong> Mutationerna analyseras på läsnivå per kodon men är inte fasade över hela RT-regionen. Ett kombinationsfynd visar att komponenterna finns i samma prov, men inte säkert på samma virusgenom.</div>

<h2>Läkemedelsöversikt</h2>
<table><thead><tr><th>Läkemedel</th><th>Bedömning</th><th>Påvisade mutationer</th><th>Kort tolkning</th></tr></thead>
<tbody>{''.join(medication_rows)}</tbody></table>

<h2>RT-resistensassocierade mutationer</h2>
<table><thead><tr><th>Mutation</th><th>Andel</th><th>Population</th><th>Mutantläsningar / djup</th><th>Lokal täckning</th><th>Typ</th><th>Läkemedelsassociation</th></tr></thead>
<tbody>{rt_mutation_rows(rt_calls, 'Inga RT-panelmutationer påvisades vid bedömbara positioner.')}</tbody></table>
<p class='small'>Mutationer med 30–99 kompletta kodonläsningar kan rapporteras, men markeras som begränsat lokalt läsdjup. Positioner med färre än {HBV_MIN_REPORTABLE_COVERAGE} kompletta kodonläsningar bedöms inte.</p>

<h2>Övriga aminosyravarianter vid RT-panelpositioner</h2>
<table><thead><tr><th>Position</th><th>Observation</th><th>Aminosyraprofil ≥{HBV_MIN_VARIANT_FREQUENCY * 100:.1f}%</th><th>Kodontäckning</th><th>Lokal täckning</th></tr></thead>
<tbody>{hbv_rt_observed_variant_rows(row)}</tbody></table>
<p class='small'>Denna tabell visar starka aminosyraobservationer som inte ingår i den aktuella resistensdefinitionen för positionen. De rapporteras som sekvensfynd och klassificeras inte automatiskt som resistensmutationer.</p>

<h2>Ytantigen- och escape-associerade mutationer</h2>
<table><thead><tr><th>Mutation</th><th>Andel</th><th>Population</th><th>Mutantläsningar / djup</th><th>Associerad escape</th><th>Panelklass</th><th>Kort tolkning</th></tr></thead>
<tbody>{escape_mutation_rows(s_calls, 'Inga escape-associerade panelmutationer påvisades.')}</tbody></table>
<p class='small'>Escape-etiketterna är mutationsspecifika och litteraturbaserade. Terapiassocierad anger en HBsAg-förändring som kan selekteras via den överlappande polymeras/S-regionen under nukleos(t)idanalog behandling; antiviral resistens bedöms separat i RT-tabellen. Etiketterna innebär inte i sig bevisad fenotypisk escape hos patienten.</p>
"""
    else:
        resistance_sections_html = f"""
<h2>Resistensbedömning</h2>
<div class='status status-unknown'><strong>{html_escape(assessment['text'])}</strong></div>
<p class='small'>RT-resistensbedömning visas inte när ingen RT-panelposition når minsta rapporterbara lokala kodondjup. Genotypning bedöms separat och kan därför fortfarande vara rapporterbar.</p>
"""

    html = f"""<!DOCTYPE html>
<html lang='sv'>
<head>
<meta charset='UTF-8'>
<meta name='viewport' content='width=device-width, initial-scale=1'>
<title>HBV-resistensrapport – {html_escape(str(row.get('sample', '')))}</title>
<style>
:root {{ --blue:#dbeeff; --blue-dark:#185b8d; --line:#c9d3dc; --text:#1e2a32; --muted:#5d6b75; }}
* {{ box-sizing:border-box; }}
body {{ font-family:Arial,Helvetica,sans-serif; color:var(--text); margin:0; background:#eef2f5; }}
.report {{ max-width:1050px; margin:24px auto; background:#fff; padding:30px 38px; box-shadow:0 3px 18px rgba(0,0,0,.12); }}
h1 {{ margin:0 0 5px; font-size:28px; color:#123e5f; }}
.header {{ display:flex; justify-content:space-between; align-items:flex-start; gap:24px; margin-bottom:16px; }}
.title-block {{ flex:1; }}
.report-logo {{ max-width:220px; max-height:105px; object-fit:contain; }}
.subtitle {{ color:var(--muted); margin-bottom:0; }}
.summary {{ background:#f4f8fb; border-left:5px solid var(--blue-dark); padding:10px 13px; margin:12px 0; }}
h2 {{ font-size:18px; margin:22px 0 8px; border-bottom:2px solid var(--blue-dark); padding-bottom:5px; }}
table {{ width:100%; border-collapse:collapse; font-size:13px; }}
th {{ background:var(--blue); text-align:left; }}
th,td {{ border:1px solid var(--line); padding:7px 8px; vertical-align:top; }}
.meta td:first-child {{ width:28%; background:#f4f8fb; font-weight:bold; }}
.status {{ padding:10px 13px; border-radius:5px; font-weight:bold; margin:10px 0 12px; }}
.status-good {{ background:#dff3e3; border-left:6px solid #3e9d54; }}
.status-medium {{ background:#fff4c7; border-left:6px solid #d0a600; }}
.status-high {{ background:#ffe1d6; border-left:6px solid #d86032; }}
.status-unknown {{ background:#e3e6e8; border-left:6px solid #6f777d; }}
.qc-badge {{ display:inline-block; min-width:82px; text-align:center; padding:4px 10px; border-radius:999px; font-weight:700; }}
.qc-high {{ background:#dff3e3; color:#28683a; }} .qc-medium {{ background:#fff1b8; color:#725b00; }} .qc-low {{ background:#f8d7cc; color:#8d2e17; }} .qc-neutral {{ background:#e3e6e8; color:#50585d; }}
.drug-status {{ font-weight:bold; }}
.drug-none {{ color:#2c7b3e; }} .drug-minor {{ color:#a17900; }} .drug-major {{ color:#b23b20; }} .drug-caution {{ color:#8a5a00; }}
.pattern-status {{ font-weight:bold; }} .pattern-high {{ color:#b23b20; }} .pattern-minor {{ color:#8a6a00; }} .pattern-caution {{ color:#8a5a00; }} .pattern-none {{ color:#2c7b3e; }} .pattern-unknown {{ color:#6f777d; }}
.compact-interpretation {{ margin:8px 0 10px; padding:8px 10px; border-left:4px solid #6f777d; background:#f7f8f9; font-size:13px; }}
.compact-high {{ border-left-color:#c94d28; }} .compact-minor {{ border-left-color:#d2a600; }} .compact-caution {{ border-left-color:#c99b00; }} .compact-none {{ border-left-color:#438f55; }}
.escape-tag {{ display:inline-block; margin:1px 3px 1px 0; padding:2px 6px; border-radius:999px; background:#eef3f7; border:1px solid #d5dee5; font-size:11px; white-space:nowrap; }}
.empty {{ color:var(--muted); font-style:italic; }}
.notice {{ background:#f5f5f5; border-left:5px solid #777; padding:10px 12px; margin-top:22px; font-size:12px; }}
.small {{ font-size:12px; color:var(--muted); }}
.report-footer {{ margin-top:24px; padding-top:10px; border-top:1px solid var(--line); color:var(--muted); font-size:10px; line-height:1.35; text-align:center; }}
@media print {{ body {{ background:#fff; }} .report {{ box-shadow:none; margin:0; max-width:none; padding:20px 24px; }} .report-logo {{ max-height:85px; }} h2 {{ break-after:avoid; }} table {{ break-inside:auto; }} tr {{ break-inside:avoid; }} }}
</style>
</head>
<body><main class='report'>
<div class='header'><div class='title-block'>
<h1>Hepatit B-virus – resistens och immunologisk escape</h1>
<div class='subtitle'>Läsnivåbaserad mutationsanalys</div>
</div>{logo_html}</div>
<table class='meta'>
<tr><td>Prov-ID</td><td>{html_escape(str(row.get('sample','')))}</td></tr>
<tr><td>HBV-genotyp</td><td><strong>{html_escape(str(row.get('hbv_curated_genotype_display','') or row.get('blast_genotype_display','') or 'Ej fastställd'))}</strong> · konfidens {html_escape(_hbv_confidence_sv(row.get('hbv_curated_genotype_confidence','') or row.get('blast_genotype_confidence','')))}</td></tr>
<tr><td>Genotypreferens</td><td><strong>{html_escape(str(row.get('hbv_curated_reference_accession','') or 'Ej tillgänglig'))}</strong> · identitet {html_escape(str(row.get('hbv_curated_reference_identity','') or '—'))} % · query coverage {html_escape(str(row.get('hbv_curated_reference_query_coverage','') or '—'))} %</td></tr>
<tr><td>Rekombinationsscreening</td><td><strong>{html_escape(_hbv_recombination_status_sv(row.get('hbv_recombination_status','')))}</strong>{(' · ' + html_escape(str(row.get('hbv_recombination_components','')))) if row.get('hbv_recombination_components') else ''}</td></tr>
<tr><td>BLAST toppträff (lokal NCBI-databas)</td><td><strong>{html_escape(str(row.get('blast_top_accession','') or 'Ej tillgänglig'))}</strong> · {html_escape(str(row.get('blast_top_organism','')))} · identitet {html_escape(str(row.get('blast_top_identity','')))} % · query coverage {html_escape(str(row.get('blast_top_query_coverage','')))} %</td></tr>
<tr><td>Bästa mappningsreferens</td><td>{html_escape(str(row.get('best_reference','')))}</td></tr>
<tr><td>Analysbedömning</td><td><strong>{html_escape(assessment['text'])}</strong></td></tr>
<tr><td>QC-status</td><td><span class='qc-badge {qc_badge_class}'>{html_escape(qc_verdict)}</span></td></tr>
<tr><td>System 1 – unik region ≥10×</td><td>{html_escape(s1_cov)}</td></tr>
<tr><td>System 2 – unik region ≥10×</td><td>{html_escape(s2_cov)}</td></tr>
<tr><td>Totalt antal mappade läsningar</td><td>{row.get('mapped_reads', 0)}</td></tr>
<tr><td>Genomtäckning</td><td>{row.get('target_or_genome_coverage_percent', 0)} %</td></tr>
<tr><td>Mediandjup</td><td>{row.get('median_depth', 0)}×</td></tr>
<tr><td>Rapportdatum</td><td>{generated_date}</td></tr>
<tr><td>Körning</td><td>{html_escape(str(run_tag))}</td></tr>
<tr><td>Rapporteringströskel</td><td>≥ {HBV_MIN_VARIANT_FREQUENCY * 100:.1f} %</td></tr>
<tr><td>Lokalt kodondjup</td><td>≥{HBV_MIN_REPORTABLE_COVERAGE} läsningar för rapportering; ≥{HBV_MIN_COVERAGE} läsningar föredragen nivå</td></tr>
</table>
<div class='summary'><strong>Sammanfattning:</strong> {html_escape(summary_text)}</div>
<div class='status {status_class()}'>{html_escape(hbv_swedish_interpretation(row.get('hbv_rt_interpretation', row.get('hbv_interpretation',''))) if assessment['resistance_ok'] else assessment['text'])}</div>

{resistance_sections_html}

<h2>HBV-genotypning och rekombinationsscreening</h2>
<table class='meta'>
<tr><td>Metod</td><td>{html_escape(str(row.get('hbv_genotyping_method','') or '—'))}</td></tr>
<tr><td>Genotypning</td><td>{html_escape(str(row.get('hbv_genotyping_status','') or '—'))}</td></tr>
<tr><td>Genotypstöd</td><td>{html_escape(str(row.get('hbv_curated_genotype_support','') or '—'))}</td></tr>
<tr><td>Analyserbar sekvens</td><td>{html_escape(str(row.get('hbv_genotype_callable_bases','') or '—'))} baser · {html_escape(str(row.get('hbv_genotype_callable_percent','') or '—'))} % · största N-gap {html_escape(str(row.get('hbv_genotype_largest_n_gap','') or '—'))} nt</td></tr>
<tr><td>Rekombinationsbedömning</td><td>{html_escape(str(row.get('hbv_recombination_summary','') or '—'))}</td></tr>
<tr><td>Fönster / steg</td><td>{html_escape(str(row.get('hbv_recombination_window_size','') or '—'))} nt / {html_escape(str(row.get('hbv_recombination_step_size','') or '—'))} nt</td></tr>
<tr><td>Genotypdatabas</td><td class='small'>{html_escape(str(row.get('hbv_genotyping_database','') or '—'))}</td></tr>
</table>
<table><thead><tr><th>Genomposition</th><th>Parental-lik genotyp</th><th>Identitet %</th><th>Query coverage %</th><th>Identitetsmarginal</th><th>Närmaste referens</th></tr></thead>
<tbody>{hbv_recombination_window_rows_html(row)}</tbody></table>
<p class='small'>Rekombinationsscreening kräver tillräckligt komplett och sammanhängande genomdata. Genotyp I karakteriseras när genomkvaliteten medger det; genotyp A–H fönsteranalyseras endast när helgenomsresultatet motiverar det. Ett möjligt rekombinant mönster kräver en uthållig alternativ genotypkomponent över flera överlappande fönster. Ovanliga eller nya mönster bör bekräftas med fylogenetisk/rekombinationsspecifik metod.</p>

<h2>Teknisk kvalitet – detaljer</h2>
<table class='meta'>
<tr><td>Kvalitetsstatus</td><td>{html_escape(str(row.get('hbv_quality','')))}</td></tr>
<tr><td>Övergripande QC</td><td>{html_escape(str(row.get('qc_verdict','')))}</td></tr>
<tr><td>QC-kommentar</td><td>{html_escape(str(row.get('qc_reason','') or '—'))}</td></tr>
<tr><td>System 1 unik region ≥10×</td><td>{row.get('hbv_system1_unique_positions_ge_10x_percent', 0)} %</td></tr>
<tr><td>System 1 unik region ≥30×</td><td>{row.get('hbv_system1_unique_positions_ge_30x_percent', 0)} %</td></tr>
<tr><td>System 1 unikt mediandjup</td><td>{row.get('hbv_system1_unique_median_depth', 0)}×</td></tr>
<tr><td>System 2 unik region ≥10×</td><td>{row.get('hbv_system2_unique_positions_ge_10x_percent', 0)} %</td></tr>
<tr><td>System 2 unik region ≥30×</td><td>{row.get('hbv_system2_unique_positions_ge_30x_percent', 0)} %</td></tr>
<tr><td>System 2 unikt mediandjup</td><td>{row.get('hbv_system2_unique_median_depth', 0)}×</td></tr>
<tr><td>RT-positioner &lt;{HBV_MIN_COVERAGE}×</td><td>{html_escape(str(row.get('hbv_rt_low_coverage_positions',' ') or '—'))}</td></tr>
<tr><td>RT-positioner &lt;{HBV_MIN_REPORTABLE_COVERAGE}× (ej bedömbara)</td><td>{html_escape(str(row.get('hbv_rt_unassessable_positions',' ') or '—'))}</td></tr>
<tr><td>S-positioner &lt;{HBV_MIN_COVERAGE}×</td><td>{html_escape(str(row.get('hbv_s_low_coverage_positions',' ') or '—'))}</td></tr>
<tr><td>S-positioner &lt;{HBV_MIN_REPORTABLE_COVERAGE}× (ej bedömbara)</td><td>{html_escape(str(row.get('hbv_s_unassessable_positions',' ') or '—'))}</td></tr>
<tr><td>Alla panelpositioner &lt;{HBV_MIN_COVERAGE}×</td><td>{html_escape(str(row.get('hbv_low_coverage_positions',' ') or '—'))}</td></tr>
<tr><td>Medeldjup</td><td>{row.get('mean_depth', 0)}×</td></tr>
<tr><td>Minimidjup</td><td>{row.get('min_depth', 0)}×</td></tr>
<tr><td>Positioner ≥10×</td><td>{row.get('positions_ge_10x_percent', 0)} %</td></tr>
<tr><td>Positioner ≥30×</td><td>{row.get('positions_ge_30x_percent', 0)} %</td></tr>
<tr><td>Konsensuslängd</td><td>{row.get('consensus_length', 0)} nt</td></tr>
<tr><td>N-baser</td><td>{row.get('n_count', 0)} ({row.get('n_percent', 0)} %)</td></tr>
<tr><td>Konsensusåtervinning</td><td>{row.get('consensus_recovery_percent', 0)} %</td></tr>
<tr><td>RT-koordinat</td><td>start {html_escape(str(row.get('hbv_rt_codon1_start','')))}, metod {html_escape(str(row.get('hbv_coordinate_method','metadata')))}</td></tr>
<tr><td>S-koordinat</td><td>start {html_escape(str(row.get('hbv_s_codon1_start','')))}, metod {html_escape(str(row.get('hbv_coordinate_method','metadata')))}</td></tr>
<tr><td>Medaka-modell</td><td>{html_escape(str(row.get('medaka_model','')))}</td></tr>
<tr><td>FASTQ-fil</td><td class='small'>{html_escape(str(row.get('input_fastq','')))}</td></tr>
<tr><td>Konsensusfil</td><td class='small'>{html_escape(str(row.get('consensus_file','')))}</td></tr>
<tr><td>BAM-fil</td><td class='small'>{html_escape(str(row.get('bam_file','')))}</td></tr>
<tr><td>Referensfil</td><td class='small'>{html_escape(str(row.get('reference_file','')))}</td></tr>
<tr><td>BLAST-databas</td><td class='small'>{html_escape(str(row.get('blast_database','')))}</td></tr>
<tr><td>BLAST-databasdatum</td><td>{html_escape(str(row.get('blast_database_date','')))}</td></tr>
<tr><td>BLAST-status</td><td>{html_escape(str(row.get('blast_status','')))} {html_escape(str(row.get('blast_message','')))}</td></tr>
<tr><td>Scriptversion</td><td>{html_escape(str(row.get('hbv_mutation_script_version','')))}</td></tr>
</table>

<h2>Paneltäckning och positionskontroll</h2>
<table><thead><tr><th>Position</th><th>Referens-AA</th><th>Kodontäckning</th><th>Panelrelevant resultat</th><th>QC</th></tr></thead>
<tbody>{hbv_panel_coverage_rows(row)}</tbody></table>
<div class='notice'><strong>Viktig begränsning:</strong> Resultatet är en litteraturbaserad screening och är inte en fristående klinisk resistensbedömning. Regelbaserade kombinationsfynd ska bedömas tillsammans med variantandelar, behandlingshistorik, HBV-DNA-utveckling, genotyp och kliniska data. Tenofovirassociationer är försiktighetsfynd och inte ensamma bevis för klinisk resistens.</div>

<footer class='report-footer'>
Vid frågor, tekniska problem eller oväntade analysresultat, kontakta skriptets författare Joakim Bedner Stenbäck (joakim.stenback@vgregion.se).
</footer>

</main></body></html>"""
    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write(html)


def _adenovirus_type_from_title(title):
    """Return an explicitly labelled HAdV type from a BLAST subject title when available."""
    text = str(title or "")
    patterns = [
        r"\b(?:human\s+)?adenovirus\s+(?:type\s*)?([A-G]?\s*\d{1,3})\b",
        r"\bHAdV[-_ ]?([A-G]?\d{1,3})\b",
        r"\btype\s+([A-G]?\d{1,3})\b",
    ]
    for pattern in patterns:
        match = re.search(pattern, text, re.I)
        if match:
            return re.sub(r"\s+", "", match.group(1)).upper()
    return ""


def write_adenovirus_html_report(row, report_path, run_tag):
    """Write a compact Adenovirus BLAST/QC report without minor-variant analysis."""
    ensure_dir(os.path.dirname(report_path), must_be_writable=True)

    top_accession = str(row.get("blast_top_accession", "") or "")
    top_organism = str(row.get("blast_top_organism", "") or "")
    top_title = str(row.get("blast_top_title", "") or "")
    top_identity = row.get("blast_top_identity", "")
    top_qcov = row.get("blast_top_query_coverage", "")
    adeno_type = _adenovirus_type_from_title(top_title) or "Ej fastställd"

    hit_rows = []
    for hit in (row.get("blast_hits", []) or [])[:5]:
        hit_title = str(hit.get("stitle", "") or "")
        hit_type = _adenovirus_type_from_title(hit_title) or "—"
        hit_rows.append(
            "<tr>"
            f"<td>{html_escape(str(hit.get('accession', '') or ''))}</td>"
            f"<td>{html_escape(hit_type)}</td>"
            f"<td>{html_escape(str(hit.get('pident', '') or ''))}</td>"
            f"<td>{html_escape(str(hit.get('qcovs', '') or ''))}</td>"
            f"<td>{html_escape(hit_title)}</td>"
            "</tr>"
        )
    if not hit_rows:
        hit_rows.append("<tr><td colspan='5'>Inga BLAST-träffar tillgängliga.</td></tr>")

    qc = str(row.get("qc_verdict", "") or "")
    qc_class = {"HIGH": "good", "MEDIUM": "medium", "LOW": "low"}.get(qc.upper(), "neutral")
    generated_date = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Use the same embedded report logo as the HCV/HBV reports.
    logo_uri = image_file_to_data_uri(find_hbv_report_logo(row, report_path))
    logo_html = f"<img class='report-logo' src='{logo_uri}' alt='Logotyp'>" if logo_uri else ""

    html = f"""<!DOCTYPE html>
<html lang='sv'><head><meta charset='utf-8'><title>Adenovirusrapport - {html_escape(str(row.get('sample','')))}</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Arial,sans-serif;background:#f5f7fa;color:#17202a;margin:0}}
.report{{max-width:980px;margin:28px auto;background:white;padding:30px 34px;border-radius:12px;box-shadow:0 2px 14px rgba(0,0,0,.08)}}
.header{{display:flex;justify-content:space-between;align-items:flex-start;gap:24px;margin-bottom:20px}}
.title-block{{flex:1}} .report-logo{{max-width:220px;max-height:105px;object-fit:contain}}
h1{{margin:0 0 4px;font-size:28px}} .subtitle{{color:#5d6d7e;margin-bottom:0}}
h2{{margin-top:28px;border-bottom:1px solid #dfe6e9;padding-bottom:7px;font-size:19px}}
table{{border-collapse:collapse;width:100%;margin:10px 0 18px}} th,td{{padding:8px 10px;border:1px solid #dfe6e9;text-align:left;vertical-align:top}} th{{background:#eef3f7}}
.meta td:first-child{{width:245px;font-weight:600;background:#f8fafb}} .small{{font-size:12px;word-break:break-all}}
.status{{display:inline-block;padding:6px 12px;border-radius:18px;font-weight:700}} .good{{background:#d5f5e3}} .medium{{background:#fcf3cf}} .low{{background:#fadbd8}} .neutral{{background:#eaecee}}
.note{{background:#f8f9f9;padding:12px 14px;border-left:4px solid #85929e;margin-top:18px}}
footer{{margin-top:28px;color:#7b7d7d;font-size:12px}}
@media print {{ body {{ background:#fff; }} .report {{ box-shadow:none; margin:0; max-width:none; }} .report-logo {{ max-height:85px; }} }}
</style></head><body><main class='report'>
<div class='header'><div class='title-block'>
<h1>Adenovirus – BLAST och teknisk kvalitetsrapport</h1>
<div class='subtitle'>Automatisk analys av Medaka-konsensus</div>
</div>{logo_html}</div>
<table class='meta'>
<tr><td>Prov-ID</td><td><strong>{html_escape(str(row.get('sample','')))}</strong></td></tr>
<tr><td>Adenovirustyp (BLAST)</td><td><strong>{html_escape(adeno_type)}</strong></td></tr>
<tr><td>BLAST toppträff</td><td><strong>{html_escape(top_accession or 'Ej tillgänglig')}</strong> · {html_escape(top_organism)} · identitet {html_escape(str(top_identity))} % · query coverage {html_escape(str(top_qcov))} %</td></tr>
<tr><td>BLAST titel</td><td>{html_escape(top_title or '—')}</td></tr>
<tr><td>Bästa mappningsreferens</td><td>{html_escape(str(row.get('best_reference','')))}</td></tr>
<tr><td>Rapportdatum</td><td>{generated_date}</td></tr>
<tr><td>Körning</td><td>{html_escape(str(run_tag))}</td></tr>
</table>

<h2>Teknisk kvalitet</h2>
<p><span class='status {qc_class}'>{html_escape(qc or 'Ej bedömd')}</span></p>
<table class='meta'>
<tr><td>QC-kommentar</td><td>{html_escape(str(row.get('qc_reason','') or '—'))}</td></tr>
<tr><td>Mappade läsningar</td><td>{html_escape(str(row.get('mapped_reads',0)))}</td></tr>
<tr><td>Mål-/genomtäckning</td><td>{html_escape(str(row.get('target_or_genome_coverage_percent',0)))} %</td></tr>
<tr><td>Medeldjup</td><td>{html_escape(str(row.get('mean_depth',0)))}×</td></tr>
<tr><td>Mediandjup</td><td>{html_escape(str(row.get('median_depth',0)))}×</td></tr>
<tr><td>Positioner ≥10×</td><td>{html_escape(str(row.get('positions_ge_10x_percent',0)))} %</td></tr>
<tr><td>Positioner ≥30×</td><td>{html_escape(str(row.get('positions_ge_30x_percent',0)))} %</td></tr>
<tr><td>Konsensuslängd</td><td>{html_escape(str(row.get('consensus_length',0)))} nt</td></tr>
<tr><td>N-baser</td><td>{html_escape(str(row.get('n_count',0)))} ({html_escape(str(row.get('n_percent',0)))} %)</td></tr>
<tr><td>Konsensusåtervinning</td><td>{html_escape(str(row.get('consensus_recovery_percent',0)))} %</td></tr>
</table>

<h2>De fem bästa BLAST-träffarna</h2>
<table><thead><tr><th>Accession</th><th>Typ</th><th>Identitet %</th><th>Query coverage %</th><th>NCBI-titel</th></tr></thead><tbody>{''.join(hit_rows)}</tbody></table>

<h2>Analysinformation</h2>
<table class='meta'>
<tr><td>BLAST-status</td><td>{html_escape(str(row.get('blast_status','')))} {html_escape(str(row.get('blast_message','')))}</td></tr>
<tr><td>BLAST-databas</td><td class='small'>{html_escape(str(row.get('blast_database','')))}</td></tr>
<tr><td>BLAST-databasdatum</td><td>{html_escape(str(row.get('blast_database_date','')))}</td></tr>
<tr><td>Medaka-modell</td><td>{html_escape(str(row.get('medaka_model','')))}</td></tr>
<tr><td>FASTQ-fil</td><td class='small'>{html_escape(str(row.get('input_fastq','')))}</td></tr>
<tr><td>Konsensusfil</td><td class='small'>{html_escape(str(row.get('consensus_file','')))}</td></tr>
<tr><td>BAM-fil</td><td class='small'>{html_escape(str(row.get('bam_file','')))}</td></tr>
<tr><td>Referensfil</td><td class='small'>{html_escape(str(row.get('reference_file','')))}</td></tr>
</table>
<div class='note'>Adenovirustyp visas endast när typen framgår explicit av BLAST-träffens NCBI-titel. Rapporten gör ingen separat rekombinations- eller minorvariantanalys.</div>
<footer>Aurora Adenovirus BLAST/QC report</footer>
</main></body></html>"""
    with open(report_path, "w", encoding="utf-8") as handle:
        handle.write(html)


def generate_adenovirus_html_reports(rows, output_folder, run_tag):
    """Generate one compact HTML BLAST/QC report per Adenovirus sample."""
    html_dir = os.path.join(run_samples_folder(output_folder, run_tag), "Adenovirus_HTML_reports")
    ensure_dir(html_dir, must_be_writable=True)
    for row in rows:
        sample = safe_report_filename(str(row.get("sample", "sample") or "sample"))
        report_path = os.path.join(html_dir, f"{sample}_adenovirus_report.html")
        write_adenovirus_html_report(row, report_path, run_tag)
    print(f"Adenovirus HTML reports saved to {html_dir}")