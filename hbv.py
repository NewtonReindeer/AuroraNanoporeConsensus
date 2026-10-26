"""HBV reverse-transcriptase resistance and surface escape analysis."""
import os
import re

from .alignment import read_codon_from_aligned_segment
from .config import (
    HBV_MEDICATIONS,
    HBV_MIN_COVERAGE,
    HBV_MIN_REPORTABLE_COVERAGE,
    HBV_MIN_VARIANT_FREQUENCY,
    HBV_MUTATION_MEDICATION_MAP,
    HBV_MUTATION_SCRIPT_VERSION,
    HBV_RT_PANEL,
    HBV_S_ESCAPE_PANEL,
    HBV_S_ESCAPE_ASSOCIATIONS,
)
from .filesystem import app_base_dir
from .fasta import circular_codon, read_named_fasta_sequence, translate_codon

try:
    import pysam
except ImportError:
    pysam = None


def hbv_medication_summary(details):
    by_medication = {drug: [] for drug in HBV_MEDICATIONS}
    unmapped = []

    for detail in details or []:
        if detail.get("gene") != "RT":
            continue
        ref_aa = hbv_panel_reference_aa(detail)
        observed_ref_aa = str(detail.get("reference_aa", ""))
        aa_pos = detail.get("aa_position", "")
        label = detail.get("position", "")
        for call in detail.get("calls", []) or []:
            alt_aa = str(call.get("aa", ""))
            if not alt_aa or alt_aa in {ref_aa, observed_ref_aa}:
                continue
            mutation = f"{ref_aa}{aa_pos}{alt_aa}"
            frequency = call.get("frequency", 0)
            display = f"rt{mutation} ({frequency}%)"
            mappings = HBV_MUTATION_MEDICATION_MAP.get(mutation, {})
            if not mappings:
                unmapped.append(f"{label}{alt_aa} ({frequency}%)")
                continue
            for drug, evidence in mappings.items():
                by_medication.setdefault(drug, []).append(f"{display} [{evidence}]")

    result = {
        f"hbv_drug_{drug.lower()}": "; ".join(by_medication.get(drug, [])) or " "
        for drug in HBV_MEDICATIONS
    }
    result["hbv_drug_unmapped"] = "; ".join(unmapped) or " "
    return result

def parse_hbv_metadata_tsv(tsv_path):
    """Read HBV reference metadata from a TSV file.

    Supported header columns include ref_name/reference/name, rt_codon1_start,
    s_codon1_start, and an optional genotype column. Older three-column TSVs
    remain fully supported.
    """
    metadata = {}
    if not tsv_path or not os.path.exists(tsv_path):
        return metadata
    try:
        with open(tsv_path, "r", encoding="utf-8", errors="replace") as handle:
            header = None
            for raw in handle:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if header is None and any(
                    x.lower() in {"ref_name", "reference", "name", "rt_codon1_start", "s_codon1_start", "genotype"}
                    for x in parts
                ):
                    header = [x.strip().lower() for x in parts]
                    continue

                genotype = ""
                if header:
                    row = dict(zip(header, parts))
                    name = row.get("ref_name") or row.get("reference") or row.get("name")
                    rt_start = row.get("rt_codon1_start") or row.get("rt_start")
                    s_start = row.get("s_codon1_start") or row.get("s_start")
                    genotype = (row.get("genotype") or row.get("hbv_genotype") or "").strip().upper()
                elif len(parts) >= 3:
                    name, rt_start, s_start = parts[:3]
                    if len(parts) >= 4:
                        genotype = parts[3].strip().upper()
                else:
                    continue

                if not name:
                    continue
                try:
                    metadata[name] = {
                        "rt_codon1_start": int(rt_start),
                        "s_codon1_start": int(s_start),
                        "genotype": genotype,
                    }
                except (TypeError, ValueError):
                    continue
    except Exception as exc:
        print(f"Error reading HBV metadata TSV {tsv_path}: {exc}")
    return metadata

def parse_hbv_metadata_from_fasta_header(fasta_path, ref_name):
    if not fasta_path or not ref_name:
        return None
    try:
        with open(fasta_path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                if header.split()[0] != str(ref_name).split()[0]:
                    continue
                rt = re.search(r"(?:rt_codon1_start|rt_start)=(\d+)", header, re.I)
                s = re.search(r"(?:s_codon1_start|s_start)=(\d+)", header, re.I)
                genotype = re.search(r"(?:hbv_genotype|genotype)=([A-J])\b", header, re.I)
                if rt and s:
                    return {
                        "rt_codon1_start": int(rt.group(1)),
                        "s_codon1_start": int(s.group(1)),
                        "genotype": genotype.group(1).upper() if genotype else "",
                    }
    except Exception as exc:
        print(f"Error parsing HBV FASTA metadata: {exc}")
    return None

def get_hbv_reference_metadata(ref_path, ref_name):
    names = [str(ref_name or ""), str(ref_name or "").split()[0], str(ref_name or "").split('.')[0]]
    for candidate in [
        os.path.join(os.path.dirname(ref_path), "hbv_reference_metadata.tsv") if ref_path else None,
        os.path.join(app_base_dir(), "hbv_reference_metadata.tsv"),
    ]:
        data = parse_hbv_metadata_tsv(candidate)
        for name in names:
            if name in data:
                return data[name]
    return parse_hbv_metadata_from_fasta_header(ref_path, ref_name)

def call_hbv_aa_position_from_bam(bam_file, ref_name, codon_start_1based, reference_aa, reportable_aas):
    if pysam is None:
        raise ImportError("pysam is required for HBV mutation calling. Install with: pip install pysam")
    from collections import Counter
    counts = Counter()
    positions = [codon_start_1based, codon_start_1based + 1, codon_start_1based + 2]
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(ref_name, codon_start_1based - 1, codon_start_1based + 2):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail:
                continue
            codon = read_codon_from_aligned_segment(read, positions)
            if not codon:
                continue
            aa = translate_codon(codon)
            if aa != "X":
                counts[aa] += 1
    depth = sum(counts.values())
    calls = []
    for aa, count in counts.most_common():
        frequency = count / depth if depth else 0.0
        if aa != reference_aa and aa in reportable_aas and frequency >= HBV_MIN_VARIANT_FREQUENCY:
            calls.append({"aa": aa, "count": count, "frequency": round(frequency * 100, 1)})
    return {"depth": depth, "calls": calls, "aa_distribution": dict(counts)}

def format_hbv_position_call(result):
    calls = result.get("calls", []) if result else []
    if not calls:
        return " "
    return " | ".join(f"{call['aa']} : {call['frequency']}%" for call in calls)

def hbv_reference_genotype(metadata, ref_name):
    """Return genotype from metadata, with a conservative ref_<letter> fallback."""
    genotype = str((metadata or {}).get("genotype", "") or "").strip().upper()
    if re.fullmatch(r"[A-J]", genotype):
        return genotype
    match = re.search(r"(?:^|[_-])ref[_-]?([a-j])(?:$|[_-])", str(ref_name or ""), re.I)
    if not match:
        match = re.fullmatch(r"ref[_-]?([a-j])", str(ref_name or ""), re.I)
    return match.group(1).upper() if match else ""


def _surface_mutation_key(reference_aa, aa_pos, alt_aa):
    return f"s{reference_aa}{aa_pos}{alt_aa}"


def _surface_call_has_curated_escape(reference_aa, aa_pos, alt_aa):
    evidence = HBV_S_ESCAPE_ASSOCIATIONS.get(_surface_mutation_key(reference_aa, aa_pos, alt_aa), {})
    labels = set(evidence.get("associations", ()) or ())
    return bool(labels - {"Terapiassocierad"})


def determine_hbv_interpretation(position_results, low_coverage_positions):
    """High-level result; therapy-associated surface changes do not equal escape."""
    if low_coverage_positions:
        return "Null"

    relevant_frequencies = []
    has_rt = False
    has_escape = False
    has_surface_other = False
    for key, result in position_results.items():
        gene, aa_pos = key
        ref_aa = str(result.get("reference_aa", "") or "")
        for call in result.get("calls", []) or []:
            try:
                frequency = float(call.get("frequency", 0)) / 100.0
            except (TypeError, ValueError):
                continue
            alt_aa = str(call.get("aa", "") or "")
            if gene == "RT":
                has_rt = True
                relevant_frequencies.append(frequency)
            elif gene == "S":
                if _surface_call_has_curated_escape(ref_aa, aa_pos, alt_aa):
                    has_escape = True
                    relevant_frequencies.append(frequency)
                else:
                    has_surface_other = True

    if not relevant_frequencies:
        if has_surface_other:
            return "Ytantigenförändring utan litteraturbaserad escape-klassificering"
        return "Inga panelmutationer påvisade"

    population = "minoritetspopulation" if max(relevant_frequencies) < 0.50 else "majoritetspopulation"
    if has_rt and has_escape:
        return f"RT-resistens och ytantigen-escape associerad {population}"
    if has_rt:
        return f"RT-resistensassocierad {population}"
    return f"Ytantigen-escape associerad {population}"

def analyse_hbv_mutations(sample, bam_file, ref_name, ref_path):
    base = {
        "hbv_quality": "✕",
        "hbv_interpretation": "Ej bedömbar",
        "hbv_rt_quality": "✕",
        "hbv_rt_assessment_status": "NOT_EVALUABLE",
        "hbv_rt_interpretation": "Ej bedömbar",
        "hbv_rt_low_coverage_positions": " ",
        "hbv_rt_unassessable_positions": " ",
        "hbv_s_quality": "✕",
        "hbv_s_assessment_status": "NOT_EVALUABLE",
        "hbv_s_low_coverage_positions": " ",
        "hbv_s_unassessable_positions": " ",
        "hbv_low_coverage_positions": " ",
        "hbv_rt_priority_flags": " ",
        "hbv_rt_caution_flags": " ",
        "hbv_drug_lamivudine": " ",
        "hbv_drug_adefovir": " ",
        "hbv_drug_entecavir": " ",
        "hbv_drug_tenofovir": " ",
        "hbv_drug_unmapped": " ",
        "hbv_mutation_details": [],
        "hbv_mutation_script_version": HBV_MUTATION_SCRIPT_VERSION,
        "hbv_rt_codon1_start": "",
        "hbv_s_codon1_start": "",
        "hbv_coordinate_method": "",
        "hbv_reference_genotype": "",
    }
    for pos in HBV_RT_PANEL:
        base[f"hbv_rt_{pos}"] = " "
    for pos in HBV_S_ESCAPE_PANEL:
        base[f"hbv_s_{pos}"] = " "

    metadata = get_hbv_reference_metadata(ref_path, ref_name)
    if not metadata:
        base["hbv_low_coverage_positions"] = (
            "HBV-koordinatmetadata saknas: lägg till rt_codon1_start och "
            "s_codon1_start i FASTA-huvudet eller hbv_reference_metadata.tsv"
        )
        return base

    base["hbv_rt_codon1_start"] = metadata.get("rt_codon1_start", "")
    base["hbv_s_codon1_start"] = metadata.get("s_codon1_start", "")
    base["hbv_coordinate_method"] = "metadata"
    base["hbv_reference_genotype"] = hbv_reference_genotype(metadata, ref_name)

    ref_seq = read_named_fasta_sequence(ref_path, ref_name)
    if not ref_seq:
        base["hbv_low_coverage_positions"] = "Den valda HBV-referenssekvensen kunde inte läsas"
        return base

    low_cov = []
    rt_low_cov = []
    rt_unassessable = []
    s_low_cov = []
    s_unassessable = []
    priority_flags = []
    caution_flags = []
    position_results = {}
    details = []

    try:
        for gene, start_key, panel in [
            ("RT", "rt_codon1_start", HBV_RT_PANEL),
            ("S", "s_codon1_start", HBV_S_ESCAPE_PANEL),
        ]:
            gene_start = int(metadata[start_key])
            for aa_pos, rule in panel.items():
                codon_start = gene_start + ((aa_pos - 1) * 3)
                ref_codon = circular_codon(ref_seq, codon_start)
                ref_aa = translate_codon(ref_codon)
                result = call_hbv_aa_position_from_bam(
                    bam_file, ref_name, codon_start, ref_aa, rule["aas"]
                )
                key = (gene, aa_pos)
                result["reference_aa"] = ref_aa
                position_results[key] = result
                label = rule["label"]
                depth = int(result.get("depth", 0) or 0)
                if depth < HBV_MIN_COVERAGE:
                    low_cov.append(label)
                    if gene == "RT":
                        rt_low_cov.append(label)
                    else:
                        s_low_cov.append(label)
                if depth < HBV_MIN_REPORTABLE_COVERAGE:
                    if gene == "RT":
                        rt_unassessable.append(label)
                    else:
                        s_unassessable.append(label)
                    # Preserve raw read evidence for audit, but do not promote a
                    # mutation call below the minimum reportable codon depth.
                    result["raw_calls"] = list(result.get("calls", []) or [])
                    result["calls"] = []

                formatted = format_hbv_position_call(result)
                base[f"hbv_{gene.lower()}_{aa_pos}"] = formatted
                flag = rule.get("flag", "ESCAPE" if gene == "S" else "")
                if formatted.strip():
                    labelled = f"[{flag}] {label}: {formatted.strip()}" if flag else f"{label}: {formatted.strip()}"
                    if flag == "PRIORITY":
                        priority_flags.append(labelled)
                    elif flag == "CAUTION":
                        caution_flags.append(labelled)

                aa_dist = result.get("aa_distribution", {}) or {}
                dist_parts = []
                threshold_profile = []
                nonpanel_variants = []
                reportable_aas = set(rule.get("aas", set()) or set())
                for aa, count in sorted(aa_dist.items(), key=lambda kv: (-kv[1], kv[0])):
                    freq = round((count / depth) * 100, 1) if depth else 0
                    dist_parts.append(f"{aa}:{count} ({freq}%)")
                    if depth >= HBV_MIN_REPORTABLE_COVERAGE and freq >= HBV_MIN_VARIANT_FREQUENCY * 100:
                        threshold_profile.append({"aa": aa, "count": count, "frequency": freq})
                        if gene == "RT" and aa != ref_aa and aa not in reportable_aas:
                            nonpanel_variants.append({"aa": aa, "count": count, "frequency": freq})

                coverage_status = (
                    "preferred" if depth >= HBV_MIN_COVERAGE
                    else "limited" if depth >= HBV_MIN_REPORTABLE_COVERAGE
                    else "insufficient"
                )

                details.append({
                    "gene": gene,
                    "position": label,
                    "aa_position": aa_pos,
                    "reference_codon": ref_codon,
                    "reference_aa": ref_aa,
                    "panel_reference_aa": re.sub(r"^(?:rt|s)([A-Z*])\d+$", r"\1", label),
                    "reportable_aas": ",".join(sorted(rule.get("aas", set()))),
                    "codon_start": codon_start,
                    "codon_range": f"{codon_start}-{codon_start + 2}",
                    "depth": depth,
                    "aa_distribution": "; ".join(dist_parts),
                    "aa_profile_above_threshold": threshold_profile,
                    "observed_nonpanel_variants": nonpanel_variants,
                    "coverage_status": coverage_status,
                    "reported_calls": formatted.strip(),
                    "calls": list(result.get("calls", []) or []),
                    "raw_calls": list(result.get("raw_calls", result.get("calls", [])) or []),
                    "flag": flag,
                    "category": rule.get("category", "immune/vaccine/diagnostic escape-associated"),
                    "drugs": rule.get("drugs", ""),
                    "reference_genotype": base.get("hbv_reference_genotype", ""),
                })

        medication_summary = hbv_medication_summary(details)
        base.update(medication_summary)

        rt_depths = [int(d.get("depth", 0) or 0) for d in details if d.get("gene") == "RT"]
        s_depths = [int(d.get("depth", 0) or 0) for d in details if d.get("gene") == "S"]
        rt_status = (
            "FULL" if rt_depths and all(d >= HBV_MIN_COVERAGE for d in rt_depths)
            else "LIMITED" if any(d >= HBV_MIN_REPORTABLE_COVERAGE for d in rt_depths)
            else "NOT_EVALUABLE"
        )
        s_status = (
            "FULL" if s_depths and all(d >= HBV_MIN_COVERAGE for d in s_depths)
            else "LIMITED" if any(d >= HBV_MIN_REPORTABLE_COVERAGE for d in s_depths)
            else "NOT_EVALUABLE"
        )
        rt_position_results = {k: v for k, v in position_results.items() if k[0] == "RT"}
        rt_interpretation = (
            determine_hbv_interpretation(rt_position_results, [])
            if rt_status != "NOT_EVALUABLE"
            else "Null"
        )

        base.update({
            # Keep the legacy all-panel quality field unchanged for backward compatibility.
            "hbv_quality": "✕" if low_cov else "✓",
            "hbv_interpretation": determine_hbv_interpretation(position_results, low_cov),
            # New HBV-only scope fields: RT resistance is no longer invalidated by S/escape coverage.
            "hbv_rt_quality": "✓" if rt_status == "FULL" else ("LIMITED" if rt_status == "LIMITED" else "✕"),
            "hbv_rt_assessment_status": rt_status,
            "hbv_rt_interpretation": rt_interpretation,
            "hbv_rt_low_coverage_positions": ", ".join(rt_low_cov) if rt_low_cov else " ",
            "hbv_rt_unassessable_positions": ", ".join(rt_unassessable) if rt_unassessable else " ",
            "hbv_s_quality": "✓" if s_status == "FULL" else ("LIMITED" if s_status == "LIMITED" else "✕"),
            "hbv_s_assessment_status": s_status,
            "hbv_s_low_coverage_positions": ", ".join(s_low_cov) if s_low_cov else " ",
            "hbv_s_unassessable_positions": ", ".join(s_unassessable) if s_unassessable else " ",
            "hbv_low_coverage_positions": ", ".join(low_cov) if low_cov else " ",
            "hbv_rt_priority_flags": "; ".join(priority_flags) if priority_flags else " ",
            "hbv_rt_caution_flags": "; ".join(caution_flags) if caution_flags else " ",
            "hbv_mutation_details": details,
        })
    except Exception as exc:
        base["hbv_quality"] = "✕"
        base["hbv_interpretation"] = "Ej bedömbar"
        base["hbv_low_coverage_positions"] = f"HBV-mutationsanalysen misslyckades: {exc}"
    return base

def hbv_swedish_evidence_label(evidence):
    labels = {
        "established": "etablerad",
        "associated": "associerad",
        "associated/background": "associerad/bakgrund",
        "resistance background": "resistensbakgrund",
        "context-dependent": "kontextberoende",
        "candidate association - caution": "kandidatassociation – försiktighet",
        "reported association - caution": "rapporterad association – försiktighet",
    }
    return labels.get(str(evidence or ""), str(evidence or ""))

def hbv_translate_free_text(value):
    text = str(value or "")
    replacements = {
        "immune/vaccine/diagnostic escape-associated": "immun-, vaccin- eller diagnostikescape-associerad",
        "lamivudine": "lamivudin",
        "adefovir": "adefovir",
        "entecavir": "entekavir",
        "tenofovir": "tenofovir",
        "compensatory/background": "kompensatorisk/bakgrund",
        "primary resistance": "primär resistens",
        "resistance background": "resistensbakgrund",
        "context-dependent": "kontextberoende",
        "associated": "associerad",
    }
    for old, new in replacements.items():
        text = re.sub(re.escape(old), new, text, flags=re.I)
    return text

def hbv_population_label(frequency):
    try:
        return "Minoritetspopulation" if float(frequency) < 50.0 else "Majoritetspopulation"
    except (TypeError, ValueError):
        return "Okänd population"

def hbv_panel_reference_aa(detail):
    """Return the canonical wild-type amino acid encoded in labels such as rtV173."""
    label = str(detail.get("position", "") or "")
    match = re.match(r"^(?:rt|s)([A-Z*])\d+$", label)
    return match.group(1) if match else str(detail.get("reference_aa", "") or "")

def hbv_detected_calls(row):
    calls = []
    for detail in row.get("hbv_mutation_details", []) or []:
        if int(detail.get("depth", 0) or 0) < HBV_MIN_REPORTABLE_COVERAGE:
            continue
        canonical_ref_aa = hbv_panel_reference_aa(detail)
        observed_ref_aa = str(detail.get("reference_aa", "") or "")
        for call in detail.get("calls", []) or []:
            alt_aa = str(call.get("aa", "") or "")
            # Mutation naming is genotype-aware for S: use the amino acid in the
            # selected genotype-specific reference, not the canonical panel label.
            if not alt_aa:
                continue
            gene = detail.get("gene", "")
            if gene == "S":
                if alt_aa == observed_ref_aa:
                    continue
                mutation = _surface_mutation_key(observed_ref_aa, detail.get("aa_position", ""), alt_aa)
            else:
                if alt_aa in {canonical_ref_aa, observed_ref_aa}:
                    continue
                mutation = f"{detail.get('position', '')}{alt_aa}"
            medication_map = {}
            if gene == "RT":
                key = f"{canonical_ref_aa}{detail.get('aa_position', '')}{alt_aa}"
                medication_map = HBV_MUTATION_MEDICATION_MAP.get(key, {})
            calls.append({
                "gene": gene,
                "mutation": mutation,
                "frequency": call.get("frequency", 0),
                "count": call.get("count", 0),
                "depth": detail.get("depth", 0),
                "population": hbv_population_label(call.get("frequency", 0)),
                "flag": detail.get("flag", ""),
                "category": detail.get("category", ""),
                "medications": medication_map,
                "reference_aa": observed_ref_aa,
                "reference_genotype": detail.get("reference_genotype", row.get("hbv_reference_genotype", "")),
                "coverage_status": detail.get("coverage_status", "preferred"),
            })
    return calls

def hbv_swedish_interpretation(value):
    translations = {
        "No panel mutation detected": "Inga panelmutationer påvisade",
        "Null": "Ej bedömbar",
        "RT resistance-associated minor population": "RT-resistensassocierad minoritetspopulation",
        "RT resistance-associated majority population": "RT-resistensassocierad majoritetspopulation",
        "Surface escape-associated minor population": "Ytantigen-escape associerad minoritetspopulation",
        "Surface escape-associated majority population": "Ytantigen-escape associerad majoritetspopulation",
        "RT resistance + surface escape-associated minor population": "RT-resistens och ytantigen-escape associerad minoritetspopulation",
        "RT resistance + surface escape-associated majority population": "RT-resistens och ytantigen-escape associerad majoritetspopulation",
    }
    return translations.get(str(value or ""), str(value or ""))

def hbv_rt_call_index(row):
    """Return detected RT calls keyed by canonical mutation name (for example rtM204V)."""
    return {
        str(item.get("mutation", "")): item
        for item in hbv_detected_calls(row)
        if item.get("gene") == "RT" and item.get("mutation")
    }

def hbv_pattern_mutation_text(index, mutations):
    parts = []
    for mutation in mutations:
        item = index.get(mutation)
        if item:
            parts.append(f"{mutation} ({float(item.get('frequency', 0) or 0):.1f} %)")
    return "; ".join(parts)

def hbv_pattern_max_frequency(pattern):
    """Return the highest reported component frequency in a combination result."""
    values = []
    for match in re.finditer(r"\(([0-9]+(?:\.[0-9]+)?)\s*%\)", str(pattern.get("mutations", "") or "")):
        try:
            values.append(float(match.group(1)))
        except ValueError:
            pass
    return max(values) if values else 0.0

def hbv_apply_population_context(patterns):
    """Separate biological rule strength from population dominance.

    An established resistance-associated pattern below 50% remains a real finding,
    but is presented as a yellow minority-population warning rather than as a
    dominant resistant population.
    """
    adjusted = []
    for original in patterns:
        pattern = dict(original)
        maximum = hbv_pattern_max_frequency(pattern)
        pattern["max_frequency"] = maximum
        if pattern.get("level") == "high":
            if maximum and maximum < 50.0:
                original_title = pattern.get("title", "Resistensassocierat mönster")
                pattern["level"] = "minor"
                pattern["title"] = "Resistensassocierat mönster i minoritetspopulation"
                pattern["interpretation"] = (
                    f"{original_title} påvisades, men den högsta rapporterade variantandelen "
                    f"är {maximum:.1f} %. Ingen dominerande resistensassocierad population "
                    "påvisades. Fyndet är en varningssignal och bör bedömas tillsammans med "
                    "behandlingshistorik och HBV-DNA-utveckling. " + pattern.get("interpretation", "")
                )
            else:
                pattern["level"] = "major"
                pattern["title"] = "Dominerande resistensassocierat mönster"
                pattern["interpretation"] = (
                    f"Ett etablerat resistensassocierat mönster påvisades i majoritetspopulation "
                    f"(högsta variantandel {maximum:.1f} %). Fyndet ger starkt genotypiskt stöd "
                    "för nedsatt läkemedelskänslighet. " + pattern.get("interpretation", "")
                )
        adjusted.append(pattern)
    return adjusted

def hbv_combination_interpretations(row):
    """Interpret established HBV RT mutation patterns without changing variant calls.

    Calls are not read-phased across distant codons. A detected combination therefore
    means that the component substitutions coexist in the sample, not necessarily on
    the same viral genome. Low-coverage samples are not assigned a definitive pattern.
    """
    rt_status = str(row.get("hbv_rt_assessment_status", "") or "").upper()
    if not rt_status:
        rt_status = "FULL" if str(row.get("hbv_quality", "")) == "✓" else "NOT_EVALUABLE"
    if rt_status == "NOT_EVALUABLE":
        return [{
            "category": "Teknisk begränsning",
            "level": "unknown",
            "title": "Kombinationsbedömning ej möjlig",
            "interpretation": (
                "Ingen RT-panelposition har tillräckligt lokalt kodondjup för rapportering. "
                "Resistensmönster ska därför inte bedömas."
            ),
            "mutations": "",
        }]

    idx = hbv_rt_call_index(row)
    present = set(idx)
    results = []

    m204 = [m for m in ("rtM204I", "rtM204V") if m in present]
    lam_background = bool(m204)
    l180 = "rtL180M" in present
    v173 = "rtV173L" in present
    a181 = [m for m in ("rtA181T", "rtA181V") if m in present]
    n236 = "rtN236T" in present
    etv_secondary = [m for m in present if re.match(r"^rt(?:T184|S202|M250)[A-Z*]$", m)]
    i169 = "rtI169T" in present

    # Lamivudine/telbivudine pathway.
    if lam_background and l180:
        muts = ["rtL180M"] + m204
        results.append({
            "category": "Lamivudin/telbivudin",
            "level": "high",
            "title": "Etablerat resistensmönster",
            "interpretation": (
                "Kombinationen rtL180M + rtM204I/V är ett klassiskt och väl etablerat "
                "resistensmönster mot lamivudin och telbivudin."
            ),
            "mutations": hbv_pattern_mutation_text(idx, muts),
        })
    elif lam_background:
        results.append({
            "category": "Lamivudin/telbivudin",
            "level": "high",
            "title": "Primär resistensmutation påvisad",
            "interpretation": (
                "rtM204I/V är en primär resistensmutation mot lamivudin/telbivudin. "
                "Avsaknad av rtL180M utesluter inte resistens."
            ),
            "mutations": hbv_pattern_mutation_text(idx, m204),
        })

    if v173 and l180 and "rtM204V" in present:
        muts = ["rtV173L", "rtL180M", "rtM204V"]
        results.append({
            "category": "Lamivudin/telbivudin",
            "level": "high",
            "title": "Kompenserat resistensmönster",
            "interpretation": (
                "rtV173L + rtL180M + rtM204V är ett etablerat lamivudinresistent mönster "
                "där rtV173L kan återställa replikationsförmåga."
            ),
            "mutations": hbv_pattern_mutation_text(idx, muts),
        })

    # Entecavir requires a lamivudine-resistant background plus secondary changes.
    if lam_background and etv_secondary:
        muts = (["rtL180M"] if l180 else []) + m204 + sorted(etv_secondary)
        results.append({
            "category": "Entekavir",
            "level": "high",
            "title": "Etablerat entekavirresistensmönster",
            "interpretation": (
                "Lamivudinresistent bakgrund med rtM204I/V och ytterligare substitution "
                "vid rtT184, rtS202 eller rtM250 är ett etablerat entekavirresistensmönster."
            ),
            "mutations": hbv_pattern_mutation_text(idx, muts),
        })
    elif etv_secondary:
        results.append({
            "category": "Entekavir",
            "level": "caution",
            "title": "Ofullständigt entekavirmönster",
            "interpretation": (
                "En entekavirassocierad sekundär mutation påvisades utan rtM204I/V. "
                "Fyndet ska inte ensamt tolkas som etablerad entekavirresistens."
            ),
            "mutations": hbv_pattern_mutation_text(idx, sorted(etv_secondary)),
        })

    if lam_background and i169:
        muts = (["rtL180M"] if l180 else []) + m204 + ["rtI169T"]
        results.append({
            "category": "Entekavir",
            "level": "caution",
            "title": "Entekavirassocierat kombinationsfynd",
            "interpretation": (
                "rtI169T tillsammans med lamivudinresistent bakgrund är associerat med "
                "minskad entekavirkänslighet, men bör bedömas tillsammans med övriga RT-fynd."
            ),
            "mutations": hbv_pattern_mutation_text(idx, muts),
        })
    elif i169:
        results.append({
            "category": "Entekavir",
            "level": "caution",
            "title": "Kontextberoende entekavirfynd",
            "interpretation": "rtI169T utan lamivudinresistent bakgrund är inte tillräckligt för etablerad entekavirresistens.",
            "mutations": hbv_pattern_mutation_text(idx, ["rtI169T"]),
        })

    # Adefovir pathway.
    adv_muts = a181 + (["rtN236T"] if n236 else [])
    if a181 and n236:
        results.append({
            "category": "Adefovir",
            "level": "high",
            "title": "Kombinerat adefovirresistensmönster",
            "interpretation": (
                "Substitution vid rtA181 tillsammans med rtN236T utgör ett starkt "
                "adefovirresistensmönster."
            ),
            "mutations": hbv_pattern_mutation_text(idx, adv_muts),
        })
    elif adv_muts:
        results.append({
            "category": "Adefovir",
            "level": "high",
            "title": "Etablerad adefovirresistensmutation",
            "interpretation": "rtA181T/V eller rtN236T är etablerat associerad med adefovirresistens.",
            "mutations": hbv_pattern_mutation_text(idx, adv_muts),
        })

    # Cross-/multidrug pattern.
    if lam_background and adv_muts:
        muts = (["rtL180M"] if l180 else []) + m204 + adv_muts
        results.append({
            "category": "Multiläkemedelsresistens",
            "level": "high",
            "title": "Mönster förenligt med resistens i flera läkemedelsklasser",
            "interpretation": (
                "Lamivudinresistent bakgrund tillsammans med rtA181T/V och/eller rtN236T "
                "är förenlig med kors- eller multiläkemedelsresistens."
            ),
            "mutations": hbv_pattern_mutation_text(idx, muts),
        })

    # Tenofovir remains cautionary.
    if adv_muts:
        results.append({
            "category": "Tenofovir",
            "level": "caution",
            "title": "Tenofovir – försiktighetsfynd",
            "interpretation": (
                "rtA181T/V och/eller rtN236T kan förekomma i behandlade populationer med "
                "minskad känslighet, men utgör inte ensamt bevis för klinisk tenofovirresistens."
            ),
            "mutations": hbv_pattern_mutation_text(idx, adv_muts),
        })

    if not results:
        results.append({
            "category": "Samlad RT-bedömning",
            "level": "none",
            "title": "Inget etablerat kombinationsmönster påvisat",
            "interpretation": (
                "Inget av de regelbaserade resistensmönstren påvisades över rapporteringströskeln "
                + ("vid de bedömbara RT-positionerna. " if rt_status == "LIMITED" else ". ")
                + "Enskilda associationer redovisas separat i mutationstabellen."
            ),
            "mutations": "",
        })
    if rt_status == "LIMITED":
        low = str(row.get("hbv_rt_low_coverage_positions", "") or "").strip()
        unassessable = str(row.get("hbv_rt_unassessable_positions", "") or "").strip()
        limitation = "RT-bedömningen är begränsad av lokalt kodondjup under den föredragna nivån på 100 läsningar."
        if unassessable:
            limitation += f" Ej bedömbara RT-positioner (<{HBV_MIN_REPORTABLE_COVERAGE} läsningar): {unassessable}."
        elif low:
            limitation += f" Positioner med 30–99 läsningar: {low}."
        results.append({
            "category": "Teknisk begränsning",
            "level": "caution",
            "title": "Begränsad lokal RT-täckning",
            "interpretation": limitation,
            "mutations": "",
        })
    return hbv_apply_population_context(results)

def hbv_combination_summary_text(row):
    patterns = hbv_combination_interpretations(row)
    important = [p for p in patterns if p.get("level") == "major"]
    selected = important or [p for p in patterns if p.get("level") == "minor"] or [p for p in patterns if p.get("level") == "caution"] or patterns
    return " | ".join(f"{p['category']}: {p['title']}" for p in selected)

def hbv_status_for_medication(row, medication):
    calls = []
    for item in hbv_detected_calls(row):
        evidence = item.get("medications", {}).get(medication)
        if evidence:
            calls.append((item, evidence))
    if not calls:
        return "Inga associerade mutationer påvisade", "", "none"
    max_frequency = max(float(item.get("frequency", 0) or 0) for item, _ in calls)
    level = "major" if max_frequency >= 50 else "minor"
    summary = "; ".join(
        f"{item['mutation']} ({item['frequency']} %, {hbv_swedish_evidence_label(evidence)})"
        for item, evidence in calls
    )
    if medication == "Tenofovir":
        status = "Försiktighet – litteraturassociation, ej ensam grund för resistensbedömning"
        level = "caution"
    elif level == "major":
        status = "Associerad mutation i majoritetspopulation"
    else:
        status = "Associerad mutation i minoritetspopulation"
    return status, summary, level
