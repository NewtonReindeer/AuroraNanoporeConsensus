"""HCV amplicon QC and NS5A resistance analysis."""
import gzip
import os
import re
import statistics
import subprocess

from .alignment import read_codon_from_aligned_segment
from .config import (
    HCV_MIN_VARIANT_FREQUENCY,
    HCV_MINOR_RAS_MAX_FREQUENCY,
    HCV_NS5A_PANEL_BY_SUBTYPE,
    HCV_NS5A_REFERENCE_MOTIF_BY_SUBTYPE,
    HCV_REFERENCE_METADATA,
    HCV_SCRIPT_VERSION,
)
# Minimum complete-codon read depth for HCV NS5A assessment.
# Keep the variant-frequency threshold from config.py unchanged.
HCV_NS5A_MIN_CODON_DEPTH = 80

from .filesystem import app_base_dir
from .fasta import circular_codon, read_named_fasta_sequence, read_single_fasta_sequence, translate_codon

try:
    import pysam
except ImportError:
    pysam = None


def calculate_hcv_amplicon_breadth_coverage_percent(bam_file, ref_name, start_1based, end_1based):
    """Calculate HCV breadth only across the detected fixed-primer amplicon."""
    if not bam_file or not os.path.exists(bam_file) or not ref_name:
        return 0.0
    if not start_1based or not end_1based or end_1based < start_1based:
        return 0.0

    target_length = end_1based - start_1based + 1
    covered_positions = 0
    region_string = f"{ref_name}:{start_1based}-{end_1based}"

    try:
        cmd = ["samtools", "depth", "-a", "-r", region_string, bam_file]
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    if int(fields[2]) > 0:
                        covered_positions += 1
                except ValueError:
                    continue
            _, stderr = proc.communicate()
            if proc.returncode != 0:
                print(f"samtools depth warning/error for HCV target {bam_file}: {stderr}")
                return 0.0
    except Exception as exc:
        print(f"Error calculating HCV target breadth for {bam_file}: {exc}")
        return 0.0

    return round((covered_positions / target_length) * 100, 2) if target_length else 0.0

def hcv_amplicon_consensus_stats(fasta_path, start_1based, end_1based):
    """Calculate consensus QC statistics for the detected HCV amplicon only."""
    seq = read_single_fasta_sequence(fasta_path)
    if not seq or not start_1based or not end_1based or end_1based < start_1based:
        return {
            "consensus_length": 0,
            "n_count": 0,
            "n_percent": 100.0,
            "non_n_length": 0,
        }

    seq = re.sub(r"[^ACGTN]", "", seq.upper())
    target_length = end_1based - start_1based + 1

    # samtools/Medaka commonly retain reference coordinates and fill uncovered
    # sequence with N. If a tool emits only the amplicon, use that sequence as-is.
    if len(seq) >= end_1based:
        target_seq = seq[start_1based - 1:end_1based]
    elif len(seq) <= target_length + 20:
        target_seq = seq
    else:
        target_seq = seq[:target_length]

    total = len(target_seq)
    n_count = target_seq.count("N")
    non_n = total - n_count
    return {
        "consensus_length": non_n,
        "n_count": n_count,
        "n_percent": round((n_count / total) * 100, 2) if total else 100.0,
        "non_n_length": non_n,
    }

def hcv_subtype_from_name(name):
    text = str(name or "").lower().replace("_", " ").replace("-", " ")
    compact = text.replace(" ", "")
    subtype_patterns = [
        ("1A", ["1a", "gt1a", "genotype1a", "subtype=1a"]),
        ("1B", ["1b", "gt1b", "genotype1b", "subtype=1b"]),
        ("2A", ["2a", "gt2a", "genotype2a", "subtype=2a"]),
        ("3A", ["3a", "gt3a", "genotype3a", "subtype=3a"]),
        ("3B", ["3b", "gt3b", "genotype3b", "subtype=3b"]),
        ("4A", ["4a", "gt4a", "genotype4a", "subtype=4a"]),
        ("4B", ["4b", "gt4b", "genotype4b", "subtype=4b"]),
    ]
    for subtype, patterns in subtype_patterns:
        if any(pattern in compact for pattern in patterns):
            return subtype
    return "Unknown"

def read_fastq_stats(fastq_gz_path):
    read_count = 0
    total_bases = 0
    try:
        with gzip.open(fastq_gz_path, "rt", encoding="utf-8", errors="replace") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    seq_len = len(line.strip())
                    read_count += 1
                    total_bases += seq_len
    except Exception as e:
        print(f"Error reading FASTQ stats in {fastq_gz_path}: {e}")

    mean_len = round(total_bases / read_count, 1) if read_count else 0
    return read_count, mean_len

def calculate_median_nonzero_depth(bam_file, ref_name, region=None):
    depths = []
    if not bam_file or not os.path.exists(bam_file) or not ref_name:
        return 0

    try:
        region_string = (
            f"{ref_name}:{region[0]}-{region[1]}"
            if region and region[0] and region[1]
            else ref_name
        )
        cmd = ["samtools", "depth", "-r", region_string, bam_file]
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    depth = int(fields[2])
                    if depth > 0:
                        depths.append(depth)
                except ValueError:
                    continue
            _, stderr = proc.communicate()
            if proc.returncode != 0:
                print(f"samtools depth warning/error for {bam_file}: {stderr}")
    except Exception as e:
        print(f"Error calculating non-zero median depth for {bam_file}: {e}")

    return round(statistics.median(depths), 1) if depths else 0

def parse_hcv_metadata_tsv(tsv_path):
    metadata = {}
    if not tsv_path or not os.path.exists(tsv_path):
        return metadata

    try:
        with open(tsv_path, "r", encoding="utf-8", errors="replace") as f:
            header = None
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if header is None and any(x.lower() in ["ref_name", "reference", "subtype"] for x in parts):
                    header = [x.strip().lower() for x in parts]
                    continue

                if header:
                    row = dict(zip(header, parts))
                    ref_name = row.get("ref_name") or row.get("reference") or row.get("name")
                    subtype = row.get("subtype") or row.get("genotype") or "Unknown"
                    start = row.get("ns5a_codon1_start") or row.get("ns5a_start") or row.get("codon1_start")
                else:
                    if len(parts) < 3:
                        continue
                    ref_name, subtype, start = parts[0], parts[1], parts[2]

                if not ref_name or not start:
                    continue
                try:
                    metadata[ref_name] = {
                        "subtype": str(subtype).upper(),
                        "ns5a_codon1_start": int(start),
                    }
                except ValueError:
                    continue
    except Exception as e:
        print(f"Error reading HCV metadata TSV {tsv_path}: {e}")

    return metadata

def parse_hcv_metadata_from_fasta_header(fasta_path, ref_name):
    if not fasta_path or not os.path.exists(fasta_path) or not ref_name:
        return None

    try:
        with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                header = line[1:].strip()
                name_only = header.split()[0]
                if name_only != ref_name:
                    continue

                subtype_match = re.search(r"(?:subtype|genotype|gt)=([A-Za-z0-9]+)", header, re.IGNORECASE)
                start_match = re.search(r"(?:ns5a_codon1_start|ns5a_start|codon1_start)=(\d+)", header, re.IGNORECASE)
                if start_match:
                    subtype = subtype_match.group(1).upper() if subtype_match else hcv_subtype_from_name(header)
                    return {
                        "subtype": subtype,
                        "ns5a_codon1_start": int(start_match.group(1)),
                    }
    except Exception as e:
        print(f"Error parsing HCV FASTA header metadata from {fasta_path}: {e}")

    return None

def get_hcv_reference_metadata(ref_path, ref_name):
    """Return HCV reference metadata for the selected SAM reference name.

    Lookup is deliberately tolerant: exact accession, accession without version,
    local TSV next to the reference folder, FASTA header tokens, then built-ins.
    """
    names_to_try = []
    if ref_name:
        names_to_try.append(ref_name)
        names_to_try.append(str(ref_name).split()[0])
        names_to_try.append(str(ref_name).split('.')[0])

    candidates = []
    if ref_path:
        candidates.append(os.path.join(os.path.dirname(ref_path), "hcv_reference_metadata.tsv"))
    candidates.append(os.path.join(app_base_dir(), "hcv_reference_metadata.tsv"))

    for candidate in candidates:
        metadata = parse_hcv_metadata_tsv(candidate)
        for name in names_to_try:
            if name in metadata:
                return metadata[name]

    from_header = parse_hcv_metadata_from_fasta_header(ref_path, ref_name)
    if from_header:
        return from_header

    for name in names_to_try:
        if name in HCV_REFERENCE_METADATA:
            return HCV_REFERENCE_METADATA[name]

    return None

def call_hcv_ns5a_position_from_bam(
    bam_file,
    ref_name,
    codon_start_1based,
    reference_aa,
    resistance_aas,
    min_variant_frequency=HCV_MIN_VARIANT_FREQUENCY
):
    if pysam is None:
        raise ImportError("pysam is required for HCV minor-variant calling. Install with: pip install pysam")

    from collections import Counter

    codon_positions = [codon_start_1based, codon_start_1based + 1, codon_start_1based + 2]
    aa_counts = Counter()

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(ref_name, codon_start_1based - 1, codon_start_1based + 2):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate or read.is_qcfail:
                continue

            codon = read_codon_from_aligned_segment(read, codon_positions)
            if not codon:
                continue

            aa = translate_codon(codon)
            if aa != "X":
                aa_counts[aa] += 1

    depth = sum(aa_counts.values())
    calls = []

    for aa, count in aa_counts.most_common():
        freq = count / depth if depth else 0
        if aa != reference_aa and aa in resistance_aas and freq >= min_variant_frequency:
            calls.append({
                "aa": aa,
                "count": count,
                "frequency": round(freq * 100, 1),
            })

    return {
        "depth": depth,
        "calls": calls,
        "aa_distribution": dict(aa_counts),
    }

def format_hcv_position_call(call_result):
    calls = call_result.get("calls", []) if call_result else []
    if not calls:
        return " "
    return " | ".join(f"{call['aa']} : {call['frequency']}%" for call in calls)

def determine_hcv_interpretation(position_results, low_coverage_positions):
    if low_coverage_positions:
        return "Ej bedömbar"

    frequencies = []
    for result in position_results.values():
        for call in result.get("calls", []) or []:
            try:
                frequencies.append(float(call.get("frequency", 0)) / 100.0)
            except (TypeError, ValueError):
                continue

    if not frequencies:
        return "Känslig"
    if max(frequencies) < HCV_MINOR_RAS_MAX_FREQUENCY:
        return "RAS minoritetsvariant"
    return "Resistent"

def hcv_candidate_ns5a_starts(meta):
    starts = []
    try:
        if meta and meta.get("ns5a_codon1_start"):
            starts.append(int(meta.get("ns5a_codon1_start")))
    except Exception:
        pass

    # Common coordinate systems in this project:
    # 6258 = complete-genome-style HCV references with 5'UTR included.
    # 5917 = CDS/polyprotein-style references starting at the polyprotein ATG.
    for candidate in [6258, 5917]:
        if candidate not in starts:
            starts.append(candidate)
    return starts

def call_hcv_panel_for_start(
    bam_file,
    ref_name,
    ref_path,
    subtype,
    ns5a_codon1_start
):
    position_results = {}
    low_coverage = []
    details = []

    subtype = str(subtype or "Unknown").upper().strip()
    subtype_panel = HCV_NS5A_PANEL_BY_SUBTYPE.get(subtype, {})
    if not subtype_panel:
        raise ValueError(f"Subtype {subtype} does not have a validated NS5A resistance panel")

    reference_sequence = read_named_fasta_sequence(ref_path, ref_name)
    if not reference_sequence:
        raise ValueError(f"Could not read selected HCV reference sequence '{ref_name}' from {ref_path}")

    for aa_pos, rule in subtype_panel.items():
        codon_start = int(ns5a_codon1_start) + ((aa_pos - 1) * 3)
        reference_codon = circular_codon(reference_sequence, codon_start)
        reference_aa = translate_codon(reference_codon)
        label = rule.get("label") or (f"{reference_aa}{aa_pos}" if reference_aa != "X" else f"AA{aa_pos}")

        result = call_hcv_ns5a_position_from_bam(
            bam_file=bam_file,
            ref_name=ref_name,
            codon_start_1based=codon_start,
            reference_aa=reference_aa,
            resistance_aas=rule.get("reportable_aas", set()),
            min_variant_frequency=HCV_MIN_VARIANT_FREQUENCY,
        )

        # Keep numeric keys internally so output columns remain stable even when
        # the subtype reference residue is A28/D58 rather than M28/H58.
        position_results[aa_pos] = result
        depth = int(result.get("depth", 0) or 0)
        if depth < HCV_NS5A_MIN_CODON_DEPTH:
            low_coverage.append(label)

        aa_dist = result.get("aa_distribution", {}) or {}
        dist_parts = []
        for aa, count in sorted(aa_dist.items(), key=lambda kv: (-kv[1], kv[0])):
            freq = round((count / depth) * 100, 1) if depth else 0
            dist_parts.append(f"{aa}:{count} ({freq}%)")

        details.append({
            "position": label,
            "aa_position": aa_pos,
            "reference_codon": reference_codon,
            "reference_aa": reference_aa,
            "reportable_aas": ",".join(sorted(rule.get("reportable_aas", set()))),
            "codon_start": codon_start,
            "codon_range": f"{codon_start}-{codon_start + 2}",
            "depth": depth,
            "depth_threshold": HCV_NS5A_MIN_CODON_DEPTH,
            "depth_flag": ("Below 80x" if depth < HCV_NS5A_MIN_CODON_DEPTH else "80-99x: lower depth" if depth < 100 else "100x+"),
            "aa_distribution": "; ".join(dist_parts),
            "reported_calls": format_hcv_position_call(result).strip(),
        })

    positions_passing = sum(1 for d in details if d["depth"] >= HCV_NS5A_MIN_CODON_DEPTH)
    total_depth = sum(d["depth"] for d in details)
    min_depth = min([d["depth"] for d in details], default=0)

    return {
        "ns5a_start": int(ns5a_codon1_start),
        "position_results": position_results,
        "low_coverage": low_coverage,
        "details": details,
        "positions_passing": positions_passing,
        "total_depth": total_depth,
        "min_depth": min_depth,
    }

def score_hcv_ns5a_start(reference_sequence, subtype, candidate_start):
    motif = HCV_NS5A_REFERENCE_MOTIF_BY_SUBTYPE.get(str(subtype or "").upper(), {})
    residues = {}
    score = 0
    for aa_pos, expected_aa in motif.items():
        codon_start = int(candidate_start) + ((aa_pos - 1) * 3)
        aa = translate_codon(circular_codon(reference_sequence, codon_start))
        residues[aa_pos] = aa
        if aa == expected_aa:
            score += 1

    # Y93 is essential for a valid NS5A resistance coordinate system.
    y93_ok = residues.get(93) == "Y"
    return score, y93_ok, residues

def locate_hcv_ns5a_start(ref_path, ref_name, subtype, metadata_start):
    reference_sequence = read_named_fasta_sequence(ref_path, ref_name)
    if not reference_sequence:
        raise ValueError(f"Could not read selected HCV reference sequence '{ref_name}' from {ref_path}")

    try:
        metadata_start = int(metadata_start)
    except (TypeError, ValueError):
        metadata_start = 0
    if metadata_start <= 0:
        raise ValueError("Missing or invalid NS5A start coordinate")

    motif = HCV_NS5A_REFERENCE_MOTIF_BY_SUBTYPE.get(str(subtype or "").upper(), {})
    if not motif:
        return metadata_start, {"score": 0, "residues": {}, "offset": 0, "method": "metadata-only"}

    candidates = []
    # Search around the supplied coordinate. This resolves accession-specific
    # offsets and frame errors while retaining the expected genomic neighbourhood.
    for candidate in range(max(1, metadata_start - (750 if str(subtype).upper() == "1B" else 120)), metadata_start + (751 if str(subtype).upper() == "1B" else 121)):
        score, y93_ok, residues = score_hcv_ns5a_start(reference_sequence, subtype, candidate)
        candidates.append((y93_ok, score, -abs(candidate - metadata_start), candidate, residues))

    candidates.sort(reverse=True)
    y93_ok, score, _, selected_start, residues = candidates[0]
    minimum_score = max(3, len(motif) - 1)
    if not y93_ok or score < minimum_score:
        observed = ", ".join(f"{p}:{aa}" for p, aa in sorted(residues.items()))
        raise ValueError(
            f"Could not validate NS5A frame for subtype {subtype} near {metadata_start}; "
            f"best start {selected_start} matched {score}/{len(motif)} residues ({observed})"
        )

    return selected_start, {
        "score": score,
        "residues": residues,
        "offset": selected_start - metadata_start,
        "method": "reference-motif-autoframe",
    }

def hcv_panel_labels_for_subtype(subtype):
    """Return the clinically relevant report labels for the selected subtype panel."""
    subtype = str(subtype or "Unknown").upper().strip()
    subtype_panel = HCV_NS5A_PANEL_BY_SUBTYPE.get(subtype, {})
    return [rule.get("label", f"AA{aa_pos}") for aa_pos, rule in subtype_panel.items()]

def choose_hcv_panel_result(bam_file, ref_name, ref_path, meta):
    try:
        metadata_start = int(meta.get("ns5a_codon1_start", 0) or 0)
    except Exception:
        metadata_start = 0

    if metadata_start <= 0:
        return None

    subtype = str(meta.get("subtype", "Unknown")).upper().strip()

    try:
        selected_start, validation = locate_hcv_ns5a_start(
            ref_path, ref_name, subtype, metadata_start
        )
        result = call_hcv_panel_for_start(
            bam_file, ref_name, ref_path, subtype, selected_start
        )
        result["metadata_ns5a_start"] = metadata_start
        result["coordinate_validation"] = validation
        for detail in result.get("details", []):
            detail["metadata_ns5a_start"] = metadata_start
            detail["selected_offset"] = validation.get("offset", 0)
            detail["coordinate_method"] = validation.get("method", "")
            detail["motif_score"] = validation.get("score", "")
        return result
    except Exception as validation_error:
        # Never call resistance codons at an unvalidated metadata coordinate.
        # An incorrect NS5A start can produce zero depth even when amplicon QC
        # is HIGH, or worse, produce misleading amino-acid calls.
        failed_positions = hcv_panel_labels_for_subtype(subtype)
        error_text = f"NS5A coordinate validation failed: {validation_error}"
        return {
            "ns5a_start": "",
            "metadata_ns5a_start": metadata_start,
            "position_results": {},
            "low_coverage": failed_positions,
            "details": [{
                "position": "panel", "aa_position": "",
                "reference_codon": "", "reference_aa": "",
                "reportable_aas": "", "metadata_ns5a_start": metadata_start,
                "selected_offset": "", "coordinate_method": "validation-failed",
                "motif_score": "", "codon_range": "", "depth": "",
                "aa_distribution": "", "reported_calls": error_text,
            }],
            "positions_passing": 0, "total_depth": 0,
            "min_depth": 0, "error": error_text,
        }


def analyse_hcv_resistance(
    sample,
    input_fastq,
    bam_file,
    ref_name,
    ref_path,
    mapped_reads_1a,
    mapped_reads_1b,
    median_nonzero_coverage,
    blast_subtype="",
):
    post_trim_reads, post_trim_mean_len = read_fastq_stats(input_fastq)
    meta = get_hcv_reference_metadata(ref_path, ref_name)

    blast_subtype = str(blast_subtype or "").strip().upper()
    if blast_subtype and not re.fullmatch(r"[1-8][A-Z]", blast_subtype):
        blast_subtype = ""

    base_row = {
        "hcv_quality": "✕",
        "hcv_interpretation": "Null",
        "hcv_M28": " ",
        "hcv_Q30": " ",
        "hcv_L31": " ",
        "hcv_H58": " ",
        "hcv_Y93": " ",
        "hcv_low_coverage": " ",
        "hcv_post_trim": f"{post_trim_reads} | {post_trim_mean_len}",
        "hcv_mapping_distribution": f"{mapped_reads_1a} | {mapped_reads_1b}",
        "hcv_median_nonzero_coverage": median_nonzero_coverage,
        "hcv_subtype": blast_subtype or hcv_subtype_from_name(ref_name),
        "hcv_panel_subtype": "",
        "hcv_ns5a_start_used": "",
        "hcv_debug_details": [],
        "hcv_min_codon_depth": HCV_NS5A_MIN_CODON_DEPTH,
        "hcv_script_version": HCV_SCRIPT_VERSION,
    }

    if not ref_name or str(ref_name).upper() == "NA":
        base_row.update({
            "hcv_quality": "✕",
            "hcv_interpretation": "Null",
            "hcv_subtype": "Unknown",
            "hcv_low_coverage": "No valid HCV reference mapping found; resistance analysis not possible.",
            "hcv_ns5a_start_used": "not scanned",
            "hcv_debug_details": [],
        })
        return base_row

    if not meta:
        base_row["hcv_low_coverage"] = (
            "Missing HCV metadata: add subtype and ns5a_codon1_start "
            "to FASTA header or hcv_reference_metadata.tsv"
        )
        return base_row

    mapping_subtype = str(meta.get("subtype", hcv_subtype_from_name(ref_name))).upper().strip()

    # Sequence typing gates subtype-specific resistance interpretation.  A mapping
    # reference is still used for coordinates, but it must never silently select
    # a different subtype panel from the subtype identified by BLAST.
    if not blast_subtype:
        base_row.update({
            "hcv_quality": "Not Included",
            "hcv_interpretation": "Not Included",
            "hcv_low_coverage": "NS5A panel not run: no explicit HCV subtype was established by BLAST.",
            "hcv_ns5a_start_used": "not scanned",
        })
        return base_row

    base_row["hcv_subtype"] = blast_subtype

    if blast_subtype not in HCV_NS5A_PANEL_BY_SUBTYPE:
        base_row.update({
            "hcv_quality": "Not Included",
            "hcv_interpretation": "Not Included",
            "hcv_low_coverage": f"NS5A panel not run: BLAST subtype {blast_subtype} is not included in Aurora's validated HCV panel.",
            "hcv_ns5a_start_used": "not scanned",
        })
        return base_row

    if mapping_subtype != blast_subtype:
        base_row.update({
            "hcv_quality": "Not Included",
            "hcv_interpretation": "Not Included",
            "hcv_low_coverage": (
                f"NS5A panel not run: BLAST subtype {blast_subtype} does not match "
                f"mapping-reference subtype {mapping_subtype}."
            ),
            "hcv_ns5a_start_used": "not scanned",
        })
        return base_row

    base_row["hcv_panel_subtype"] = blast_subtype
    panel_meta = dict(meta)
    panel_meta["subtype"] = blast_subtype

    try:
        selected = choose_hcv_panel_result(bam_file, ref_name, ref_path, panel_meta)
        if not selected:
            base_row["hcv_low_coverage"] = "HCV resistance analysis failed: no coordinate candidates"
            return base_row

        position_results = selected.get("position_results", {})
        low_coverage = selected.get("low_coverage", [])
        interpretation = ("Ej bedömbar" if selected.get("error")
                          else determine_hcv_interpretation(position_results, low_coverage))

        low_coverage_text = ", ".join(low_coverage) if low_coverage else " "
        if selected.get("error"):
            low_coverage_text = f"HCV resistance analysis failed: {selected.get('error')}"

        base_row.update({
            "hcv_quality": "✕" if low_coverage or selected.get("error") else "✓",
            "hcv_interpretation": interpretation,
            "hcv_M28": format_hcv_position_call(position_results.get(28)),
            "hcv_Q30": format_hcv_position_call(position_results.get(30)),
            "hcv_L31": format_hcv_position_call(position_results.get(31)),
            "hcv_H58": format_hcv_position_call(position_results.get(58)),
            "hcv_Y93": format_hcv_position_call(position_results.get(93)),
            "hcv_low_coverage": low_coverage_text,
            "hcv_ns5a_start_used": selected.get("ns5a_start", ""),
            "hcv_debug_details": selected.get("details", []),
        })

    except Exception as e:
        base_row["hcv_quality"] = "✕"
        base_row["hcv_interpretation"] = "Null"
        base_row["hcv_low_coverage"] = f"HCV resistance analysis failed: {e}"

    return base_row
