"""Top-level sample processing and run orchestration."""
import os
import subprocess

from .blast import blast_consensus, empty_blast_result
from .cancellation import check_cancelled, run_cancellable
from .alignment import (
    calculate_depth_stats,
    calculate_reference_breadth_coverage_percent,
    detect_aligned_target_region_from_sam,
    detect_depth_supported_target_region,
    empty_depth_stats,
    evaluate_alignment_quality,
)
from .config import HBV_EXPECTED_GENOME_LENGTH, HBV_MIN_VARIANT_FREQUENCY, HCV_SCRIPT_VERSION
from .consensus import choose_medaka_model, generate_consensus, generate_iub_consensus_from_bam
from .external_tools import check_required_tools, run_ncbi_database_script, safe_remove
from .fasta import fasta_sequence_stats, amplicon_fasta_sequence_stats
from .filesystem import (
    allocate_run_identity,
    concatenate_fasta_files,
    ensure_dir,
    export_run_consensus_files,
    find_files_recursively,
    is_real_fasta,
    is_real_fastq_gz,
    run_samples_folder,
    split_fasta,
)
from .hbv import analyse_hbv_mutations, get_hbv_reference_metadata
from .cmv import analyse_cmv_mutations
from .cmv_qc import cmv_covered_depth_stats
from .cmv_qc import add_cmv_read_yield
from .hcv import (
    analyse_hcv_resistance,
    calculate_hcv_amplicon_breadth_coverage_percent,
    calculate_median_nonzero_depth,
    hcv_amplicon_consensus_stats,
    hcv_subtype_from_name,
)
from .qc import (
    calculate_hbv_system_depth_stats,
    determine_qc_verdict,
    empty_hbv_system_stats,
    expected_target_length_for_qc,
    hbv_system_qc_flags,
)
from .reporting_excel import log_hbv_mutation_results, log_hcv_resistance_results, log_reference_usage
from .reporting_html import generate_hcv_html_reports, generate_adenovirus_html_reports, write_hbv_html_report
from .filesystem import safe_report_filename


def add_hbv_system_stats_to_row(row, hbv_system_stats, hbv_flags):
    s1_full = hbv_system_stats.get("system1_full", empty_depth_stats())
    s2_full = hbv_system_stats.get("system2_full", empty_depth_stats())
    s1_unique = hbv_system_stats.get("system1_unique", empty_depth_stats())
    s2_unique = hbv_system_stats.get("system2_unique", empty_depth_stats())
    overlap = hbv_system_stats.get("overlap_only", empty_depth_stats())

    row["hbv_system1_full_median_depth"] = s1_full.get("median_depth", 0)
    row["hbv_system1_full_positions_ge_10x_percent"] = s1_full.get("positions_ge_10x_percent", 0)
    row["hbv_system1_full_positions_ge_30x_percent"] = s1_full.get("positions_ge_30x_percent", 0)

    row["hbv_system2_full_median_depth"] = s2_full.get("median_depth", 0)
    row["hbv_system2_full_positions_ge_10x_percent"] = s2_full.get("positions_ge_10x_percent", 0)
    row["hbv_system2_full_positions_ge_30x_percent"] = s2_full.get("positions_ge_30x_percent", 0)

    row["hbv_system1_unique_median_depth"] = s1_unique.get("median_depth", 0)
    row["hbv_system1_unique_positions_ge_10x_percent"] = s1_unique.get("positions_ge_10x_percent", 0)
    row["hbv_system1_unique_positions_ge_30x_percent"] = s1_unique.get("positions_ge_30x_percent", 0)

    row["hbv_system2_unique_median_depth"] = s2_unique.get("median_depth", 0)
    row["hbv_system2_unique_positions_ge_10x_percent"] = s2_unique.get("positions_ge_10x_percent", 0)
    row["hbv_system2_unique_positions_ge_30x_percent"] = s2_unique.get("positions_ge_30x_percent", 0)

    row["hbv_overlap_median_depth"] = overlap.get("median_depth", 0)
    row["hbv_overlap_positions_ge_10x_percent"] = overlap.get("positions_ge_10x_percent", 0)
    row["hbv_overlap_positions_ge_30x_percent"] = overlap.get("positions_ge_30x_percent", 0)

    row["hbv_system_balance_ratio"] = hbv_flags.get("balance_ratio", 0)
    row["hbv_system_qc_warning"] = hbv_flags.get("warning", "")

    return row

def process_sample(
    input_file,
    reference_files,
    output_folder,
    medaka_model,
    virus_type,
    final_consensus_list,
    result_rows,
    auto_detect_models=False
):
    try:
        base_name = os.path.splitext(os.path.basename(input_file))[0].replace(".fastq", "").replace(".gz", "")
        sample_dir = os.path.join(output_folder, base_name)
        ensure_dir(sample_dir, must_be_writable=True)
        selected_medaka_model = choose_medaka_model(input_file, medaka_model, auto_detect_models)

        best = {
            "ref_path": None,
            "ref_sn": None,
            "approx_coverage": -1.0,
            "sam_path": None,
            "mapped": 0,
            "aligned_query_bases": 0,
            "ref_length": 0
        }

        hcv_best_mapped_by_subtype = {"1A": 0, "1B": 0}

        for reference_file in reference_files:
            ref_name = os.path.splitext(os.path.basename(reference_file))[0]
            sam_tmp = os.path.join(sample_dir, f"{ref_name}.sam")

            cmd = ["minimap2", "-ax", "map-ont", reference_file, input_file]
            with open(sam_tmp, "w", encoding="utf-8") as sam_out:
                run_cancellable(cmd, check=True, stdout=sam_out)

            reference_selection_mode = (
                "support_first" if virus_type in ("Adenovirus", "Hepatitis C Virus") else "identity_first"
            )
            approx_cov, mapped, sn, aligned_query_bases, ref_length = evaluate_alignment_quality(
                sam_tmp,
                selection_mode=reference_selection_mode,
            )

            if virus_type == "Hepatitis C Virus":
                subtype = hcv_subtype_from_name(sn or ref_name)
                if subtype in hcv_best_mapped_by_subtype:
                    hcv_best_mapped_by_subtype[subtype] = max(hcv_best_mapped_by_subtype[subtype], mapped)

            if approx_cov > best["approx_coverage"]:
                safe_remove(best["sam_path"])
                best.update({
                    "ref_path": reference_file,
                    "ref_sn": sn,
                    "approx_coverage": approx_cov,
                    "sam_path": sam_tmp,
                    "mapped": mapped,
                    "aligned_query_bases": aligned_query_bases,
                    "ref_length": ref_length
                })
            else:
                safe_remove(sam_tmp)

        if not (best["ref_path"] and best["sam_path"] and best["ref_sn"] and best["mapped"] > 0):
            print(f"No suitable reference found for {input_file}. Skipping.")
            return

        # HBV genotype is metadata attached to the selected best reference.
        # Reference selection itself is unchanged.
        hbv_reference_metadata = {}
        hbv_genotype = ""
        if virus_type == "Hepatitis B Virus":
            hbv_reference_metadata = (
                get_hbv_reference_metadata(best["ref_path"], best["ref_sn"]) or {}
            )
            hbv_genotype = str(
                hbv_reference_metadata.get("genotype", "") or ""
            ).strip().upper()

            if hbv_genotype:
                print(
                    f"HBV genotype from selected reference {best['ref_sn']}: "
                    f"{hbv_genotype}"
                )
            else:
                print(
                    f"HBV genotype unavailable for selected reference {best['ref_sn']}. "
                    "Add a genotype column to hbv_reference_metadata.tsv."
                )

        detected_target_start, detected_target_end, detected_target_length = detect_aligned_target_region_from_sam(
            best["sam_path"],
            best["ref_sn"]
        )

        tmp_ref = os.path.join(sample_dir, f"{base_name}_best_reference.fasta")
        found = False

        with open(best["ref_path"], "r", encoding="utf-8", errors="replace") as ref_in, \
             open(tmp_ref, "w", encoding="utf-8") as ref_out:
            write = False
            for line in ref_in:
                if line.startswith(">"):
                    header = line[1:].strip()
                    name_only = header.split()[0]
                    write = (name_only == best["ref_sn"])
                if write:
                    ref_out.write(line)
                    found = True

        if not found:
            print(f"Could not find sequence '{best['ref_sn']}' in {best['ref_path']}. Using full reference instead.")
            tmp_ref = best["ref_path"]

        bam_sorted = os.path.join(sample_dir, f"{base_name}.sorted.bam")
        run_cancellable(["samtools", "sort", "-O", "bam", "-o", bam_sorted, best["sam_path"]], check=True)
        run_cancellable(["samtools", "index", bam_sorted], check=True)

        # Adenovirus target refinement: define the amplicon from the largest
        # depth-supported block, but do not require a fixed 100x at every
        # position. A fixed 100x gate was stricter than the downstream MEDIUM
        # QC definition and could discard otherwise usable amplicons.
        #
        # The threshold scales with the number of reads mapped to the selected
        # reference: 25% of mapped reads, bounded to 10..100x. High-depth
        # samples therefore retain the original ~100x stringency, while
        # moderate samples can still define an amplicon for honest QC. Very
        # low-read samples (for example 4-5 reads) still cannot reach 10x and
        # remain LOW without being artificially rescued.
        adeno_target_min_depth = None
        if virus_type == "Adenovirus":
            adeno_target_min_depth = max(
                10,
                min(100, int(round(float(best.get("mapped", 0) or 0) * 0.25))),
            )
            depth_target_start, depth_target_end, depth_target_length = (
                detect_depth_supported_target_region(
                    bam_file=bam_sorted,
                    ref_name=best["ref_sn"],
                    min_depth=adeno_target_min_depth,
                    max_gap=10,
                )
            )

            if depth_target_start and depth_target_end:
                detected_target_start = depth_target_start
                detected_target_end = depth_target_end
                detected_target_length = depth_target_length
                print(
                    "Adenovirus target refined using adaptive depth threshold "
                    f">={adeno_target_min_depth}x: "
                    f"{detected_target_start}-{detected_target_end} "
                    f"({detected_target_length} bp)."
                )
            else:
                detected_target_start = None
                detected_target_end = None
                detected_target_length = 0
                print(
                    "WARNING: No Adenovirus region reached the adaptive target "
                    f"threshold (>={adeno_target_min_depth}x); no Adenovirus "
                    "amplicon will be exported or BLASTed."
                )

        reference_coverage_percent = calculate_reference_breadth_coverage_percent(
            bam_file=bam_sorted,
            ref_name=best["ref_sn"],
            ref_length=best["ref_length"]
        )

        hbv_system_stats = empty_hbv_system_stats()
        hbv_flags = {
            "high_ok": True,
            "medium_ok": True,
            "warning": "",
            "balance_ratio": 0,
        }

        if virus_type == "Adenovirus":
            if detected_target_start and detected_target_end:
                depth_region = (detected_target_start, detected_target_end)
                coverage_type = "Detected Adenovirus depth-supported Target/Amplicon Coverage"
                expected_target_length = expected_target_length_for_qc(virus_type, detected_target_length)

                depth_stats = calculate_depth_stats(
                    bam_file=bam_sorted,
                    ref_name=best["ref_sn"],
                    region=depth_region
                )
                target_or_genome_coverage_percent = 100.0
            else:
                depth_region = None
                coverage_type = "No Adenovirus depth-supported Target/Amplicon"
                expected_target_length = 0
                depth_stats = empty_depth_stats()
                target_or_genome_coverage_percent = 0.0

        elif virus_type == "Hepatitis B Virus":
            depth_region = None
            coverage_type = "Approximate HBV Whole-Genome Coverage"
            expected_target_length = HBV_EXPECTED_GENOME_LENGTH

            depth_stats = calculate_depth_stats(
                bam_file=bam_sorted,
                ref_name=best["ref_sn"],
                region=None
            )

            target_or_genome_coverage_percent = reference_coverage_percent

            hbv_system_stats = calculate_hbv_system_depth_stats(
                bam_file=bam_sorted,
                ref_name=best["ref_sn"]
            )
            hbv_flags = hbv_system_qc_flags(hbv_system_stats)

        elif virus_type == "Hepatitis C Virus":
            depth_region = (
                (detected_target_start, detected_target_end)
                if detected_target_start and detected_target_end
                else None
            )
            coverage_type = "Detected HCV Target/Amplicon Coverage"
            expected_target_length = detected_target_length or 0

            depth_stats = calculate_depth_stats(
                bam_file=bam_sorted,
                ref_name=best["ref_sn"],
                region=depth_region
            )

            target_or_genome_coverage_percent = calculate_hcv_amplicon_breadth_coverage_percent(
                bam_file=bam_sorted,
                ref_name=best["ref_sn"],
                start_1based=detected_target_start,
                end_1based=detected_target_end,
            )

        else:
            depth_region = None
            coverage_type = "Reference Coverage"
            expected_target_length = expected_target_length_for_qc(virus_type, detected_target_length)

            depth_stats = calculate_depth_stats(
                bam_file=bam_sorted,
                ref_name=best["ref_sn"],
                region=None
            )

            target_or_genome_coverage_percent = reference_coverage_percent

        cmv_depth = {}
        if virus_type == "Cytomegalovirus (CMV)":
            cmv_depth = cmv_covered_depth_stats(bam_sorted, best["ref_sn"])
            # The min-to-max alignment span is NOT the number of covered bases.
            # Keep original span coordinates for compatibility, but report the
            # observed covered-base count as the detected sequenced length.
            detected_target_length = cmv_depth["cmv_covered_bases"]
            depth_stats = {**depth_stats,
                           "median_depth": cmv_depth["cmv_covered_median_depth"],
                           "mean_depth": cmv_depth["cmv_covered_mean_depth"]}
            coverage_type = "CMV whole-reference breadth (target depth excludes uncovered bases)"
            expected_target_length = 0  # unknown without amplicon BED coordinates
            print("CMV observed covered bases: {cmv_covered_bases}; median covered-base depth: {cmv_covered_median_depth}x; >=30x: {cmv_covered_bases_30x}".format(**cmv_depth))

        print(
            f"Best reference for {input_file}: {best['ref_sn']} "
            f"with {best['mapped']} mapped reads."
        )

        if virus_type == "Adenovirus":
            print(
                f"Detected adenovirus target region: "
                f"{detected_target_start}-{detected_target_end} "
                f"({detected_target_length} bp)."
            )
        elif virus_type == "Hepatitis C Virus":
            print(
                f"Detected HCV target amplicon: "
                f"{detected_target_start}-{detected_target_end} "
                f"({detected_target_length} bp)."
            )

        print(
            f"{coverage_type}: {target_or_genome_coverage_percent}% | "
            f"Median depth: {depth_stats.get('median_depth', 0)}x"
        )

        if virus_type == "Hepatitis B Virus":
            print(
                "HBV unique-region QC: "
                f"System 1 median {hbv_system_stats['system1_unique'].get('median_depth', 0)}x, "
                f"System 2 median {hbv_system_stats['system2_unique'].get('median_depth', 0)}x | "
                f"{hbv_flags.get('warning', '')}"
            )

        # If no Adenovirus region reaches the 100x support threshold, do not
        # produce/BLAST a whole-reference consensus. Record a LOW-QC result and
        # stop this sample here.
        if virus_type == "Adenovirus" and not (detected_target_start and detected_target_end):
            fasta_stats = {
                "consensus_length": 0,
                "n_count": 0,
                "n_percent": 100.0,
                "non_n_length": 0,
            }
            qc_verdict, qc_reason, consensus_recovery_percent = determine_qc_verdict(
                virus_type=virus_type,
                mapped_reads=best["mapped"],
                coverage_percent=0.0,
                median_depth=0,
                n_percent=100.0,
                consensus_length=0,
                expected_length=0,
            )
            qc_reason = (
                "no Adenovirus region reached the adaptive amplicon depth threshold; "
                + qc_reason
            )
            row = {
                "sample": base_name,
                "virus_type": virus_type,
                "best_reference": best["ref_sn"] or "NA",
                "hbv_genotype": "",
                "mapped_reads": best["mapped"],
                "reference_length": best["ref_length"],
                "detected_target_start": "",
                "detected_target_end": "",
                "detected_target_length": 0,
                "expected_target_length": 0,
                "coverage_type": coverage_type,
                "target_or_genome_coverage_percent": 0.0,
                "mean_depth": 0,
                "median_depth": 0,
                "min_depth": 0,
                "positions_ge_10x_percent": 0,
                "positions_ge_30x_percent": 0,
                "consensus_length": 0,
                "n_count": 0,
                "n_percent": 100.0,
                "consensus_recovery_percent": consensus_recovery_percent,
                "qc_verdict": qc_verdict,
                "qc_reason": qc_reason,
                "reference_selection_mode": (
                    "support_first" if virus_type in ("Adenovirus", "Hepatitis C Virus") else "identity_first"
                ),
                "adeno_target_min_depth": adeno_target_min_depth if virus_type == "Adenovirus" else "",
                "bam_file": bam_sorted,
                "input_fastq": input_file,
                "consensus_file": "",
                "iub_consensus_file": "",
                "reference_file": best["ref_path"],
                "medaka_model": "Ej körd",
            }
            row.update(empty_blast_result(
                status="Skipped",
                message=(
                    "No Adenovirus region reached the adaptive depth threshold; "
                    "consensus export and BLAST were skipped."
                ),
            ))
            result_rows.append(row)
            print(f"QC Verdict for {base_name}: {qc_verdict} ({qc_reason})")
            safe_remove(best["sam_path"])
            if tmp_ref != best["ref_path"]:
                safe_remove(tmp_ref)
            return

        # CMV pre-consensus indel assessment: observation only, before draft and Medaka.
        cmv_preconsensus_qc = {}
        if virus_type == "Cytomegalovirus (CMV)":
            try:
                from .cmv_preconsensus_qc import audit_alignment
                cmv_preconsensus_qc = audit_alignment(
                    bam_sorted, best["ref_path"], best["ref_sn"],
                    os.path.join(sample_dir, base_name),
                )
                print(f"[CMV PRE-CONSENSUS] {base_name}: "
                      f"{cmv_preconsensus_qc['frame_disrupting_events']} non-3n CIGAR indel events; "
                      "review-only QC saved before Medaka")
            except Exception as exc:
                cmv_preconsensus_qc = {"assessment": "QC_FAILED", "error": str(exc)}
                print(f"[CMV PRE-CONSENSUS] {base_name}: QC FAILED: {exc}")

        draft_consensus = os.path.join(sample_dir, f"{base_name}.draft.fasta")
        # For amplicon workflows, build the draft from exactly the detected
        # reference interval.  A whole-alignment samtools consensus can include
        # lower-depth flanks outside the depth-supported target (e.g. 291 nt for
        # a 272-nt Adenovirus amplicon), which then makes Medaka polish the wrong
        # span and breaks downstream length/QC logic.
        consensus_cmd = ["samtools", "consensus", "-f", "fasta"]
        if (
            virus_type in {"Adenovirus", "Hepatitis C Virus"}
            and detected_target_start
            and detected_target_end
            and best.get("ref_sn")
        ):
            consensus_cmd.extend([
                "-r",
                f"{best['ref_sn']}:{detected_target_start}-{detected_target_end}",
            ])
        consensus_cmd.extend(["-o", draft_consensus, bam_sorted])
        run_cancellable(consensus_cmd, check=True)

        if not (os.path.exists(draft_consensus) and os.path.getsize(draft_consensus) > 0):
            print(f"Draft consensus missing/empty for {base_name}; skipping Medaka.")

            fasta_stats = {
                "consensus_length": 0,
                "n_count": 0,
                "n_percent": 100.0,
                "non_n_length": 0,
            }

            qc_verdict, qc_reason, consensus_recovery_percent = determine_qc_verdict(
                virus_type=virus_type,
                mapped_reads=best["mapped"],
                coverage_percent=target_or_genome_coverage_percent,
                median_depth=depth_stats.get("median_depth", 0),
                n_percent=fasta_stats["n_percent"],
                consensus_length=fasta_stats["consensus_length"],
                expected_length=expected_target_length,
                hbv_system_high_ok=hbv_flags.get("high_ok", True),
                hbv_system_medium_ok=hbv_flags.get("medium_ok", True),
                hbv_system_warning=hbv_flags.get("warning", "")
            )

            row = {
                "sample": base_name,
                "virus_type": virus_type,
                "best_reference": best["ref_sn"] or "NA",
                "hbv_genotype": hbv_genotype if virus_type == "Hepatitis B Virus" else "",
                "mapped_reads": best["mapped"],
                "reference_length": best["ref_length"],
                "detected_target_start": detected_target_start or "",
                "detected_target_end": detected_target_end or "",
                "detected_target_length": detected_target_length,
                "expected_target_length": expected_target_length,
                "coverage_type": coverage_type,
                "target_or_genome_coverage_percent": target_or_genome_coverage_percent,
                "mean_depth": depth_stats.get("mean_depth", 0),
                "median_depth": depth_stats.get("median_depth", 0),
                "min_depth": depth_stats.get("min_depth", 0),
                "positions_ge_10x_percent": depth_stats.get("positions_ge_10x_percent", 0),
                "positions_ge_30x_percent": depth_stats.get("positions_ge_30x_percent", 0),
                "consensus_length": fasta_stats["consensus_length"],
                "n_count": fasta_stats["n_count"],
                "n_percent": fasta_stats["n_percent"],
                "consensus_recovery_percent": consensus_recovery_percent,
                "qc_verdict": qc_verdict,
                "qc_reason": qc_reason,
                "reference_selection_mode": (
                    "support_first" if virus_type in ("Adenovirus", "Hepatitis C Virus") else "identity_first"
                ),
                "adeno_target_min_depth": adeno_target_min_depth if virus_type == "Adenovirus" else "",
                "bam_file": bam_sorted,
                "input_fastq": input_file,
                "consensus_file": draft_consensus,
                "reference_file": best["ref_path"],
                "medaka_model": "Ej körd",
            }

            # No usable consensus was produced, so automatic local BLAST is skipped
            # but the standard BLAST result fields are still present in the row.
            row.update(empty_blast_result(
                status="Skipped",
                message="Consensus FASTA missing or empty; automatic local BLAST was not run.",
            ))

            if virus_type == "Hepatitis B Virus":
                row = add_hbv_system_stats_to_row(row, hbv_system_stats, hbv_flags)
                row.update(analyse_hbv_mutations(base_name, bam_sorted, best["ref_sn"], best["ref_path"]))

            if virus_type == "Hepatitis C Virus":
                row.update(analyse_hcv_resistance(
                    sample=base_name,
                    input_fastq=input_file,
                    bam_file=bam_sorted,
                    ref_name=best["ref_sn"],
                    ref_path=best["ref_path"],
                    mapped_reads_1a=hcv_best_mapped_by_subtype.get("1A", 0),
                    mapped_reads_1b=hcv_best_mapped_by_subtype.get("1B", 0),
                    median_nonzero_coverage=calculate_median_nonzero_depth(
                        bam_sorted, best["ref_sn"], depth_region
                    ),
                    blast_subtype=row.get("blast_subtype", ""),
                ))

            if virus_type == "Cytomegalovirus (CMV)":
                row.update(cmv_depth)
                add_cmv_read_yield(row)
                row["cmv_qc_status"] = "NO_CONSENSUS"
                row["cmv_qc_note"] = "Draft consensus missing or empty; no CMV target-level QC or protein assessment possible."
            result_rows.append(row)

            safe_remove(best["sam_path"])
            if tmp_ref != best["ref_path"]:
                safe_remove(tmp_ref)
            return

        iub_consensus_path = generate_iub_consensus_from_bam(
            bam_file=bam_sorted,
            ref_name=best["ref_sn"],
            ref_path=best["ref_path"],
            output_folder=sample_dir,
            sample_name=base_name,
            variant_frequency=HBV_MIN_VARIANT_FREQUENCY if virus_type == "Hepatitis B Virus" else 0.15,
            trim_start_1based=(
                detected_target_start
                if virus_type in {"Adenovirus", "Hepatitis C Virus"}
                else None
            ),
            trim_end_1based=(
                detected_target_end
                if virus_type in {"Adenovirus", "Hepatitis C Virus"}
                else None
            ),
        )

        # Medaka may emit multiple contigs when the draft/reference is multi-FASTA.
        # Always retain only the contig corresponding to the selected mapping
        # reference. For Adenovirus and HCV, also crop the exported consensus
        # to the detected read-supported amplicon interval. Internal Ns are preserved.
        medaka_result = generate_consensus(
            input_file,
            draft_consensus,
            selected_medaka_model,
            sample_dir,
            base_name,
            ref_name=best["ref_sn"],
            trim_start_1based=(
                detected_target_start
                if virus_type in {"Adenovirus", "Hepatitis C Virus"}
                else None
            ),
            trim_end_1based=(
                detected_target_end
                if virus_type in {"Adenovirus", "Hepatitis C Virus"}
                else None
            ),
        )
        if virus_type == "Cytomegalovirus (CMV)":
            try:
                from .cmv_preconsensus_qc import write_consensus_comparison
                cmv_preconsensus_qc["draft_medaka_summary_path"] = write_consensus_comparison(
                    draft_consensus, medaka_result, os.path.join(sample_dir, base_name))
            except Exception as exc:
                cmv_preconsensus_qc["comparison_error"] = str(exc)

        final_consensus_path_for_analysis = medaka_result if medaka_result else draft_consensus
        consensus_source = "Medaka" if medaka_result else "samtools draft consensus"
        consensus_fallback_reason = ""

        # Adenovirus safety gate:
        # A successfully-created Medaka FASTA is not necessarily a usable
        # consensus. In validation we observed short, deeply covered amplicons
        # where Medaka returned an amplicon consisting entirely of Ns while the
        # BAM-derived IUPAC consensus contained high-confidence A/C/G/T calls at
        # every position. Do not let such a technical polishing failure erase a
        # valid read-supported sequence.
        #
        # The fallback is intentionally conservative: Medaka must fail the
        # ordinary MEDIUM-QC N threshold (>10% N), while the IUPAC consensus
        # must have <=10% N, <=2% ambiguity, >=90% unambiguous A/C/G/T, and
        # approximately the expected amplicon length.
        if virus_type == "Adenovirus" and iub_consensus_path:
            medaka_stats_for_gate = amplicon_fasta_sequence_stats(
                final_consensus_path_for_analysis,
                detected_target_start,
                detected_target_end,
            )
            iub_stats_for_gate = amplicon_fasta_sequence_stats(
                iub_consensus_path,
                detected_target_start,
                detected_target_end,
            )

            expected_len = int(detected_target_length or 0)
            iub_len = int(iub_stats_for_gate.get("consensus_length", 0))
            length_ok = (
                expected_len > 0
                and iub_len >= max(1, int(round(expected_len * 0.90)))
                and iub_len <= expected_len + 20
            )
            medaka_poor = (
                float(medaka_stats_for_gate.get("n_percent", 100.0)) > 10.0
                or int(medaka_stats_for_gate.get("consensus_length", 0))
                    < max(1, int(round(expected_len * 0.90)))
            )
            iub_good = (
                length_ok
                and float(iub_stats_for_gate.get("n_percent", 100.0)) <= 10.0
                and float(iub_stats_for_gate.get("ambiguity_percent", 100.0)) <= 2.0
                and float(iub_stats_for_gate.get("acgt_percent", 0.0)) >= 90.0
            )

            if medaka_poor and iub_good:
                final_consensus_path_for_analysis = iub_consensus_path
                consensus_source = "BAM-derived IUPAC consensus (Medaka fallback)"
                consensus_fallback_reason = (
                    "Medaka consensus failed sequence-quality validation "
                    f"({medaka_stats_for_gate.get('n_percent', 100.0)}% N); "
                    "the BAM-derived amplicon consensus passed the fallback criteria."
                )
                print(
                    f"[ADENOVIRUS CONSENSUS FALLBACK] {base_name}: "
                    f"{consensus_fallback_reason}"
                )
            elif medaka_poor:
                # Keep the LOW/failed Medaka result rather than silently rescuing
                # it with a weak BAM consensus, but record exactly why fallback
                # was rejected so problematic samples can be audited.
                consensus_fallback_reason = (
                    "Medaka consensus failed sequence-quality validation, but "
                    "the BAM-derived consensus did not meet fallback criteria: "
                    f"length {iub_len}/{expected_len} nt, "
                    f"N {iub_stats_for_gate.get('n_percent', 100.0)}%, "
                    f"IUPAC ambiguity {iub_stats_for_gate.get('ambiguity_percent', 100.0)}%, "
                    f"ACGT {iub_stats_for_gate.get('acgt_percent', 0.0)}%."
                )
                print(
                    f"[ADENOVIRUS FALLBACK NOT USED] {base_name}: "
                    f"{consensus_fallback_reason}"
                )

        # Export/collect the consensus actually used for downstream analysis.
        if final_consensus_path_for_analysis and os.path.isfile(final_consensus_path_for_analysis):
            final_consensus_list.append(final_consensus_path_for_analysis)

        fasta_stats = fasta_sequence_stats(final_consensus_path_for_analysis)

        if virus_type == "Adenovirus":
            fasta_stats = amplicon_fasta_sequence_stats(
                final_consensus_path_for_analysis,
                detected_target_start,
                detected_target_end,
            )

        elif virus_type == "Hepatitis C Virus":
            fasta_stats = hcv_amplicon_consensus_stats(
                final_consensus_path_for_analysis,
                detected_target_start,
                detected_target_end,
            )

        qc_verdict, qc_reason, consensus_recovery_percent = determine_qc_verdict(
            virus_type=virus_type,
            mapped_reads=best["mapped"],
            coverage_percent=target_or_genome_coverage_percent,
            median_depth=depth_stats.get("median_depth", 0),
            n_percent=fasta_stats["n_percent"],
            consensus_length=fasta_stats["consensus_length"],
            expected_length=expected_target_length,
            hbv_system_high_ok=hbv_flags.get("high_ok", True),
            hbv_system_medium_ok=hbv_flags.get("medium_ok", True),
            hbv_system_warning=hbv_flags.get("warning", "")
        )

        row = {
            "sample": base_name,
            "virus_type": virus_type,
            "best_reference": best["ref_sn"] or "NA",
            "hbv_genotype": hbv_genotype if virus_type == "Hepatitis B Virus" else "",
            "mapped_reads": best["mapped"],
            "reference_length": best["ref_length"],
            "detected_target_start": detected_target_start or "",
            "detected_target_end": detected_target_end or "",
            "detected_target_length": detected_target_length,
            "expected_target_length": expected_target_length,
            "coverage_type": coverage_type,
            "target_or_genome_coverage_percent": target_or_genome_coverage_percent,
            "mean_depth": depth_stats.get("mean_depth", 0),
            "median_depth": depth_stats.get("median_depth", 0),
            "min_depth": depth_stats.get("min_depth", 0),
            "positions_ge_10x_percent": depth_stats.get("positions_ge_10x_percent", 0),
            "positions_ge_30x_percent": depth_stats.get("positions_ge_30x_percent", 0),
            "consensus_length": fasta_stats["consensus_length"],
            "n_count": fasta_stats["n_count"],
            "n_percent": fasta_stats["n_percent"],
            "consensus_recovery_percent": consensus_recovery_percent,
            "qc_verdict": qc_verdict,
            "qc_reason": qc_reason,
            "bam_file": bam_sorted,
            "input_fastq": input_file,
            "consensus_file": final_consensus_path_for_analysis,
            "medaka_consensus_file": medaka_result or "",
            "iub_consensus_file": iub_consensus_path or "",
            "consensus_source": consensus_source,
            "consensus_fallback_reason": consensus_fallback_reason,
            "reference_selection_mode": (
                "support_first" if virus_type in ("Adenovirus", "Hepatitis C Virus") else "identity_first"
            ),
            "adeno_target_min_depth": adeno_target_min_depth if virus_type == "Adenovirus" else "",
            "iub_ambiguity_count": fasta_stats.get("ambiguity_count", 0),
            "iub_ambiguity_percent": fasta_stats.get("ambiguity_percent", 0.0),
            "reference_file": best["ref_path"],
            "medaka_model": selected_medaka_model,
        }

        if virus_type == "Cytomegalovirus (CMV)":
            row.update(cmv_depth)
            row["cmv_preconsensus_qc_status"] = cmv_preconsensus_qc.get("assessment", "NOT_RUN")
            row["cmv_preconsensus_indel_events"] = cmv_preconsensus_qc.get("indel_events", 0)
            row["cmv_preconsensus_frameshift_candidate_events"] = cmv_preconsensus_qc.get("frame_disrupting_events", 0)
            row["cmv_preconsensus_qc_json"] = cmv_preconsensus_qc.get("json_path", "")
            row["cmv_preconsensus_qc_tsv"] = cmv_preconsensus_qc.get("tsv_path", "")
            row["cmv_draft_medaka_summary"] = cmv_preconsensus_qc.get("draft_medaka_summary_path", "")
            add_cmv_read_yield(row)
            # Existing QC rules were designed around genome-wide recovery.
            # Do not silently reinterpret their verdict as target QC.
            row["cmv_qc_status"] = ("NO_COVERAGE" if not cmv_depth.get("cmv_covered_bases")
                                    else "TARGET_QC_NOT_VALIDATED")
            row["cmv_qc_note"] = (
                "Whole-genome QC verdict is retained separately and is not an amplicon QC verdict. "
                "Covered-base median excludes zero-depth genome positions. "
                "Expected amplicon intervals are unknown, so target completeness and clinical QC cannot be determined."
            )
            row["cmv_assessment_status"] = (
                "INSUFFICIENT_DATA" if cmv_depth.get("cmv_covered_bases_30x", 0) == 0
                else "EXPLORATORY_GENE_SPECIFIC_ONLY"
            )
            row["cmv_legacy_whole_genome_qc"] = qc_verdict
            print(f"[CMV QC] {base_name}: covered-base median {cmv_depth['cmv_covered_median_depth']}x; "
                  f"{cmv_depth['cmv_covered_bases_30x']} bases >=30x; "
                  "target QC not validated (no assay BED).")

        # Automatic local NCBI BLAST for every successfully generated consensus.
        # blast_consensus() routes HBV/HCV/Adenovirus/CMV to their dedicated
        # small databases and uses the full nt_viruses database for Other.
        print(f"Running automatic local NCBI BLAST for {base_name} ({virus_type})...")
        row.update(blast_consensus(final_consensus_path_for_analysis, virus_type))
        print(
            f"[BLAST] {base_name}: status={row.get('blast_status', '')}; "
            f"database={row.get('blast_database', '')}; "
            f"hits={len(row.get('blast_hits') or [])}; "
            f"subtype={row.get('blast_subtype', '')}; "
            f"message={row.get('blast_message', '')}"
        )

        if virus_type == "Hepatitis B Virus":
            row = add_hbv_system_stats_to_row(row, hbv_system_stats, hbv_flags)
            row.update(analyse_hbv_mutations(base_name, bam_sorted, best["ref_sn"], best["ref_path"]))

        if virus_type == "Cytomegalovirus (CMV)":
            # CMV resistance screening uses a fixed Merlin coordinate system;
            # the original best reference still controls strain/consensus selection.
            from .cmv import analyse_cmv_standard
            try:
                cmv_result = analyse_cmv_standard(
                    base_name, input_file, os.path.dirname(best["ref_path"]), sample_dir
                )
                row.update(cmv_result)
                print(f"[CMV STANDARD] {base_name}: UL54/UL97 screening on NC_006273.2")
            except Exception as exc:
                # Never report missing calls as mutation-negative.
                row.update({
                    "cmv_protein_mutation_matches": [],
                    "cmv_protein_variant_calls": [],
                    "cmv_gene_coverage": [],
                    "cmv_published_mutation_status": f"NOT ASSESSED: {exc}",
                    "cmv_interpretation": f"NOT ASSESSED: {exc}",
                    "cmv_resistance_mapping_status": "FAILED",
                    "cmv_variant_calls": [],
                    "cmv_indel_calls": [],
                    "cmv_assessed_positions": 0,
                })
                print(f"[CMV STANDARD] {base_name}: NOT ASSESSED: {exc}")

            # Independent CMV-only reference mapping / BCFtools caller.
            # It never overwrites existing Aurora/Medaka calls and its failure
            # is written separately rather than terminating other samples.
            from .cmv_bcftools import optional_run as cmv_bcftools_review
            row.update(cmv_bcftools_review(
                base_name, input_file, os.path.dirname(best["ref_path"]), sample_dir
            ))

        if virus_type == "Hepatitis C Virus":
            row.update(analyse_hcv_resistance(
                sample=base_name,
                input_fastq=input_file,
                bam_file=bam_sorted,
                ref_name=best["ref_sn"],
                ref_path=best["ref_path"],
                mapped_reads_1a=hcv_best_mapped_by_subtype.get("1A", 0),
                mapped_reads_1b=hcv_best_mapped_by_subtype.get("1B", 0),
                median_nonzero_coverage=calculate_median_nonzero_depth(
                        bam_sorted, best["ref_sn"], depth_region
                    ),
                blast_subtype=row.get("blast_subtype", ""),
            ))

        result_rows.append(row)

        print(f"QC Verdict for {base_name}: {qc_verdict} ({qc_reason})")

        safe_remove(best["sam_path"])
        safe_remove(draft_consensus)

        if tmp_ref != best["ref_path"]:
            safe_remove(tmp_ref)

    except subprocess.CalledProcessError as e:
        print(f"Subprocess error while processing {input_file}: {e}")
    except Exception as e:
        print(f"Unexpected error in processing sample {input_file}: {e}")

def run_full_pipeline(input_folder, output_folder, reference_folder, medaka_model, virus_type, run_ncbi=False, auto_detect_models=False):
    # CMV-only macOS tool resolution: Medaka's mini_align inherits this PATH.
    # Do not alter executable resolution for HBV, HCV or Adenovirus runs.
    if virus_type == "Cytomegalovirus (CMV)" and os.path.isfile("/opt/homebrew/bin/minimap2"):
        import shutil
        native_bin = "/opt/homebrew/bin"
        os.environ["PATH"] = native_bin + os.pathsep + os.environ.get("PATH", "")
        resolved = shutil.which("minimap2")
        print(f"[CMV TOOLS] minimap2: {resolved}", flush=True)
        if resolved != "/opt/homebrew/bin/minimap2":
            print("[CMV TOOLS] WARNING: native minimap2 was not selected", flush=True)
    print("\nRunning full pipeline...\n")
    print(f"Aurora script file: {os.path.abspath(__file__)}")
    print(f"HCV script version: {HCV_SCRIPT_VERSION if virus_type == 'Hepatitis C Virus' else 'not HCV run'}")
    print(f"Virus type: {virus_type}")
    print(f"Input: {input_folder}")
    print(f"Output: {output_folder}")
    print(f"Reference folder: {reference_folder}")
    print(f"Medaka fallback model: {medaka_model}")
    print(f"Auto-detect Medaka models: {auto_detect_models}")
    print("Keep BAM files: True\n")

    try:
        check_required_tools()
    except Exception as e:
        print(f"[ERROR] {e}")
        return

    ensure_dir(output_folder, must_be_writable=True)
    try:
        run_date, run_number, run_tag, run_marker_path = allocate_run_identity(output_folder)
    except FileExistsError:
        # Extremely unlikely collision (for example two launches at the same time).
        run_date, run_number, run_tag, run_marker_path = allocate_run_identity(output_folder)
    print(f"Run identifier: {run_tag}")
    print(f"Run marker: {run_marker_path}\n")

    input_files = find_files_recursively(input_folder, is_real_fastq_gz)
    reference_files = [
        os.path.join(reference_folder, f)
        for f in os.listdir(reference_folder)
        if is_real_fasta(os.path.join(reference_folder, f))
    ]

    if not input_files:
        print("No .fastq.gz input files found.")
        return

    if not reference_files:
        print("No reference .fa/.fasta files found.")
        return

    result_rows = []
    final_consensus_list = []
    samples_output_folder = run_samples_folder(output_folder, run_tag)
    ensure_dir(samples_output_folder, must_be_writable=True)

    for input_file in input_files:
        check_cancelled()
        print(f"\n--- Processing sample: {input_file} ---")
        process_sample(
            input_file=input_file,
            reference_files=reference_files,
            output_folder=samples_output_folder,
            medaka_model=medaka_model,
            virus_type=virus_type,
            final_consensus_list=final_consensus_list,
            result_rows=result_rows,
            auto_detect_models=auto_detect_models
        )

    check_cancelled()
    if result_rows:
        if virus_type == "Hepatitis C Virus":
            log_hcv_resistance_results(result_rows, output_folder, run_tag)
            log_reference_usage(result_rows, output_folder, virus_type, run_tag)
            generate_hcv_html_reports(result_rows, output_folder, run_tag)
        elif virus_type == "Hepatitis B Virus":
            log_hbv_mutation_results(result_rows, output_folder, run_tag)
            log_reference_usage(result_rows, output_folder, virus_type, run_tag)
            hbv_html_dir = os.path.join(samples_output_folder, "HBV_HTML_reports")
            ensure_dir(hbv_html_dir, must_be_writable=True)
            for row in result_rows:
                sample = safe_report_filename(str(row.get("sample", "sample") or "sample"))
                report_path = os.path.join(hbv_html_dir, f"{sample}_{run_tag}_HBV.html")
                write_hbv_html_report(row, report_path, run_tag)
                row["hbv_html_report"] = report_path
            print(f"HBV HTML reports saved to {hbv_html_dir}")
        elif virus_type == "Cytomegalovirus (CMV)":
            log_reference_usage(result_rows, output_folder, virus_type, run_tag)
            from .cmv_reporting import generate_cmv_reports
            generate_cmv_reports(result_rows, output_folder, run_tag)
        else:
            log_reference_usage(result_rows, output_folder, virus_type, run_tag)
            if virus_type == "Adenovirus":
                generate_adenovirus_html_reports(result_rows, output_folder, run_tag)
    else:
        print("No reference usage data recorded.")

    if result_rows:
        export_run_consensus_files(result_rows, output_folder, run_tag)

    safe_virus_name = virus_type.replace(" ", "_")
    final_fasta_path = os.path.join(
        output_folder,
        f"{safe_virus_name}_{run_tag}_{len(final_consensus_list)}_consensus.fasta"
    )

    if final_consensus_list:
        success = concatenate_fasta_files(final_consensus_list, final_fasta_path)
        if success and run_ncbi:
            split_folder = os.path.join(samples_output_folder, "split_fasta")
            split_fasta(final_fasta_path, split_folder)
            run_ncbi_database_script(split_folder)
    else:
        print("No consensus files produced; skipping concatenation and optional NCBI step.")