"""CMV-only workflow hooks. Independent of HBV/HCV/Adenovirus dispatch.

Clinical research only. This module intentionally preserves the October 20
CMV mutation caller and does not change Medaka or variant thresholds.
"""
import os


def before_consensus(bam_sorted, selected_reference_fasta, selected_contig, output_prefix):
    """Read-only alignment audit; do not edit reads or draft consensus."""
    try:
        from .cmv_preconsensus_qc import audit_alignment
        result = audit_alignment(bam_sorted, selected_reference_fasta, selected_contig, output_prefix)
        print(f"[CMV PRE-CONSENSUS] {os.path.basename(output_prefix)}: "
              f"{result['frame_disrupting_events']} non-3n CIGAR events (review only)")
        return result
    except Exception as exc:
        print(f"[CMV PRE-CONSENSUS] QC FAILED: {exc}")
        return {"assessment": "QC_FAILED", "error": str(exc)}


def after_polishing(qc_record, draft_fasta, medaka_fasta, output_prefix):
    """Preserve existing draft/Medaka length comparison, not a local variant caller."""
    try:
        from .cmv_preconsensus_qc import write_consensus_comparison
        qc_record["draft_medaka_summary_path"] = write_consensus_comparison(
            draft_fasta, medaka_fasta, output_prefix)
    except Exception as exc:
        qc_record["comparison_error"] = str(exc)
    return qc_record


def resistance_screen(sample, fastq, reference_directory, sample_dir):
    """Fixed Merlin remapping and existing exploratory UL54/UL97 caller."""
    from .cmv import analyse_cmv_standard
    try:
        result = analyse_cmv_standard(sample, fastq, reference_directory, sample_dir)
        print(f"[CMV STANDARD] {sample}: UL54/UL97 screening on NC_006273.2")
        from .cmv_viralrecon_bridge import safe_optional_review
        from .cmv_independent_caller import optional_run
        result = optional_run(sample, fastq, reference_directory, sample_dir, result)
        return safe_optional_review(sample, result, sample_dir)
    except Exception as exc:
        print(f"[CMV STANDARD] {sample}: NOT ASSESSED: {exc}")
        return {
            "cmv_protein_mutation_matches": [], "cmv_protein_variant_calls": [],
            "cmv_gene_coverage": [], "cmv_published_mutation_status": f"NOT ASSESSED: {exc}",
            "cmv_interpretation": f"NOT ASSESSED: {exc}",
            "cmv_resistance_mapping_status": "FAILED", "cmv_variant_calls": [],
            "cmv_indel_calls": [], "cmv_assessed_positions": 0,
        }
