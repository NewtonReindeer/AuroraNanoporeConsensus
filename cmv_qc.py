"""CMV-specific exploratory read-yield classification.

Read-yield categories are not an assay-wide pass/fail decision and do not
replace target breadth, gene coverage or per-codon evidence requirements.
"""


def classify_cmv_read_yield(mapped_reads):
    """Classify mapped CMV reads: HIGH >=100, MEDIUM >=50, LOW >=10."""
    try:
        count = int(mapped_reads)
    except (TypeError, ValueError):
        return "UNKNOWN"
    if count < 0:
        return "UNKNOWN"
    if count >= 100:
        return "HIGH"
    if count >= 50:
        return "MEDIUM"
    if count >= 10:
        return "LOW"
    return "INSUFFICIENT"


def add_cmv_read_yield(row):
    row["cmv_read_yield_qc"] = classify_cmv_read_yield(row.get("mapped_reads"))
    row["cmv_read_yield_note"] = (
        "Mapped-read yield only (HIGH >=100, MEDIUM 50-99, LOW 10-49, "
        "INSUFFICIENT 0-9). Not overall assay QC; interpret gene coverage, "
        "depth and variant evidence separately."
    )
    return row


# Covered-base depth metrics (formerly cmv_depth.py)
"""CMV-only covered-base depth statistics (no expected amplicon inference)."""
from statistics import median

def cmv_covered_depth_stats(bam_path, reference_name, min_mapq=20, min_baseq=20):
    import pysam
    depths = []
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for col in bam.pileup(reference_name, stepper='samtools', truncate=True,
                              min_mapping_quality=min_mapq, min_base_quality=min_baseq,
                              max_depth=100000):
            depth = sum(1 for p in col.pileups if not p.is_del and not p.is_refskip
                        and not p.alignment.is_secondary and not p.alignment.is_supplementary
                        and not p.alignment.is_duplicate and not p.alignment.is_qcfail
                        and p.query_position is not None)
            if depth > 0:
                depths.append(depth)
    covered = len(depths)
    ge30 = sum(d >= 30 for d in depths)
    return {'cmv_covered_bases': covered,
            'cmv_covered_median_depth': float(median(depths)) if depths else 0.0,
            'cmv_covered_mean_depth': round(sum(depths)/covered, 2) if covered else 0.0,
            'cmv_covered_bases_30x': ge30,
            'cmv_covered_bases_30x_percent': round(100*ge30/covered, 2) if covered else 0.0}
