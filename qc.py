"""Coverage and assay-specific QC calculations."""
from .alignment import calculate_multi_region_depth_stats, empty_depth_stats
from .config import (
    HBV_EXPECTED_GENOME_LENGTH,
    HBV_SYSTEM1_FULL_REGIONS,
    HBV_SYSTEM1_UNIQUE_REGIONS,
    HBV_SYSTEM2_FULL_REGIONS,
    HBV_SYSTEM2_UNIQUE_REGIONS,
    HBV_SYSTEM_OVERLAP_REGIONS,
)


def expected_target_length_for_qc(virus_type, detected_target_length):
    if virus_type == "Hepatitis B Virus":
        return HBV_EXPECTED_GENOME_LENGTH

    if virus_type == "Adenovirus":
        return detected_target_length or 0

    return detected_target_length or 0

def calculate_hbv_system_depth_stats(bam_file, ref_name):
    return {
        "system1_full": calculate_multi_region_depth_stats(bam_file, ref_name, HBV_SYSTEM1_FULL_REGIONS),
        "system2_full": calculate_multi_region_depth_stats(bam_file, ref_name, HBV_SYSTEM2_FULL_REGIONS),
        "system1_unique": calculate_multi_region_depth_stats(bam_file, ref_name, HBV_SYSTEM1_UNIQUE_REGIONS),
        "system2_unique": calculate_multi_region_depth_stats(bam_file, ref_name, HBV_SYSTEM2_UNIQUE_REGIONS),
        "overlap_only": calculate_multi_region_depth_stats(bam_file, ref_name, HBV_SYSTEM_OVERLAP_REGIONS),
    }

def empty_hbv_system_stats():
    return {
        "system1_full": empty_depth_stats(),
        "system2_full": empty_depth_stats(),
        "system1_unique": empty_depth_stats(),
        "system2_unique": empty_depth_stats(),
        "overlap_only": empty_depth_stats(),
    }

def hbv_system_qc_flags(hbv_system_stats):
    s1 = hbv_system_stats.get("system1_unique", empty_depth_stats())
    s2 = hbv_system_stats.get("system2_unique", empty_depth_stats())
    overlap = hbv_system_stats.get("overlap_only", empty_depth_stats())

    s1_median = s1.get("median_depth", 0)
    s2_median = s2.get("median_depth", 0)
    s1_10x = s1.get("positions_ge_10x_percent", 0)
    s2_10x = s2.get("positions_ge_10x_percent", 0)
    s1_30x = s1.get("positions_ge_30x_percent", 0)
    s2_30x = s2.get("positions_ge_30x_percent", 0)

    high_ok = (
        s1_median >= 30 and
        s2_median >= 30 and
        s1_10x >= 90 and
        s2_10x >= 90 and
        s1_30x >= 70 and
        s2_30x >= 70
    )

    medium_ok = (
        s1_median >= 10 and
        s2_median >= 10 and
        s1_10x >= 70 and
        s2_10x >= 70
    )

    warnings = []

    if s1_median < 10:
        warnings.append("System 1 unique region has low median depth <10x")
    elif s1_median < 30:
        warnings.append("System 1 unique region has moderate median depth <30x")

    if s2_median < 10:
        warnings.append("System 2 unique region has low median depth <10x")
    elif s2_median < 30:
        warnings.append("System 2 unique region has moderate median depth <30x")

    if s1_10x < 70:
        warnings.append(f"System 1 unique region has poor >=10x breadth ({s1_10x}%)")
    elif s1_10x < 90:
        warnings.append(f"System 1 unique region has partial >=10x breadth ({s1_10x}%)")

    if s2_10x < 70:
        warnings.append(f"System 2 unique region has poor >=10x breadth ({s2_10x}%)")
    elif s2_10x < 90:
        warnings.append(f"System 2 unique region has partial >=10x breadth ({s2_10x}%)")

    min_median = max(min(s1_median, s2_median), 1)
    max_median = max(s1_median, s2_median)
    balance_ratio = round(max_median / min_median, 2) if max_median > 0 else 0.0

    if s1_median == 0 or s2_median == 0:
        warnings.append("One HBV system unique region has zero median depth")
    elif balance_ratio >= 10:
        warnings.append(f"Strong imbalance between HBV systems based on unique regions ({balance_ratio}x)")
    elif balance_ratio >= 5:
        warnings.append(f"Moderate imbalance between HBV systems based on unique regions ({balance_ratio}x)")

    overlap_median = overlap.get("median_depth", 0)
    if overlap_median >= 30 and (s1_median < 30 or s2_median < 30):
        warnings.append("Overlap has good depth but at least one unique system region is weaker; QC uses unique regions")

    return {
        "high_ok": high_ok,
        "medium_ok": medium_ok,
        "warning": "; ".join(warnings) if warnings else "Both HBV systems acceptable by unique-region QC",
        "balance_ratio": balance_ratio,
    }

def determine_qc_verdict(
    virus_type,
    mapped_reads,
    coverage_percent,
    median_depth,
    n_percent,
    consensus_length,
    expected_length,
    hbv_system_high_ok=True,
    hbv_system_medium_ok=True,
    hbv_system_warning=""
):
    reasons = []

    if expected_length and expected_length > 0:
        if virus_type == "Adenovirus":
            # For short amplicons, sequence length alone can be misleading:
            # reference-coordinate trimming may retain insertions (raw length >
            # expected span), while a mostly-N sequence can still have the full
            # expected length. Define recovery as the fraction of the expected
            # amplicon represented by non-N bases and cap it at 100%.
            callable_bases = max(0.0, float(consensus_length) * (1.0 - (float(n_percent) / 100.0)))
            recovery_percent = round(min(100.0, (callable_bases / expected_length) * 100.0), 2)
        else:
            # Recovery is a percentage of the expected sequence span and should
            # never exceed 100%, even when the consensus contains insertions.
            recovery_percent = round(min(100.0, (consensus_length / expected_length) * 100.0), 2)
    else:
        recovery_percent = 0.0

    if virus_type == "Hepatitis B Virus":
        coverage_name = "HBV genome recovery"
    elif virus_type == "Adenovirus":
        coverage_name = "Adenovirus target/amplicon recovery"
    else:
        coverage_name = "target recovery"

    high = (
        mapped_reads >= 100 and
        coverage_percent >= 90 and
        median_depth >= 30 and
        n_percent <= 2 and
        recovery_percent >= 90 and
        hbv_system_high_ok
    )

    medium = (
        mapped_reads >= 50 and
        coverage_percent >= 70 and
        median_depth >= 10 and
        n_percent <= 10 and
        recovery_percent >= 70 and
        hbv_system_medium_ok
    )

    if high:
        verdict = "HIGH"
    elif medium:
        verdict = "MEDIUM"
    else:
        verdict = "LOW"

    if mapped_reads < 50:
        reasons.append("low mapped read count")
    elif mapped_reads < 100:
        reasons.append("moderate mapped read count")
    else:
        reasons.append("good mapped read count")

    if coverage_percent < 70:
        reasons.append(f"low {coverage_name}")
    elif coverage_percent < 90:
        reasons.append(f"partial {coverage_name}")
    else:
        reasons.append(f"strong {coverage_name}")

    if median_depth < 10:
        reasons.append("low median depth")
    elif median_depth < 30:
        reasons.append("moderate median depth")
    else:
        reasons.append("good median depth")

    if n_percent > 10:
        reasons.append("high N content")
    elif n_percent > 2:
        reasons.append("moderate N content")
    else:
        reasons.append("low N content")

    if expected_length and expected_length > 0:
        if recovery_percent < 70:
            if virus_type == "Adenovirus":
                reasons.append(f"low callable consensus recovery ({recovery_percent}%)")
            else:
                reasons.append(f"short consensus vs expected length ({recovery_percent}%)")
        elif recovery_percent < 90:
            if virus_type == "Adenovirus":
                reasons.append(f"partial callable consensus recovery ({recovery_percent}%)")
            else:
                reasons.append(f"partial consensus length recovery ({recovery_percent}%)")
        else:
            if virus_type == "Adenovirus":
                reasons.append(f"good callable consensus recovery ({recovery_percent}%)")
            else:
                reasons.append(f"good consensus length recovery ({recovery_percent}%)")

    if virus_type == "Hepatitis B Virus":
        if hbv_system_high_ok:
            reasons.append("HBV System 1 and System 2 unique-region depth supports HIGH QC")
        elif hbv_system_medium_ok:
            reasons.append("HBV System 1 and System 2 unique-region depth supports MEDIUM QC only")
        else:
            reasons.append("HBV System 1/System 2 unique-region depth fails medium threshold")

        if hbv_system_warning:
            reasons.append(hbv_system_warning)

    return verdict, "; ".join(reasons), recovery_percent
