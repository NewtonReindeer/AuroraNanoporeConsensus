"""SAM/BAM alignment parsing and depth calculations."""
import os
import re
import statistics
import subprocess


LN_RE = re.compile(r"\bLN:(\d+)")

SN_RE = re.compile(r"\bSN:([^\t]+)")

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

def reference_consumed_by_cigar(cigar):
    ref_consumed = 0
    for length_str, op in CIGAR_RE.findall(cigar):
        length = int(length_str)
        if op in ["M", "=", "X", "D", "N"]:
            ref_consumed += length
    return ref_consumed

def query_aligned_consumed_by_cigar(cigar):
    query_consumed = 0
    for length_str, op in CIGAR_RE.findall(cigar):
        length = int(length_str)
        if op in ["M", "=", "X", "I"]:
            query_consumed += length
    return query_consumed

def evaluate_alignment_quality(output_sam, selection_mode="identity_first"):
    """Score a mapping reference from primary alignments.

    ``selection_mode`` controls how competing references are ranked:

    - ``identity_first`` preserves the HBV-oriented behavior introduced to
      distinguish very similar whole-genome references.
    - ``support_first`` is used for Adenovirus amplicons, where a reference
      supported by hundreds of concordant reads must not lose to a different
      reference represented by only one or a few perfect reads.

    In both modes the SAM NM edit-distance tag estimates matching query bases,
    and secondary/supplementary/QC-failed alignments are ignored.
    """
    ref_lengths = {}
    ref_mapped_reads = {}
    ref_aligned_query_bases = {}
    ref_matching_query_bases = {}
    ref_mapq_sum = {}

    try:
        with open(output_sam, "r", encoding="utf-8", errors="replace") as sam_file:
            for line in sam_file:
                if line.startswith("@SQ"):
                    ln_m = LN_RE.search(line)
                    sn_m = SN_RE.search(line)
                    if ln_m and sn_m:
                        ref_lengths[sn_m.group(1)] = int(ln_m.group(1))
                    continue

                if not line or line.startswith("@"):
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 11:
                    continue

                try:
                    flag = int(fields[1])
                    mapq = int(fields[4])
                except ValueError:
                    continue

                ref_name = fields[2]
                cigar = fields[5]

                # Skip unmapped, secondary, supplementary and QC-failed records.
                if flag & (0x4 | 0x100 | 0x200 | 0x800):
                    continue
                if ref_name == "*" or cigar == "*":
                    continue

                query_aligned = query_aligned_consumed_by_cigar(cigar)
                if query_aligned <= 0:
                    continue

                nm = 0
                for tag in fields[11:]:
                    if tag.startswith("NM:i:"):
                        try:
                            nm = int(tag[5:])
                        except ValueError:
                            nm = 0
                        break

                matching_bases = max(query_aligned - nm, 0)

                ref_mapped_reads[ref_name] = ref_mapped_reads.get(ref_name, 0) + 1
                ref_aligned_query_bases[ref_name] = ref_aligned_query_bases.get(ref_name, 0) + query_aligned
                ref_matching_query_bases[ref_name] = ref_matching_query_bases.get(ref_name, 0) + matching_bases
                ref_mapq_sum[ref_name] = ref_mapq_sum.get(ref_name, 0) + mapq

        best_ref = None
        best_rank = None

        for ref_name, matching_bases in ref_matching_query_bases.items():
            aligned_bases = ref_aligned_query_bases.get(ref_name, 0)
            identity = matching_bases / aligned_bases if aligned_bases else 0.0
            mapped_reads = ref_mapped_reads.get(ref_name, 0)
            mean_mapq = ref_mapq_sum.get(ref_name, 0) / mapped_reads if mapped_reads else 0.0

            if selection_mode == "support_first":
                # Amplicon assays: prioritize the amount of concordant read
                # support. Identity remains an important tie-breaker, but a
                # single perfect read cannot defeat hundreds of high-identity
                # reads mapping to the proper reference.
                rank = (matching_bases, mapped_reads, identity, mean_mapq, ref_name)
            else:
                # Whole-genome/default behavior: preserve identity-first
                # selection used by the HBV workflow.
                rank = (identity, matching_bases, mapped_reads, mean_mapq, ref_name)
            if best_rank is None or rank > best_rank:
                best_rank = rank
                best_ref = ref_name

        if not best_ref:
            return 0.0, 0, None, 0, 0

        aligned_bases = ref_aligned_query_bases.get(best_ref, 0)
        matching_bases = ref_matching_query_bases.get(best_ref, 0)
        identity_percent = round((matching_bases / aligned_bases) * 100, 6) if aligned_bases else 0.0

        # The caller compares this scalar across one-reference SAM files.
        if selection_mode == "support_first":
            # Matching bases dominate, followed by mapped-read count, with
            # identity as a fine tie-breaker.
            selection_score = (
                float(matching_bases) * 1_000_000.0
                + float(ref_mapped_reads.get(best_ref, 0)) * 1_000.0
                + identity_percent
            )
        else:
            # Preserve HBV/default identity-first behavior.
            selection_score = (identity_percent * 1_000_000_000.0) + float(matching_bases)

        return (
            selection_score,
            ref_mapped_reads.get(best_ref, 0),
            best_ref,
            aligned_bases,
            ref_lengths.get(best_ref, 0)
        )

    except Exception as e:
        print(f"Error evaluating alignment quality for {output_sam}: {e}")
        return 0.0, 0, None, 0, 0

def calculate_reference_breadth_coverage_percent(bam_file, ref_name, ref_length):
    if not bam_file or not os.path.exists(bam_file):
        return 0.0

    if not ref_name or not ref_length:
        return 0.0

    covered_positions = 0

    try:
        cmd = ["samtools", "depth", "-a", "-r", ref_name, bam_file]

        with subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        ) as proc:
            for line in proc.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue

                try:
                    depth = int(fields[2])
                    if depth > 0:
                        covered_positions += 1
                except ValueError:
                    continue

            _, stderr = proc.communicate()

            if proc.returncode != 0:
                print(f"samtools depth warning/error for {bam_file}: {stderr}")
                return 0.0

        return round((covered_positions / ref_length) * 100, 2) if ref_length > 0 else 0.0

    except Exception as e:
        print(f"Error calculating reference breadth coverage for {bam_file}: {e}")
        return 0.0

def empty_depth_stats():
    return {
        "mean_depth": 0.0,
        "median_depth": 0.0,
        "min_depth": 0,
        "positions_ge_10x_percent": 0.0,
        "positions_ge_30x_percent": 0.0,
    }

def calculate_depth_stats(bam_file, ref_name, region=None):
    depths = []

    if not bam_file or not os.path.exists(bam_file) or not ref_name:
        return empty_depth_stats()

    try:
        if region:
            region_string = f"{ref_name}:{region[0]}-{region[1]}"
            cmd = ["samtools", "depth", "-a", "-r", region_string, bam_file]
        else:
            cmd = ["samtools", "depth", "-a", "-r", ref_name, bam_file]

        with subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        ) as proc:
            for line in proc.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue

                try:
                    depths.append(int(fields[2]))
                except ValueError:
                    continue

            _, stderr = proc.communicate()

            if proc.returncode != 0:
                print(f"samtools depth warning/error for {bam_file}: {stderr}")

        if not depths:
            return empty_depth_stats()

        total = len(depths)

        return {
            "mean_depth": round(sum(depths) / total, 2),
            "median_depth": round(statistics.median(depths), 2),
            "min_depth": min(depths),
            "positions_ge_10x_percent": round((sum(1 for d in depths if d >= 10) / total) * 100, 2),
            "positions_ge_30x_percent": round((sum(1 for d in depths if d >= 30) / total) * 100, 2),
        }

    except Exception as e:
        print(f"Error calculating depth stats for {bam_file}: {e}")
        return empty_depth_stats()

def merge_depth_values(depth_values):
    if not depth_values:
        return empty_depth_stats()

    total = len(depth_values)

    return {
        "mean_depth": round(sum(depth_values) / total, 2),
        "median_depth": round(statistics.median(depth_values), 2),
        "min_depth": min(depth_values),
        "positions_ge_10x_percent": round((sum(1 for d in depth_values if d >= 10) / total) * 100, 2),
        "positions_ge_30x_percent": round((sum(1 for d in depth_values if d >= 30) / total) * 100, 2),
    }

def get_depth_values_for_region(bam_file, ref_name, region):
    depths = []

    if not bam_file or not os.path.exists(bam_file) or not ref_name or not region:
        return depths

    try:
        region_string = f"{ref_name}:{region[0]}-{region[1]}"
        cmd = ["samtools", "depth", "-a", "-r", region_string, bam_file]

        with subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        ) as proc:
            for line in proc.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                try:
                    depths.append(int(fields[2]))
                except ValueError:
                    continue

            _, stderr = proc.communicate()
            if proc.returncode != 0:
                print(f"samtools depth warning/error for {bam_file}: {stderr}")

    except Exception as e:
        print(f"Error getting depth values for region {region} in {bam_file}: {e}")

    return depths

def calculate_multi_region_depth_stats(bam_file, ref_name, regions):
    all_depths = []
    for region in regions:
        all_depths.extend(get_depth_values_for_region(bam_file, ref_name, region))
    return merge_depth_values(all_depths)

def detect_depth_supported_target_region(
    bam_file,
    ref_name,
    min_depth=100,
    max_gap=10,
):
    """Detect the main amplicon as the largest depth-supported block.

    Only reference positions with depth >= ``min_depth`` define the target.
    Short dips below the threshold are tolerated when consecutive supported
    positions are separated by at most ``max_gap`` bases. This prevents a few
    distant off-target reads from stretching an amplicon across most of an
    adenovirus genome.

    Returns ``(start_1based, end_1based, length)`` or ``(None, None, 0)``.
    """
    if not bam_file or not ref_name or min_depth < 1:
        return None, None, 0

    supported = []
    try:
        cmd = ["samtools", "depth", "-a", "-r", str(ref_name), str(bam_file)]
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)

        for line in result.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            try:
                pos = int(fields[1])
                depth = int(fields[2])
            except ValueError:
                continue
            if depth >= min_depth:
                supported.append((pos, depth))

        if not supported:
            return None, None, 0

        blocks = []
        current = [supported[0]]
        for item in supported[1:]:
            # Difference of 1 means directly adjacent. Allow up to max_gap
            # unsupported bases between qualifying positions.
            if item[0] - current[-1][0] <= max_gap + 1:
                current.append(item)
            else:
                blocks.append(current)
                current = [item]
        blocks.append(current)

        # Prefer the block with the greatest number of qualifying positions;
        # break ties using total depth and then genomic span.
        best = max(
            blocks,
            key=lambda block: (
                len(block),
                sum(depth for _, depth in block),
                block[-1][0] - block[0][0] + 1,
            ),
        )

        start = best[0][0]
        end = best[-1][0]
        return start, end, end - start + 1

    except Exception as e:
        print(f"Error detecting depth-supported target region for {bam_file}: {e}")
        return None, None, 0


def detect_aligned_target_region_from_sam(sam_file, ref_name):
    starts = []
    ends = []

    try:
        with open(sam_file, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("@"):
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue

                try:
                    flag = int(fields[1])
                    rname = fields[2]
                    pos = int(fields[3])
                    cigar = fields[5]
                except ValueError:
                    continue

                if (flag & 0x4) or rname != ref_name or cigar == "*":
                    continue

                ref_consumed = reference_consumed_by_cigar(cigar)
                if ref_consumed <= 0:
                    continue

                starts.append(pos)
                ends.append(pos + ref_consumed - 1)

        if not starts or not ends:
            return None, None, 0

        start = min(starts)
        end = max(ends)
        length = end - start + 1

        return start, end, length

    except Exception as e:
        print(f"Error detecting aligned target region from {sam_file}: {e}")
        return None, None, 0

def read_codon_from_aligned_segment(read, ref_positions_1based):
    pairs = read.get_aligned_pairs(matches_only=False)
    ref_to_query = {}

    for query_pos, ref_pos in pairs:
        if ref_pos is None:
            continue
        ref_to_query[ref_pos + 1] = query_pos

    bases = []
    seq = read.query_sequence or ""

    for ref_pos in ref_positions_1based:
        qpos = ref_to_query.get(ref_pos)
        if qpos is None:
            return None
        if qpos >= len(seq):
            return None
        base = seq[qpos].upper()
        if base not in "ACGT":
            return None
        bases.append(base)

    return "".join(bases)
