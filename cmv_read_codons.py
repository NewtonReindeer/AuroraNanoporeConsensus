"""Independent CMV read-level codon audit (NOT a clinical caller).

Reference-strand-aware CIGAR interpretation; keep separate clean, restored-frame,
local-indel and unresolved-frame evidence. Only clean reads form primary audit calls.
"""
import csv
from collections import Counter
from pathlib import Path
from .cmv_codon import _CODONS, _COMPLEMENT
from .cmv_thresholds import MIN_ALT_FRACTION, MIN_CLEAN_DEPTH, MIN_ALT_READS


def indel_events(cigartuples, reference_start):
    """Yield (zero-based reference anchor, signed query-minus-reference length)."""
    refpos = reference_start
    for op, length in cigartuples or ():
        if op == 1:  # insertion after preceding reference position
            yield (refpos - 1, length)
        elif op == 2:  # deletion begins at refpos
            yield (refpos, -length)
        if op in (0, 2, 3, 7, 8):
            refpos += length


def frame_context(events, positions, strand, local_window=12):
    """Classify read's CIGAR around a codon; all positions are zero-based.

    A restored frame is still a REVIEW group: intervening amino acids may differ.
    Events outside the sequenced read cannot be inferred.
    """
    lo, hi = min(positions), max(positions)
    local = any(lo-local_window <= anchor <= hi+local_window for anchor, _ in events)
    upstream = [delta for anchor, delta in events
                if (anchor < lo if strand == '+' else anchor > hi)]
    if local:
        return 'local_indel'
    if sum(upstream) % 3:
        return 'unresolved_upstream_frame'
    if upstream:
        return 'restored_or_inframe_upstream'
    return 'clean'


def audit_read_codons(bam_path, reference_path, reference_name, min_mapq=20,
                      min_baseq=20, min_depth=MIN_CLEAN_DEPTH, min_alt_reads=MIN_ALT_READS,
                      min_fraction=MIN_ALT_FRACTION, local_window=12, include_qc=False):
    import pysam
    coord_file = Path(reference_path).parent / 'cmv_gene_coordinates.tsv'
    if not coord_file.is_file():
        return ([], 'Missing gene coordinates: ' + str(coord_file), []) if include_qc else ([], 'Missing gene coordinates: ' + str(coord_file))
    with coord_file.open(encoding='utf-8-sig', newline='') as fh:
        entries = [r for r in csv.DictReader(fh, delimiter='\t')
                   if r.get('accession') == reference_name and r.get('gene') in ('UL54', 'UL97')]
    if not entries:
        return ([], 'No UL54/UL97 coordinate entries for ' + reference_name, []) if include_qc else ([], 'No UL54/UL97 coordinate entries for ' + reference_name)
    results = []
    codon_qc = []
    with pysam.AlignmentFile(bam_path, 'rb') as bam, pysam.FastaFile(reference_path) as fasta:
        for entry in entries:
            start, end = int(entry['start']), int(entry['end'])
            strand, offset = entry['strand'], int(entry['first_codon_number'])
            if start < 1 or end < start or strand not in ('+', '-') or (end-start+1) % 3 or offset < 1:
                return ([], 'Invalid CDS geometry for ' + entry['gene'], []) if include_qc else ([], 'Invalid CDS geometry for ' + entry['gene'])
            seq = fasta.fetch(reference_name, start-1, end).upper()
            if len(seq) != end-start+1 or any(base not in 'ACGT' for base in seq):
                return ([], 'Invalid or incomplete reference CDS for ' + entry['gene'], []) if include_qc else ([], 'Invalid or incomplete reference CDS for ' + entry['gene'])
            if strand == '-':
                seq = seq.translate(_COMPLEMENT)[::-1]
            n = len(seq)//3
            groups = [dict() for _ in range(n)]
            # Track reconstructed frame state independently of the conservative caller.
            # Nonzero-frame reads remain excluded; restored/in-frame reads are review-only.
            excluded = [Counter() for _ in range(n)]
            for read in bam.fetch(reference_name, start-1, end):
                if (read.is_unmapped or read.is_secondary or read.is_supplementary or
                    read.is_duplicate or read.is_qcfail or read.mapping_quality < min_mapq or
                    read.query_sequence is None or read.query_qualities is None):
                    continue
                bases = {}
                for q, rp in read.get_aligned_pairs(matches_only=True):
                    if rp is not None and q is not None and start-1 <= rp < end:
                        bases[rp] = (read.query_sequence[q].upper(), read.query_qualities[q])
                events = [(anchor, delta) for anchor, delta in indel_events(read.cigartuples, read.reference_start)
                          if start-1 <= anchor < end]
                # Iterate only codons potentially overlapping this alignment.
                for j in range(n):
                    positions = ([start-1 + j*3+k for k in range(3)] if strand == '+'
                                 else [end-1-j*3-k for k in range(3)])
                    if not (read.reference_start <= min(positions) and
                            max(positions) < read.reference_end):
                        continue
                    if not all(p in bases and bases[p][0] in 'ACGT' for p in positions):
                        excluded[j]['incomplete_or_deleted'] += 1
                        continue
                    if min(bases[p][1] for p in positions) < min_baseq:
                        excluded[j]['low_base_quality'] += 1
                        continue
                    observed = ''.join(bases[p][0] for p in positions)
                    if strand == '-':
                        observed = observed.translate(_COMPLEMENT)
                    category = frame_context(events, positions, strand, local_window)
                    groups[j].setdefault(category, Counter())[observed] += 1
            for j, by_group in enumerate(groups):
                clean = by_group.get('clean', Counter())
                depth = sum(clean.values())
                ref_codon = seq[j*3:j*3+3]
                ref_aa = _CODONS.get(ref_codon)
                if not ref_aa or ref_aa == '*':
                    codon_qc.append(dict(gene=entry['gene'], codon=offset+j,
                        reference_codon=ref_codon, assessment='INVALID_REFERENCE_CODON',
                        clean_depth=0, clean_reference_reads=0, local_indel_reads=0,
                        restored_frame_reads=0, unresolved_frame_reads=0,
                        low_baseq_reads=0, incomplete_reads=0))
                    continue
                # Independently reconstructed, fully aligned codons with net-zero frame shift.
                # Report them as supplemental evidence, not passing calls: an in-frame
                # insertion/deletion can still change the intervening protein sequence.
                reconstructed = clean + by_group.get('restored_or_inframe_upstream', Counter())
                reconstructed_depth = sum(reconstructed.values())
                excluded_frame = sum(by_group.get('unresolved_upstream_frame', {}).values())
                local_depth = sum(by_group.get('local_indel', {}).values())
                restored_depth = sum(by_group.get('restored_or_inframe_upstream', {}).values())
                missing = excluded[j]['incomplete_or_deleted']
                low_quality = excluded[j]['low_base_quality']
                if depth >= min_depth:
                    qc_status = 'ASSESSED_CLEAN_NO_ALT' if not any(
                        _CODONS.get(c) != ref_aa for c in clean if _CODONS.get(c)
                    ) else 'ASSESSED_CLEAN_WITH_ALT'
                elif depth == 0 and not (excluded_frame or local_depth or restored_depth or missing or low_quality):
                    qc_status = 'NO_ELIGIBLE_EVIDENCE'
                else:
                    qc_status = 'NOT_ASSESSABLE_LOW_CLEAN_DEPTH'
                codon_qc.append(dict(gene=entry['gene'], codon=offset+j,
                    reference_codon=ref_codon, clean_depth=depth,
                    reconstructed_depth=reconstructed_depth,
                    reconstructed_reference_reads=sum(count for c, count in reconstructed.items() if _CODONS.get(c) == ref_aa),
                    clean_reference_reads=sum(count for c, count in clean.items() if _CODONS.get(c) == ref_aa),
                    local_indel_reads=local_depth, restored_frame_reads=restored_depth,
                    unresolved_frame_reads=excluded_frame,
                    low_baseq_reads=low_quality, incomplete_reads=missing,
                    assessment=qc_status))
                # Retain audit rows at sites with alternate evidence in any category,
                # even if clean depth is below the calling threshold.
                alternate_codons = set().union(*(set(c) for c in by_group.values())) if by_group else set()
                for codon in sorted(alternate_codons):
                    aa = _CODONS.get(codon)
                    if not aa or aa == ref_aa:
                        continue
                    support = clean[codon]
                    passes = depth >= min_depth and support >= min_alt_reads and support/depth >= min_fraction
                    reconstructed_support = reconstructed[codon]
                    reconstructed_fraction = reconstructed_support/reconstructed_depth if reconstructed_depth else 0.0
                    reconstructed_pass = (reconstructed_depth >= min_depth and
                                          reconstructed_support >= min_alt_reads and
                                          reconstructed_fraction >= min_fraction)
                    results.append(dict(
                        gene=entry['gene'], mutation=f'{ref_aa}{offset+j}{aa}',
                        reference_codon=ref_codon, alternate_codon=codon,
                        depth=depth, alt_reads=support,
                        reconstructed_depth=reconstructed_depth,
                        reconstructed_alt_reads=reconstructed_support,
                        reconstructed_frequency=round(100*reconstructed_fraction, 2),
                        reconstructed_threshold_pass=reconstructed_pass,
                        frequency=round(100*support/depth, 2) if depth else 0,
                        frameshift_reads_excluded=sum(by_group.get('unresolved_upstream_frame', {}).values()),
                        local_indel_support=by_group.get('local_indel', Counter())[codon],
                        restored_frame_support=by_group.get('restored_or_inframe_upstream', Counter())[codon],
                        unresolved_frame_support=by_group.get('unresolved_upstream_frame', Counter())[codon],
                        low_baseq_reads=excluded[j]['low_base_quality'],
                        incomplete_reads=excluded[j]['incomplete_or_deleted'],
                        status='READ_PHASED_EXPLORATORY' if passes else 'READ_LEVEL_REVIEW_INSUFFICIENT_CLEAN_SUPPORT'))
    note = ('Read-level codon audit: conservative clean-read calls unchanged; '
            'net-frame-reconstructed evidence is review-only, not a clinical call')
    return (results, note, codon_qc) if include_qc else (results, note)
