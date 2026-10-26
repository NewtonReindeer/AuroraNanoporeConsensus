"""Read-only CMV pre-consensus indel audit. Never edits BAM, FASTQ, or consensus."""
import csv
import json
from collections import defaultdict
from pathlib import Path


def classify_indel(length, homopolymer_length):
    if length <= 0:
        raise ValueError('Indel length must be positive')
    return ('in_frame' if length % 3 == 0 else 'frame_disrupting',
            'homopolymer_context' if homopolymer_length >= 4 else 'other_context')


def homopolymer_run(sequence, position):
    if position < 0 or position >= len(sequence):
        return 0
    base = sequence[position].upper()
    if base not in 'ACGT':
        return 0
    lo = position
    hi = position + 1
    while lo > 0 and sequence[lo - 1].upper() == base:
        lo -= 1
    while hi < len(sequence) and sequence[hi].upper() == base:
        hi += 1
    return hi - lo


def audit_alignment(bam_path, reference_fasta, reference_name, out_prefix, min_mapq=20):
    import pysam
    prefix = Path(out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    observations = defaultdict(lambda: {'read_ids': set(), 'forward': set(), 'reverse': set(), 'low_mapq': set()})
    stats = {'reads_seen': 0, 'primary_mapped_reads': 0, 'reads_with_indel': 0,
             'indel_events': 0, 'frame_disrupting_events': 0, 'in_frame_events': 0,
             'homopolymer_indel_events': 0, 'reference_name': reference_name,
             'reference_fasta': str(reference_fasta), 'bam': str(bam_path),
             'min_mapq': min_mapq, 'assessment': 'SCREENING_ONLY_NOT_VARIANT_CALLS'}
    with pysam.FastaFile(str(reference_fasta)) as fasta, pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        if reference_name not in bam.references or reference_name not in fasta.references:
            raise ValueError('Reference contig must be present in BAM and FASTA with the same name')
        if bam.get_reference_length(reference_name) != fasta.get_reference_length(reference_name):
            raise ValueError('BAM reference length differs from FASTA reference length')
        refseq = fasta.fetch(reference_name).upper()
        for read in bam.fetch(reference_name):
            stats['reads_seen'] += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail:
                continue
            stats['primary_mapped_reads'] += 1
            refpos = read.reference_start
            had_indel = False
            # CIGAR op codes: 0 M, 1 I, 2 D, 3 N, 4 S, 5 H, 6 P, 7 =, 8 X
            for op, length in read.cigartuples or ():
                if op in (1, 2):
                    had_indel = True
                    stats['indel_events'] += 1
                    anchor = max(0, refpos - 1) if op == 1 else refpos
                    hp = max(homopolymer_run(refseq, anchor),
                             homopolymer_run(refseq, min(anchor + 1, len(refseq) - 1)))
                    frame, context = classify_indel(length, hp)
                    stats[frame + '_events'] += 1
                    if context == 'homopolymer_context':
                        stats['homopolymer_indel_events'] += 1
                    key = (reference_name, anchor + 1, 'insertion' if op == 1 else 'deletion', length, frame, context)
                    item = observations[key]
                    # One molecule contributes once per identical event; avoid counting duplicated reads.
                    item['read_ids'].add(read.query_name)
                    item['reverse' if read.is_reverse else 'forward'].add(read.query_name)
                    if read.mapping_quality < min_mapq:
                        item['low_mapq'].add(read.query_name)
                if op in (0, 2, 3, 7, 8):
                    refpos += length
            if had_indel:
                stats['reads_with_indel'] += 1
    rows = []
    for key, item in sorted(observations.items()):
        row = dict(zip(('contig', 'anchor_position_1based', 'event', 'length', 'frame_class', 'context'), key))
        row.update({'supporting_distinct_read_names': len(item['read_ids']),
                    'forward_distinct_read_names': len(item['forward']),
                    'reverse_distinct_read_names': len(item['reverse']),
                    'low_mapq_distinct_read_names': len(item['low_mapq']),
                    'interpretation': 'REVIEW_ONLY_CIGAR_SIGNAL'})
        rows.append(row)
    stats['distinct_event_types'] = len(rows)
    json_path = str(prefix) + '.preconsensus_indel_qc.json'
    tsv_path = str(prefix) + '.preconsensus_indel_events.tsv'
    with open(json_path, 'w', encoding='utf-8') as out:
        json.dump(stats, out, indent=2)
    headers = ['contig', 'anchor_position_1based', 'event', 'length', 'frame_class', 'context',
               'supporting_distinct_read_names', 'forward_distinct_read_names',
               'reverse_distinct_read_names', 'low_mapq_distinct_read_names', 'interpretation']
    with open(tsv_path, 'w', newline='', encoding='utf-8') as out:
        writer = csv.DictWriter(out, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    return {**stats, 'json_path': json_path, 'tsv_path': tsv_path}


def write_consensus_comparison(draft_fasta, polished_fasta, out_prefix):
    """Metadata-only comparison. Length difference is NOT an indel/frameshift call."""
    from .fasta import fasta_sequence_stats
    draft = fasta_sequence_stats(draft_fasta)
    polished = fasta_sequence_stats(polished_fasta) if polished_fasta else None
    result = {
        'draft_length': draft['consensus_length'],
        'polished_length': polished['consensus_length'] if polished else None,
        'length_difference': polished['consensus_length'] - draft['consensus_length'] if polished else None,
        'warning': ('Global consensus length differences cannot localize indels, infer CDS frame, '
                    'or prove Medaka caused an error. A coordinate-resolved alignment is required.'),
        'status': 'COMPARISON_METADATA_ONLY' if polished else 'POLISHING_NOT_AVAILABLE',
    }
    path = str(out_prefix) + '.draft_medaka_summary.json'
    with open(path, 'w', encoding='utf-8') as out:
        json.dump(result, out, indent=2)
    return path
