"""Standalone, read-only CMV UL54 codon-coordinate and BAM evidence validator.

Requires pysam for real BAM/FASTA; does not change Aurora variant calls.
Run: python cmv_locus_validation.py --bam FILE.bam --fasta refs.fa
     --accession NC_006273.2 --coordinates cmv_gene_coordinates.tsv
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

COMPLEMENT = str.maketrans('ACGT', 'TGCA')
CODON_TO_AA = {}
# Standard genetic code, in T/C/A/G order per base.
_aa = ('FFLLSSSSYY**CC*W' 'LLLLPPPPHHQQRRRR' 'IIIMTTTTNNKKSSRR' 'VVVVAAAADDEEGGGG')
for i, a in enumerate('TCAG'):
    for j, b in enumerate('TCAG'):
        for k, c in enumerate('TCAG'):
            CODON_TO_AA[a+b+c] = _aa[i*16+j*4+k]


def codon_positions(entry, number):
    """Return ordered, zero-based reference bases in transcription orientation."""
    start, end = int(entry['start']), int(entry['end'])
    first = int(entry['first_codon_number'])
    strand = entry['strand']
    if start < 1 or end < start or strand not in ('+', '-') or (end-start+1) % 3:
        raise ValueError('Invalid CDS interval/strand: ' + str(entry))
    index = number - first
    if index < 0 or (index+1)*3 > end-start+1:
        raise ValueError(f'Codon {number} is outside the declared CDS')
    if strand == '+':
        return tuple(start-1+index*3+i for i in range(3))
    return tuple(end-1-index*3-i for i in range(3))


def gene_entry(coordinates_path, accession, gene='UL54'):
    with open(coordinates_path, newline='', encoding='utf-8-sig') as f:
        rows = [r for r in csv.DictReader(f, delimiter='\t')
                if r.get('accession') == accession and r.get('gene') == gene]
    if len(rows) != 1:
        raise ValueError(f'Expected exactly one coordinate entry for {accession}/{gene}, got {len(rows)}')
    return rows[0]


def reference_codon(fasta, accession, entry, number):
    positions = codon_positions(entry, number)
    seq = ''.join(fasta.fetch(accession, p, p+1).upper() for p in positions)
    if len(seq) != 3 or any(c not in 'ACGT' for c in seq):
        raise ValueError(f'Incomplete or ambiguous reference sequence at codon {number}: {seq}')
    if entry['strand'] == '-':
        seq = seq.translate(COMPLEMENT)
    return seq


def classify_read(read, positions, strand, min_baseq=20, local_window=12):
    """Individually classify codon evidence; no sample-wide frame assumption.

    Non-primary/duplicate/low-MAPQ filtering happens in audit_loci().
    A covered codon adjacent to an indel is retained as diagnostic-only.
    """
    bases = {}
    for q, rp in read.get_aligned_pairs(matches_only=True):
        if q is not None and rp in positions:
            bases[rp] = (read.query_sequence[q].upper(), read.query_qualities[q])
    if not all(p in bases for p in positions):
        return 'incomplete_or_deleted', None
    if any(bases[p][0] not in 'ACGT' for p in positions):
        return 'ambiguous_base', None
    if min(bases[p][1] for p in positions) < min_baseq:
        return 'low_base_quality', None
    codon = ''.join(bases[p][0] for p in positions)
    if strand == '-':
        codon = codon.translate(COMPLEMENT)
    from math import fmod
    lo, hi = min(positions), max(positions)
    cursor = read.reference_start
    events = []
    for op, length in (read.cigartuples or []):
        if op == 1:
            events.append((cursor-1, length, 'insertion'))
        elif op == 2:
            events.append((cursor, -length, 'deletion'))
        if op in (0, 2, 3, 7, 8):
            cursor += length
    if any(lo-local_window <= anchor <= hi+local_window for anchor, _, _ in events):
        return 'local_indel', codon
    upstream = [(p, delta) for p, delta, _ in events
                if (p < lo if strand == '+' else p > hi)]
    if sum(delta for _, delta in upstream) % 3:
        return 'unresolved_frame', codon
    if upstream:
        return 'restored_or_inframe', codon
    return 'clean', codon


def audit_loci(bam_path, fasta_path, coordinates_path, accession, loci=(522, 552),
               min_mapq=20, min_baseq=20, local_window=12):
    import pysam
    entry = gene_entry(coordinates_path, accession)
    output = {'bam': str(bam_path), 'fasta': str(fasta_path),
              'coordinates': str(coordinates_path), 'accession': accession,
              'gene': 'UL54', 'strand': entry['strand'], 'first_codon_number': int(entry['first_codon_number']),
              'reference_interval_1based': [int(entry['start']), int(entry['end'])], 'loci': []}
    with pysam.FastaFile(str(fasta_path)) as fasta, pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        if accession not in fasta.references or accession not in bam.references:
            raise ValueError(f'Accession {accession} missing from FASTA or BAM; verify the actual mapping reference')
        for number in loci:
            pos = codon_positions(entry, number)
            ref = reference_codon(fasta, accession, entry, number)
            refaa = CODON_TO_AA[ref]
            status = 'REF_MATCHES_P' if refaa == 'P' else 'REF_IS_NOT_P'
            counters = Counter()
            evidence = Counter()
            records = []
            for read in bam.fetch(accession, min(pos), max(pos)+1):
                if (read.is_unmapped or read.is_secondary or read.is_supplementary or
                    read.is_duplicate or read.is_qcfail or read.mapping_quality < min_mapq or
                    read.query_sequence is None or read.query_qualities is None):
                    counters['excluded_alignment'] += 1
                    continue
                if not (read.reference_start <= min(pos) and max(pos) < read.reference_end):
                    continue
                category, observed = classify_read(read, pos, entry['strand'], min_baseq, local_window)
                counters[category] += 1
                if observed:
                    alt = CODON_TO_AA.get(observed, '?')
                    evidence[(category, observed, alt)] += 1
                    records.append({'read': read.query_name, 'category': category, 'codon': observed,
                                    'aa': alt, 'mapq': read.mapping_quality,
                                    'cigar': read.cigarstring, 'reverse_alignment': bool(read.is_reverse)})
            clean_depth = sum(n for (cat, _, _), n in evidence.items() if cat == 'clean')
            restored_depth = sum(n for (cat, _, _), n in evidence.items() if cat == 'restored_or_inframe')
            clean_A = sum(n for (cat, _, aa), n in evidence.items() if cat == 'clean' and aa == 'A')
            restored_A = sum(n for (cat, _, aa), n in evidence.items()
                             if cat == 'restored_or_inframe' and aa == 'A')
            output['loci'].append({'codon_number': number,
                'reference_positions_1based_transcription_order': [p+1 for p in pos],
                'reference_codon': ref, 'reference_amino_acid': refaa,
                'proline_variant_coordinate_status': status,
                'read_category_counts': dict(sorted(counters.items())),
                'codon_evidence': [dict(category=c, codon=co, amino_acid=a, reads=n)
                                   for (c, co, a), n in sorted(evidence.items())],
                'clean_depth': clean_depth, 'clean_alanine_reads': clean_A,
                'reconstructed_review_depth': clean_depth+restored_depth,
                'reconstructed_review_alanine_reads': clean_A+restored_A,
                'read_details': records})
    return output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bam', required=True)
    p.add_argument('--fasta', required=True)
    p.add_argument('--coordinates', required=True)
    p.add_argument('--accession', default='NC_006273.2')
    p.add_argument('--loci', nargs='+', type=int, default=[522, 552])
    p.add_argument('--out-prefix', default='CMV_UL54_P522_P552_validation')
    args = p.parse_args()
    result = audit_loci(args.bam, args.fasta, args.coordinates, args.accession, args.loci)
    target = Path(args.out_prefix)
    target.parent.mkdir(parents=True, exist_ok=True)
    summary = {**result, 'loci': [{k:v for k,v in locus.items() if k != 'read_details'} for locus in result['loci']]}
    target.with_suffix('.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    with target.with_suffix('.tsv').open('w', encoding='utf-8', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['codon_number','read','category','codon','aa','mapq','cigar','reverse_alignment'], delimiter='\t')
        w.writeheader()
        for locus in result['loci']:
            for read in locus['read_details']:
                w.writerow({'codon_number': locus['codon_number'], **read})
    print(json.dumps(summary, indent=2))
    print('Wrote:', target.with_suffix('.json'), target.with_suffix('.tsv'))


if __name__ == '__main__':
    main()
