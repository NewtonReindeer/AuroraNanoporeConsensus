"""Optional, CMV-only nucleotide VCF concordance. Research use only.

No Nextflow execution, no variant promotion, no interpretation of VCF as truth.
A viralrecon-generated VCF must use exactly the fixed Merlin coordinate system.
"""
import csv
import gzip
import json
from pathlib import Path

ACCESSION = 'NC_006273.2'


def _open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'rt')


def load_vcf(path, accession=ACCESSION):
    """Parse PASS SNP/indel alleles, rejecting incompatible contigs and dubious records."""
    variants = set()
    seen_header = False
    rejected = {'filtered': 0, 'symbolic_or_invalid': 0, 'other_contig': 0}
    with _open(path) as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                seen_header = True
                continue
            if line.startswith('#') or not line.strip():
                continue
            if not seen_header:
                raise ValueError('VCF header is missing')
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                raise ValueError('Invalid VCF record')
            contig, pos, _, ref, alts, _, filt, _ = fields[:8]
            if contig != accession:
                rejected['other_contig'] += 1
                continue
            if filt not in ('PASS', '.'):
                rejected['filtered'] += 1
                continue
            pos = int(pos)
            if pos <= 0:
                raise ValueError('Invalid VCF position')
            for alt in alts.split(','):
                if not ref or not alt or set(ref.upper()) - set('ACGT') or set(alt.upper()) - set('ACGT') or alt.upper() == ref.upper():
                    rejected['symbolic_or_invalid'] += 1
                else:
                    variants.add((pos, ref.upper(), alt.upper()))
    if not seen_header:
        raise ValueError('VCF header is missing')
    if rejected['other_contig'] and not variants:
        raise ValueError('No variants on %s: incompatible reference or empty target' % accession)
    return variants, rejected


def compare_to_aurora(vcf_path, aurora_calls, out_prefix, accession=ACCESSION):
    """Compare only reference-validated nucleotide alleles; retain unmatched as review only.

    Aurora nucleotide SNP rows often have 'major', not a validated 'ref': these
    MUST NOT be silently treated as reference-anchored VCF calls.
    """
    variants, rejected = load_vcf(vcf_path, accession)
    qualified = set()
    skipped = 0
    for call in aurora_calls:
        try:
            pos = int(call['position'])
            ref = str(call['ref']).upper()
            alt = str(call['alt']).upper()
        except (KeyError, ValueError, TypeError):
            skipped += 1
            continue
        if pos > 0 and ref != alt and ref and alt and not(set(ref + alt) - set('ACGT')):
            qualified.add((pos, ref, alt))
        else:
            skipped += 1
    entries = []
    for pos, ref, alt in sorted(variants | qualified):
        key = (pos, ref, alt)
        entries.append({'chrom': accession, 'position': pos, 'ref': ref, 'alt': alt,
                        'viralrecon': key in variants, 'aurora': key in qualified,
                        'status': 'MATCH' if key in variants and key in qualified else ('VIRALRECON_ONLY' if key in variants else 'AURORA_ONLY')})
    output = Path(out_prefix)
    output.parent.mkdir(parents=True, exist_ok=True)
    tsv = str(output) + '.tsv'
    with open(tsv, 'w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=['chrom', 'position', 'ref', 'alt', 'viralrecon', 'aurora', 'status'], delimiter='\t')
        writer.writeheader()
        writer.writerows(entries)
    summary = {'status': 'REVIEW_ONLY', 'source_vcf': str(vcf_path), 'reference': accession,
               'viralrecon_alleles': len(variants), 'aurora_reference_anchored_alleles': len(qualified),
               'aurora_unverifiable_alleles': skipped,
               'matches': sum(x['status'] == 'MATCH' for x in entries),
               'viralrecon_only': sum(x['status'] == 'VIRALRECON_ONLY' for x in entries),
               'aurora_only': sum(x['status'] == 'AURORA_ONLY' for x in entries),
               'vcf_rejections': rejected,
               'warning': 'VCF normalization, reference identity, primer scheme and platform suitability require independent verification. No clinical conclusions.'}
    js = str(output) + '.json'
    with open(js, 'w') as fp:
        json.dump(summary, fp, indent=2)
    return summary, js, tsv


def safe_optional_review(sample, aurora_result, sample_dir):
    """Only opt-in; local VCF must be explicitly supplied by a CMV operator."""
    import os
    directory = os.environ.get('AURORA_CMV_VIRALRECON_VCF_DIR', '').strip()
    if not directory:
        return aurora_result
    # Avoid sample IDs being interpreted as paths.
    if Path(sample).name != sample or '/' in sample or '\\' in sample:
        raise ValueError('Invalid CMV sample identifier')
    directory = Path(directory)
    matches = [directory / (sample + '.vcf'), directory / (sample + '.vcf.gz')]
    present = [p for p in matches if p.is_file()]
    if len(present) != 1:
        aurora_result['cmv_viralrecon_status'] = 'NOT_ASSESSED: expected one sample-named VCF in configured directory'
        return aurora_result
    try:
        # The current Aurora SNP rows report a major allele, not REF.
        # Do not convert major into REF; only compare truly reference-anchored entries.
        alleles = [x for x in aurora_result.get('cmv_variant_calls', []) + aurora_result.get('cmv_indel_calls', []) if 'ref' in x]
        summary, js, tsv = compare_to_aurora(present[0], alleles, Path(sample_dir) / (sample + '.cmv_viralrecon_review'))
        aurora_result['cmv_viralrecon_status'] = 'REVIEW_ONLY: %s matched alleles; see %s' % (summary['matches'], js)
        aurora_result['cmv_viralrecon_review_json'] = js
        aurora_result['cmv_viralrecon_review_tsv'] = tsv
    except Exception as exc:
        aurora_result['cmv_viralrecon_status'] = 'REVIEW_FAILED: %s' % exc
    return aurora_result
