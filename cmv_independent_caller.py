"""Optional CMV-only independent minimap2/Clair3/bcftools evidence branch.

Research-use review; never modifies Aurora's primary mutation calls.
No ARTIC, Nextflow, automatic resistance grading, or primer trimming.
"""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

ACCESSION = 'NC_006273.2'


class ReviewError(RuntimeError):
    pass


def _run(argv, stdout_path=None):
    if stdout_path:
        with open(stdout_path, 'wb') as out:
            subprocess.run(argv, check=True, stdout=out)
    else:
        subprocess.run(argv, check=True)


def _exists(path, label):
    p = Path(path)
    if not p.is_file():
        raise ReviewError(f'{label} file not found: {p}')
    return p.resolve()


def _tool(name):
    # Prefer BCFtools from the active Python/Conda environment.
    if name == "bcftools":
        import sys
        candidate = Path(sys.prefix) / "bin" / "bcftools"
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)

    path = shutil.which(name)
    if not path:
        raise ReviewError(f'Required executable missing: {name}')
    return path


def _contigs(fasta):
    result = {}
    name, seq = None, []
    with open(fasta) as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    if name in result:
                        raise ReviewError('Duplicate FASTA contig')
                    result[name] = ''.join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        if name in result:
            raise ReviewError('Duplicate FASTA contig')
        result[name] = ''.join(seq).upper()
    return result


def validate_reference(fasta, accession=ACCESSION):
    sequences = _contigs(fasta)
    if len(sequences) != 1 or accession not in sequences:
        raise ReviewError(f'Expected single-contig FASTA named {accession}, found {list(sequences)}')
    seq = sequences[accession]
    if not seq or set(seq) - set('ACGTN'):
        raise ReviewError('Reference contains invalid or empty sequence')
    return {'contig': accession, 'length': len(seq), 'sha256': hashlib.sha256(seq.encode()).hexdigest()}


def validate_vcf(vcf, fasta, accession=ACCESSION):
    """Reject bad REF, non-normalized output and incompatible CHROM (including empty VCF)."""
    seq = _contigs(fasta)[accession]
    opener = gzip.open if str(vcf).endswith('.gz') else open
    total = 0
    header = False
    with opener(vcf, 'rt') as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                header = True
                continue
            if line.startswith('#') or not line.strip():
                continue
            if not header:
                raise ReviewError('VCF lacks #CHROM header')
            fields = line.rstrip().split('\t')
            if len(fields) < 8:
                raise ReviewError('Malformed VCF row')
            chrom, pos, _, ref, alts = fields[:5]
            if chrom != accession:
                raise ReviewError(f'Unexpected VCF contig {chrom}; expected {accession}')
            start = int(pos) - 1
            if start < 0 or start + len(ref) > len(seq) or seq[start:start + len(ref)] != ref.upper():
                raise ReviewError(f'VCF reference mismatch at {chrom}:{pos} {ref}')
            for alt in alts.split(','):
                if not alt or set(alt.upper()) - set('ACGTN'):
                    raise ReviewError(f'Unsupported allele at {chrom}:{pos}: {alt}')
            total += 1
    if not header:
        raise ReviewError('VCF lacks #CHROM header')
    return total


def compare_vcfs(primary_vcf, independent_vcf, out_prefix, fasta):
    """Compare verified allele keys; missing callers/empty output are not negative diagnoses."""
    from .cmv_viralrecon_bridge import load_vcf
    validate_vcf(independent_vcf, fasta)
    external, _ = load_vcf(independent_vcf)
    if primary_vcf:
        validate_vcf(primary_vcf, fasta)
        primary, _ = load_vcf(primary_vcf)
    else:
        primary = set()
    rows = []
    for pos, ref, alt in sorted(external | primary):
        rows.append({'chrom': ACCESSION, 'pos': pos, 'ref': ref, 'alt': alt,
                     'independent': int((pos, ref, alt) in external),
                     'aurora_vcf': int((pos, ref, alt) in primary),
                     'status': ('CONCORDANT' if (pos, ref, alt) in external and (pos, ref, alt) in primary
                                else 'INDEPENDENT_ONLY' if (pos, ref, alt) in external else 'AURORA_ONLY') if primary_vcf else 'INDEPENDENT_ONLY_UNCOMPARABLE'})
    out = Path(out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(str(out) + '.tsv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['chrom','pos','ref','alt','independent','aurora_vcf','status'], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    summary = {'status': 'REVIEW_ONLY', 'independent_alleles': len(external),
               'aurora_alleles': len(primary) if primary_vcf else None,
               'concordant': sum(row['status'] == 'CONCORDANT' for row in rows),
               'primary_vcf_available': bool(primary_vcf),
               'warning': 'Not clinically validated. Caller agreement is not ground truth. Low-frequency mixtures may be missed.'}
    with open(str(out) + '.json', 'w') as handle:
        json.dump(summary, handle, indent=2)
    return summary


def run(sample, fastq, merlin_fasta, model_dir, sample_dir, *, threads=4, bed=None, primary_vcf=None):
    """Explicitly opt-in; isolated output folder. Requires local Clair3 runner + tools."""
    if Path(sample).name != sample or '/' in sample or '\\' in sample:
        raise ReviewError('Invalid sample ID')
    if not (1 <= threads <= 128):
        raise ReviewError('threads must be in 1..128')
    fastq = _exists(fastq, 'Original FASTQ')
    merlin_fasta = _exists(merlin_fasta, 'Merlin FASTA')
    model_dir = Path(model_dir).resolve()
    if not model_dir.is_dir():
        raise ReviewError(f'Clair3 model directory not found: {model_dir}')
    reference = validate_reference(merlin_fasta)
    if bed:
        bed = _exists(bed, 'Custom-target BED')
    minimap2, samtools, bcftools = (_tool(x) for x in ('minimap2', 'samtools', 'bcftools'))
    clair3 = _tool('run_clair3.sh')
    root = Path(sample_dir).resolve() / 'cmv_independent_review'
    root.mkdir(parents=True, exist_ok=True)
    bam = root / (sample + '.merlin.independent.sorted.bam')
    # Each subprocess has checked return codes; no shell command interpolation.
    p1 = subprocess.Popen([minimap2, '-ax', 'map-ont', '-t', str(threads), str(merlin_fasta), str(fastq)], stdout=subprocess.PIPE)
    try:
        subprocess.run([samtools, 'sort', '-@', str(threads), '-o', str(bam), '-'], stdin=p1.stdout, check=True)
    finally:
        p1.stdout.close()
    if p1.wait() != 0:
        raise ReviewError('minimap2 alignment failed')
    _run([samtools, 'index', str(bam)])
    _run([samtools, 'faidx', str(merlin_fasta)])
    clair_dir = root / 'clair3'
    cmd = [clair3, '--bam_fn=' + str(bam), '--ref_fn=' + str(merlin_fasta),
           '--threads=' + str(threads), '--platform=ont', '--model_path=' + str(model_dir),
           '--output=' + str(clair_dir), '--include_all_ctgs', '--haploid_precise',
           '--no_phasing_for_fa', '--enable_variant_calling_at_sequence_head_and_tail']
    if bed:
        cmd.append('--bed_fn=' + str(bed))
    _run(cmd)
    raw = _exists(clair_dir / 'merge_output.vcf.gz', 'Clair3 output')
    normalized = root / (sample + '.merlin.norm.vcf.gz')
    _run([bcftools, 'norm', '-f', str(merlin_fasta), '-m', '-any', '-Oz', '-o', str(normalized), str(raw)])
    _run([bcftools, 'index', '-t', str(normalized)])
    validate_vcf(normalized, merlin_fasta)
    if primary_vcf:
        primary_vcf = _exists(primary_vcf, 'Aurora normalized primary VCF')
    summary = compare_vcfs(primary_vcf, normalized, root / (sample + '.concordance'), merlin_fasta)
    summary.update({'sample':sample, 'reference':reference, 'fastq':str(fastq),
                    'model_path':str(model_dir),'bam':str(bam),'normalized_vcf':str(normalized),
                    'target_bed':str(bed) if bed else None})
    with open(root / (sample + '.manifest.json'), 'w') as handle:
        json.dump(summary, handle, indent=2)
    return summary


def optional_run(sample, fastq, reference_dir, sample_dir, current_result):
    """Never enter branch unless explicit enable flag is set; never change variant calls."""
    if os.environ.get('AURORA_CMV_INDEPENDENT_CALLER', '').strip() != '1':
        return current_result
    result = dict(current_result)
    try:
        fasta = os.environ.get('AURORA_CMV_MERLIN_FASTA', '').strip()
        model = os.environ.get('AURORA_CMV_CLAIR3_MODEL', '').strip()
        if not fasta or not model:
            raise ReviewError('Set AURORA_CMV_MERLIN_FASTA and AURORA_CMV_CLAIR3_MODEL')
        threads = int(os.environ.get('AURORA_CMV_CLAIR3_THREADS', '4'))
        bed = os.environ.get('AURORA_CMV_TARGET_BED') or None
        primary = os.environ.get('AURORA_CMV_PRIMARY_VCF') or None
        summary = run(sample, fastq, fasta, model, sample_dir, threads=threads, bed=bed, primary_vcf=primary)
        result['cmv_independent_status'] = 'REVIEW_ONLY'
        result['cmv_independent_vcf'] = summary['normalized_vcf']
        result['cmv_independent_alleles'] = summary['independent_alleles']
    except Exception as exc:
        result['cmv_independent_status'] = 'REVIEW_FAILED: ' + str(exc)
    return result
