"""CMV exploratory SNV analysis. Clinical interpretation requires a validated panel."""
import csv
import os
from collections import defaultdict

CMV_GENES = ('UL97', 'UL54', 'UL56', 'UL89', 'UL51')
from .cmv_thresholds import MIN_ALT_FRACTION, MIN_CLEAN_DEPTH, MIN_ALT_READS
MIN_DEPTH = MIN_CLEAN_DEPTH
MIN_INDEL_READS = 3


def load_cmv_panel(reference_path):
    """Optional coordinate-specific panel beside the reference FASTA.

    TSV columns: position, ref, alt, gene, mutation, interpretation, evidence.
    Positions are 1-based coordinates on the selected reference sequence.
    """
    path = os.path.join(os.path.dirname(reference_path), 'cmv_mutation_panel.tsv')
    if not os.path.isfile(path):
        return {}, 'No reference-specific CMV mutation panel supplied; resistance interpretation unavailable.'
    panel = {}
    with open(path, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'position', 'ref', 'alt', 'gene', 'mutation', 'interpretation', 'evidence'}
        if not required.issubset(set(reader.fieldnames or [])):
            return {}, 'CMV mutation panel missing required columns; interpretation unavailable.'
        for item in reader:
            try:
                pos = int(item['position'])
                ref, alt = item['ref'].upper(), item['alt'].upper()
                if pos < 1 or len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
                    continue
                panel[(pos, ref, alt)] = item
            except (ValueError, TypeError):
                continue
    return panel, 'Reference-coordinate panel loaded; verify reference identity and clinical validation.'


def analyse_cmv_mutations(sample, bam_file, ref_name, ref_path):
    import pysam
    panel, note = load_cmv_panel(ref_path)
    calls = []
    assessed = 0
    site_counts = {}
    covered_positions = set()
    indels = []
    with pysam.FastaFile(ref_path) as fasta, pysam.AlignmentFile(bam_file, 'rb') as bam:
        if ref_name not in fasta.references or ref_name not in bam.references:
            raise ValueError('BAM/FASTA reference accession mismatch: ' + ref_name)
        if fasta.get_reference_length(ref_name) != bam.get_reference_length(ref_name):
            raise ValueError('BAM/FASTA reference length mismatch: ' + ref_name)
        for column in bam.pileup(ref_name, truncate=True, stepper='samtools', min_base_quality=20, min_mapping_quality=20, max_depth=100000):
            counts = {base: 0 for base in 'ACGT'}
            # Indels are anchored to the preceding reference base (VCF convention).
            # A single read contributes at most one event at this anchor.
            events = {}
            informative = 0
            for read in column.pileups:
                if read.is_del or read.is_refskip or read.query_position is None:
                    continue
                aln = read.alignment
                if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.is_duplicate or aln.is_qcfail:
                    continue
                informative += 1
                if read.indel != 0:
                    anchor = column.reference_pos
                    try:
                        anchor_base = fasta.fetch(ref_name, anchor, anchor + 1).upper()
                        if len(anchor_base) != 1 or anchor_base not in 'ACGT':
                            raise ValueError('invalid anchor')
                        if read.indel < 0:
                            length = -read.indel
                            deleted = fasta.fetch(ref_name, anchor + 1, anchor + 1 + length).upper()
                            if len(deleted) != length or any(c not in 'ACGT' for c in deleted):
                                raise ValueError('invalid deletion')
                            key = ('DEL', anchor_base + deleted, anchor_base, length)
                        else:
                            length = read.indel
                            start = read.query_position + 1
                            inserted = aln.query_sequence[start:start + length].upper()
                            if len(inserted) != length or any(c not in 'ACGT' for c in inserted):
                                raise ValueError('invalid insertion')
                            key = ('INS', anchor_base, anchor_base + inserted, length)
                        events[key] = events.get(key, 0) + 1
                    except (ValueError, KeyError, IndexError):
                        pass
                base = aln.query_sequence[read.query_position].upper()
                if base in counts:
                    counts[base] += 1
            # The denominator includes reads spanning the anchor with no indel.
            if informative >= MIN_DEPTH:
                for (kind, ref, alt, length), support in sorted(events.items()):
                    fraction = support / informative
                    if support < MIN_INDEL_READS or fraction < MIN_ALT_FRACTION:
                        continue
                    indels.append({'position': column.reference_pos + 1,
                                   'type': 'Deletion' if kind == 'DEL' else 'Insertion',
                                   'ref': ref, 'alt': alt, 'length': length,
                                   'depth': informative, 'alt_reads': support,
                                   'frequency': round(fraction * 100, 2),
                                   'interpretation': 'Unclassified; protein consequence not assessed'})
            depth = sum(counts.values())
            if depth < MIN_DEPTH:
                continue
            covered_positions.add(column.reference_pos + 1)
            site_counts[column.reference_pos + 1] = (counts.copy(), depth)
            assessed += 1
            major = max(counts, key=counts.get)
            for alt, count in counts.items():
                fraction = count / depth
                if alt == major or fraction < MIN_ALT_FRACTION:
                    continue
                # A major-base inference is not a verified reference allele.
                # Only annotate panel calls after checking the FASTA reference.
                calls.append({'position': column.reference_pos + 1, 'major': major,
                              'alt': alt, 'depth': depth, 'alt_reads': count,
                              'frequency': round(fraction * 100, 2),
                              'gene': '', 'mutation': '', 'interpretation': 'Unclassified', 'evidence': ''})
    # Read-level evidence must be available BEFORE annotating panel mutations.
    try:
        from .cmv_read_codons import audit_read_codons
        phased_calls, phased_note, codon_qc = audit_read_codons(bam_file, ref_path, ref_name, include_qc=True)
    except Exception as exc:
        phased_calls, phased_note, codon_qc = [], 'READ_CODON_AUDIT_FAILED: ' + str(exc), []
    from .cmv_codon import assess_protein_variants
    protein_result = assess_protein_variants(ref_path, ref_name, site_counts, covered_positions, indels,
                                             codon_qc=codon_qc, phased_calls=phased_calls)
    if not codon_qc:
        protein_result['cmv_protein_mutation_matches'] = []
        protein_result['cmv_interpretation'] = 'NOT ASSESSABLE: read-level codon evidence unavailable.'
    from .cmv_comparison import compare_calls
    summary, details = compare_calls(protein_result.get('cmv_protein_variant_calls', []), phased_calls)
    protein_result['cmv_read_codon_calls'] = summary
    protein_result['cmv_read_codon_diagnostics'] = details
    protein_result['cmv_read_codon_status'] = phased_note
    protein_result['cmv_read_codon_qc'] = codon_qc
    return {**protein_result,'cmv_variant_calls': calls, 'cmv_indel_calls': indels,
            'cmv_indel_status': 'Exploratory read-supported indels; no protein consequence or resistance interpretation.',
            'cmv_assessed_positions': assessed,
            'cmv_panel_status': note, 'cmv_panel_entries': len(panel)}


# Standardized Merlin remapping and reference validation (formerly cmv_standard.py)
"""CMV-only research screening against a fixed, annotated Merlin reference.

Strain selection/consensus remains independent. This module remaps raw reads to
NC_006273.2 solely for exploratory gene-specific variant screening.
"""
import csv
import os
from pathlib import Path
import subprocess

MERLIN = 'NC_006273.2'
GENES = {'UL54', 'UL97'}


def _fasta_ids(path):
    with open(path, encoding='ascii') as handle:
        for line in handle:
            if line.startswith('>'):
                yield line[1:].split()[0]


def _find_merlin(reference_folder):
    hits = []
    for path in sorted(Path(reference_folder).iterdir()):
        if path.name.startswith('.') or path.suffix.lower() not in ('.fa', '.fasta', '.fna'):
            continue
        try:
            if MERLIN in set(_fasta_ids(path)):
                hits.append(path)
        except (UnicodeError, OSError):
            continue
    if len(hits) != 1:
        raise RuntimeError(f'Expected exactly one Merlin FASTA ({MERLIN}), found {len(hits)} in {reference_folder}')
    return hits[0]


def _check_coordinates(reference_folder):
    path = Path(reference_folder) / 'cmv_gene_coordinates.tsv'
    if not path.is_file():
        raise RuntimeError(f'Coordinate table missing: {path}')
    with path.open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    genes = {r.get('gene') for r in rows if r.get('accession') == MERLIN}
    if not GENES.issubset(genes):
        raise RuntimeError(f'Merlin coordinates missing: {sorted(GENES - genes)}')


def analyse_cmv_standard(sample, input_fastq, reference_folder, sample_dir):
    """Return CMV analysis fields; raises on missing tools/annotation or alignment failure."""
    ref = _find_merlin(reference_folder)
    _check_coordinates(reference_folder)
    prefix = Path(sample_dir) / f'{sample}.cmv_merlin_resistance'
    bam = str(prefix) + '.sorted.bam'
    align = subprocess.Popen(['minimap2', '-ax', 'map-ont', str(ref), str(input_fastq)],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        sorted_result = subprocess.run(['samtools', 'sort', '-o', bam, '-'],
                                       stdin=align.stdout, capture_output=True)
    finally:
        if align.stdout:
            align.stdout.close()
    minimap_err = align.stderr.read().decode('utf-8', 'replace') if align.stderr else ''
    rc = align.wait()
    if rc or sorted_result.returncode:
        raise RuntimeError('Merlin alignment failed: ' + (minimap_err + ' ' + sorted_result.stderr.decode('utf-8', 'replace'))[-1200:])
    subprocess.run(['samtools', 'index', bam], check=True)
    result = analyse_cmv_mutations(sample, bam, MERLIN, str(ref))
    result['cmv_resistance_reference'] = MERLIN
    result['cmv_resistance_reference_file'] = str(ref)
    result['cmv_resistance_bam'] = bam
    result['cmv_resistance_mapping_status'] = 'STANDARD_MERLIN_MAPPING'
    return result
