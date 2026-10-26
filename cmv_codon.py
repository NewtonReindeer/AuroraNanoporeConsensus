"""Reference-specific, conservative CMV single-nucleotide codon analysis.

Requires cmv_gene_coordinates.tsv alongside the selected reference FASTA.
Columns: accession, gene, start, end, strand, first_codon_number.
Coordinates are 1-based inclusive; first_codon_number is normally 1.
Only contiguous, complete CDS with verified coordinates are supported.
No clinical resistance call is made without verified gene mapping.
"""
import csv
from pathlib import Path
from collections import defaultdict
from .cmv_thresholds import MIN_ALT_FRACTION, MIN_CLEAN_DEPTH, MIN_ALT_READS

_COMPLEMENT = str.maketrans('ACGT', 'TGCA')
_BASES = 'ACGT'
_CODONS = {c: aa for aa, codons in {
 'F':['TTT','TTC'],'L':['TTA','TTG','CTT','CTC','CTA','CTG'],
 'I':['ATT','ATC','ATA'],'M':['ATG'],'V':['GTT','GTC','GTA','GTG'],
 'S':['TCT','TCC','TCA','TCG','AGT','AGC'],'P':['CCT','CCC','CCA','CCG'],
 'T':['ACT','ACC','ACA','ACG'],'A':['GCT','GCC','GCA','GCG'],
 'Y':['TAT','TAC'],'H':['CAT','CAC'],'Q':['CAA','CAG'],
 'N':['AAT','AAC'],'K':['AAA','AAG'],'D':['GAT','GAC'],
 'E':['GAA','GAG'],'C':['TGT','TGC'],'W':['TGG'],
 'R':['CGT','CGC','CGA','CGG','AGA','AGG'],'G':['GGT','GGC','GGA','GGG'],
 '*':['TAA','TAG','TGA']}.items() for c in codons}

def _default(note):
    return {'cmv_protein_mutation_matches': [], 'cmv_protein_variant_calls': [],
            'cmv_gene_coverage': [], 'cmv_published_mutation_status': note,
            'cmv_interpretation': 'NOT ASSESSED: '+note}

def assess_protein_variants(reference_path, reference_name, site_counts, covered_positions, indels=None, codon_qc=None, phased_calls=None):
    """Translate single alternate bases on complete, correctly annotated codons.

    Independent SNVs in one codon are not phased: compound codon changes
    cannot be inferred and must be manually reviewed.
    """
    import pysam
    ref = Path(reference_path)
    coords = ref.parent / 'cmv_gene_coordinates.tsv'
    if not coords.is_file():
        return _default('Missing accession-specific cmv_gene_coordinates.tsv beside reference FASTA.')
    try:
        qc_index = {(str(r.get("gene")), int(r.get("codon"))): r for r in (codon_qc or [])}
        phased_index = {(str(r.get("gene")), str(r.get("mutation")), str(r.get("alternate_codon"))): r for r in (phased_calls or [])}
        with coords.open(encoding='utf-8-sig', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            required = {'accession','gene','start','end','strand','first_codon_number'}
            if not required.issubset(reader.fieldnames or []):
                return _default('Gene-coordinate TSV missing required columns.')
            entries = [r for r in reader if r['accession'].strip() == reference_name]
        if not entries:
            return _default('No coordinate entries match selected reference sequence ID '+reference_name)
        gene_rows, variants, issues = [], [], []
        indels = indels or []
        with pysam.FastaFile(str(ref)) as fasta:
            for r in entries:
                gene = r['gene'].strip().upper()
                if gene not in ('UL54','UL97','UL56','UL89','UL51'):
                    continue
                try:
                    start, end = int(r['start']), int(r['end'])
                    offset = int(r['first_codon_number'])
                    strand = r['strand'].strip()
                    if start < 1 or end < start or (end-start+1)%3 or strand not in ('+','-') or offset < 1:
                        raise ValueError('invalid CDS geometry')
                    seq = fasta.fetch(reference_name, start-1, end).upper()
                    if len(seq) != end-start+1 or any(c not in _BASES for c in seq):
                        raise ValueError('missing/ambiguous reference bases')
                    if strand == '-': seq = seq.translate(_COMPLEMENT)[::-1]
                    covered = sum(start <= p <= end for p in covered_positions)
                    assessable_codons = 0
                    covered_codons = 0
                    excluded_codons = 0
                    multi_snv_codons = 0
                    mixed_frame_codons = 0
                    for j in range(0, len(seq), 3):
                        codon = seq[j:j+3]
                        pos = [start+j+k if strand == '+' else end-j-k for k in range(3)]
                        if not all(p in site_counts for p in pos): continue
                        covered_codons += 1
                        # Distinguish local indel ambiguity from downstream mixed
                        # reading-frame uncertainty. A minority frameshift does
                        # not imply every downstream molecule is frameshifted.
                        # However, unphased downstream SNVs must not be treated
                        # as independently verified resistance substitutions.
                        codon_lo, codon_hi = min(pos), max(pos)
                        local_indel = False
                        downstream_mixed_frame = False
                        for event in indels:
                            anchor = int(event['position'])
                            length = int(event['length'])
                            if not (start <= anchor <= end):
                                continue
                            # Indels at/near this codon invalidate its direct
                            # isolated-SNV translation regardless of frequency.
                            if codon_lo - 3 <= anchor <= codon_hi + 3:
                                local_indel = True
                                break
                            if length % 3:
                                downstream = ((strand == '+' and anchor < codon_lo)
                                              or (strand == '-' and anchor > codon_hi))
                                if downstream:
                                    downstream_mixed_frame = True
                        if local_indel:
                            excluded_codons += 1
                            continue
                        if downstream_mixed_frame:
                            mixed_frame_codons += 1
                        # Sample-wide indels are NOT evidence that every read at this codon
                        # is out of frame. Use read-specific clean-depth when available.
                        qc = qc_index.get((gene, offset+j//3))
                        if qc is not None:
                            clean_assessable = int(qc.get("clean_depth", 0)) >= 30
                        else:
                            clean_assessable = not downstream_mixed_frame
                        # Independent alternate alleles at two positions of one
                        # codon cannot be translated as an observed haplotype.
                        alternate_sites = sum(
                            any((b if strand == '+' else b.translate(_COMPLEMENT)) != codon[k]
                                and n >= MIN_ALT_READS and n / depth >= MIN_ALT_FRACTION
                                for b, n in site_counts[pos[k]][0].items())
                            for k, depth in enumerate(site_counts[q][1] for q in pos)
                        )
                        if alternate_sites > 1:
                            multi_snv_codons += 1
                            excluded_codons += 1
                            continue
                        if clean_assessable:
                            assessable_codons += 1
                        ref_aa = _CODONS.get(codon)
                        if ref_aa is None or ref_aa == '*': continue
                        for k, p in enumerate(pos):
                            counts, depth = site_counts[p]
                            ref_base = codon[k]
                            for observed_base, n in counts.items():
                                coding_base = observed_base if strand == '+' else observed_base.translate(_COMPLEMENT)
                                if coding_base == ref_base or n < MIN_ALT_READS or n/depth < MIN_ALT_FRACTION: continue
                                changed = codon[:k]+coding_base+codon[k+1:]
                                aa = _CODONS.get(changed)
                                if aa is None or aa == ref_aa: continue
                                variants.append({'gene':gene, 'mutation':f'{ref_aa}{offset+j//3}{aa}',
                                  'position':p, 'depth':depth, 'alt_reads':n,
                                  'frequency':round(100*n/depth,2), 'allele_class':('DOMINANT_REFERENCE_DIFFERENCE' if n == max(counts.values()) else 'MINORITY_ALLELE'), 'reference_codon':codon,
                                  'alternate_codon':changed, 'status':('EXPLORATORY_UNPHASED_SNV' if clean_assessable else 'INSUFFICIENT_CLEAN_READ_EVIDENCE') })
                    gene_rows.append({'gene':gene,'covered_bases_30x':covered,'length':len(seq),
                       'coverage_pct':round(100*covered/len(seq),2),
                       'assessable_codons_30x':assessable_codons,'covered_codons_30x':covered_codons,
                       'excluded_codons':excluded_codons,'multi_snv_codons':multi_snv_codons,
                       'mixed_frame_codons':mixed_frame_codons,
                       'total_codons':len(seq)//3,
                       'status':('NOT_COVERED' if covered_codons == 0 else 'EXCLUDED_BY_VARIANT_COMPLEXITY' if assessable_codons == 0 else 'PARTIAL' if assessable_codons < len(seq)//3 else 'COVERED')})
                except (ValueError, KeyError, IndexError) as e:
                    issues.append(f'{gene}: {e}')
        annotated = annotate_protein_variants(variants)
        # A panel match requires a matching, read-phased codon with sufficient clean
        # evidence. Unphased site-count translations remain exploratory diagnostics.
        for item in annotated:
            evidence = phased_index.get((item['gene'], item['mutation'], item['alternate_codon']))
            verified = (evidence is not None and
                        evidence.get('status') == 'READ_PHASED_EXPLORATORY' and
                        int(evidence.get('depth', 0)) >= MIN_CLEAN_DEPTH and
                        int(evidence.get('alt_reads', 0)) >= MIN_ALT_READS and
                        float(evidence.get('frequency', 0)) >= 100 * MIN_ALT_FRACTION)
            item['read_phased_support'] = bool(verified)
            if verified:
                item['depth'] = evidence['depth']
                item['alt_reads'] = evidence['alt_reads']
                item['frequency'] = evidence['frequency']
                item['status'] = 'READ_PHASED_PANEL_CANDIDATE' if item.get('panel_match') else 'READ_PHASED_EXPLORATORY'
            else:
                item['status'] = 'NOT_CONFIRMED_BY_CLEAN_READS'
        matched = [v for v in annotated if v.get('panel_match') and v.get('read_phased_support')]
        return {'cmv_protein_mutation_matches':matched, 'cmv_protein_variant_calls':annotated,
            'cmv_gene_coverage':gene_rows,
            'cmv_published_mutation_status':('Exploratory codon-level SNV screening; coordinate TSV must be independently verified. '
                'Local indel-affected codons excluded; downstream mixed-frame candidates flagged and withheld from published mutation matches. No phased haplotype or indel protein effects. '+('; '.join(issues) if issues else '')),
            'cmv_interpretation':'EXPLORATORY ONLY: protein SNVs assessed where complete codons have >=30x; no clinical resistance conclusion.'}
    except (OSError, ValueError, KeyError) as e:
        return _default('Unable to load CMV gene mapping: '+str(e))


# Published mutation catalog lookup (formerly cmv_resistance.py)
"""CMV published amino-acid mutation lookup. No third-party dependencies.

This module matches independently called, gene-specific protein substitutions.
It deliberately does not infer protein positions from genomic SNV positions.
"""
import csv
from pathlib import Path

PANEL = Path(__file__).resolve().parent / 'cmv_published_mutations.tsv'


def load_published_mutations(path=None):
    path = Path(path) if path else PANEL
    with path.open(encoding='utf-8-sig', newline='') as f:
        return {(r['gene'].upper(), r['mutation'].upper()): r
                for r in csv.DictReader(f, delimiter='\t')}


def annotate_protein_variants(variants, path=None):
    """Input: iterable of dicts with gene, mutation, and optionally frequency/depth.

    Caller MUST supply independently verified protein substitutions with correct
    accession-specific codon numbering. No interpretation of indels or mixtures.
    """
    panel = load_published_mutations(path)
    results = []
    for variant in variants:
        gene = str(variant.get('gene', '')).upper().strip()
        mutation = str(variant.get('mutation', '')).upper().strip()
        record = panel.get((gene, mutation))
        if record:
            results.append({**variant, 'panel_match': True,
                            'drug_codes': record['drug_codes'],
                            'source_row': record['source_row'],
                            'comment': record['comment'],
                            'source_status': record['status'], 'status': variant.get('status', 'UNREVIEWED')})
        else:
            results.append({**variant, 'panel_match': False, 'drug_codes': '',
                            'source_row': '', 'comment': '',
                            'source_status': 'NOT_LISTED_IN_UPLOADED_WORKBOOK', 'status': variant.get('status', 'UNREVIEWED')})
    return results
