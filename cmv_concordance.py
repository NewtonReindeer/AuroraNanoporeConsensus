#!/usr/bin/env python3
"""Independent, standard-library CMV concordance audit. Does not change variant calls."""
import argparse, csv, re
from collections import defaultdict
from pathlib import Path

MUT = re.compile(r'^([A-Z*])(\d+)([A-Z*])$')

def rows(path):
    with open(path, newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def normalize(value):
    return str(value or '').strip().upper().replace(' ', '')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--expected', required=True, help='Manually verified comparison results TSV')
    p.add_argument('--observed', required=True, help='Aurora protein variant TSV export')
    p.add_argument('--coverage', help='Optional gene coverage TSV with sample,gene,status')
    p.add_argument('--out', required=True)
    a=p.parse_args()
    expected=rows(a.expected); observed=rows(a.observed)
    coverage={(normalize(r.get('sample')),normalize(r.get('gene'))):normalize(r.get('status')) for r in rows(a.coverage)} if a.coverage else {}
    by_gene=defaultdict(list)
    for r in observed:
        key=(normalize(r.get('sample')),normalize(r.get('gene')))
        if not MUT.fullmatch(normalize(r.get('mutation'))):
            raise ValueError('Invalid observed mutation: '+repr(r))
        by_gene[key].append(r)
    results=[]
    for e in expected:
        sample,gene=normalize(e.get('sample')),normalize(e.get('gene'))
        mutation=normalize(e.get('mutation'))
        if not sample or not gene or not mutation:
            raise ValueError('Expected rows require sample, gene, mutation')
        if mutation not in ('NONE','UNKNOWN') and not MUT.fullmatch(mutation):
            raise ValueError('Expected mutation must be a single amino acid substitution or NONE/UNKNOWN: '+mutation)
        calls=by_gene.get((sample,gene),[])
        matches=[r for r in calls if normalize(r.get('mutation'))==mutation]
        cov=coverage.get((sample,gene),'UNKNOWN')
        if mutation=='UNKNOWN': verdict='NOT_COMPARABLE'
        elif mutation=='NONE': verdict='UNEXPECTED_CALLS' if calls else ('NO_CALLS_REPORTED' if cov=='UNKNOWN' else 'NO_CALLS_'+cov)
        elif matches: verdict='CONCORDANT_CANDIDATE'
        else: verdict='NOT_DETECTED_COVERAGE_UNKNOWN' if cov=='UNKNOWN' else 'NOT_DETECTED_'+cov
        results.append(dict(sample=sample,gene=gene,expected_mutation=mutation,expected_fraction=e.get('expected_fraction',''),observed_mutations=';'.join(sorted(set(normalize(r.get('mutation')) for r in calls))),observed_frequency=';'.join(str(r.get('frequency','')) for r in matches),gene_status=cov,comparison=verdict,notes=e.get('notes','')))
    out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(results[0]) if results else ['sample','gene','comparison'],delimiter='\t');writer.writeheader();writer.writerows(results)
    print('Wrote',len(results),'comparisons to',out)
    print('Candidate matches:',sum(r['comparison']=='CONCORDANT_CANDIDATE' for r in results))
    print('Not detected:',sum(r['comparison'].startswith('NOT_DETECTED') for r in results))

if __name__=='__main__':main()
