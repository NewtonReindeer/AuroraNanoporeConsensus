"""Conservative, versioned CMV annotation overlay for reports only.

Never use this overlay to change variant calling or QC. Evidence must be curated
before promoting any record to clinically validated status.
"""
import json
from pathlib import Path

ANNOTATION_PATH = Path(__file__).with_name('cmv_clinical_annotations.json')


def load_annotations(path=None):
    with open(path or ANNOTATION_PATH, encoding='utf-8') as handle:
        data = json.load(handle)
    for x in data['mutations']:
        if not x.get('sources') or any(not (src.get('pmid') or src.get('doi')) for src in x['sources']):
            raise ValueError('Annotation lacks traceable PMID/DOI: %s %s' % (x.get('gene'), x.get('mutation')))
    return {(x['gene'].upper(), x['mutation'].upper()): x for x in data['mutations']}


def classify(gene, mutation, annotations=None):
    records = annotations if annotations is not None else load_annotations()
    record = records.get((str(gene).upper(), str(mutation).upper()))
    if record is None:
        return {'category': 'unknown', 'label': 'Klinisk betydelse ej fastställd',
                'drugs': 'Ej fastställt', 'source': 'Ingen granskad annotation',
                'note': 'Ingen kuraterad mutationsspecifik evidens tillgänglig.'}
    return {'category': record['category'], 'label': record['label'],
            'drugs': ', '.join(record.get('drugs', [])) or 'Ej fastställt',
            'source': '; '.join(str(src.get('pmid') or src.get('doi')) for src in record.get('sources', [])), 'note': record.get('note', '')}
