"""Read-only CMV report discrepancy triage. No clinical calls are changed.

Usage: python cmv_discrepancy_audit.py --reports /path/to/CMV_reports --out audit.tsv
"""
import argparse
import csv
from pathlib import Path
from html.parser import HTMLParser

class TableReader(HTMLParser):
    def __init__(self):
        super().__init__(); self.tables=[]; self.table=None; self.row=None; self.cell=None; self.section=''
        self.heading=False; self.heading_text=''
    def handle_starttag(self, tag, attrs):
        if tag=='h2': self.heading=True; self.heading_text=''
        if tag=='table': self.table=[]
        if tag=='tr' and self.table is not None: self.row=[]
        if tag in ('td','th') and self.row is not None: self.cell=[]
    def handle_data(self,data):
        if self.heading: self.heading_text+=data
        if self.cell is not None: self.cell.append(data)
    def handle_endtag(self,tag):
        if tag=='h2': self.heading=False; self.section=self.heading_text.strip()
        if tag in ('td','th') and self.cell is not None:
            self.row.append(''.join(self.cell).strip()); self.cell=None
        if tag=='tr' and self.row is not None:
            self.table.append(self.row); self.row=None
        if tag=='table' and self.table is not None:
            self.tables.append((self.section,self.table)); self.table=None

def numeric(s):
    try: return float(str(s).replace('%','').replace(',','.').replace('×','').strip())
    except (ValueError,TypeError): return None

def inspect(path):
    parser=TableReader(); parser.feed(path.read_text(encoding='utf-8'))
    result=[]
    for section,table in parser.tables:
        if not table: continue
        hdr=table[0]
        if 'Gen' not in hdr or not any('Kodon' in x for x in hdr): continue
        for cells in table[1:]:
            if len(cells)!=len(hdr): continue
            r=dict(zip(hdr,cells)); gene=r.get('Gen','')
            if gene not in ('UL54','UL97','UL56'): continue
            total_field=next((v for k,v in r.items() if 'Bedömbara / totalt' in k),None)
            if total_field and '/' in total_field:
                try: assessed,total=[int(x.strip()) for x in total_field.split('/')]
                except ValueError: continue
                unclear=numeric(r.get('Oklar läsram',r.get('Oklart läsramsstöd')))
                covered=numeric(r.get('Kodon ≥30×'))
                if assessed>total: result.append((path.stem,gene,'ERROR','Bedömbara kodon > totalt',f'{assessed}/{total}'))
                if unclear is not None and unclear>0:
                    result.append((path.stem,gene,'REVIEW','Kodon med oklart läsramsstöd',f'{int(unclear)}; bedömbara {assessed}/{total}; kodon ≥30× {covered}'))
                if covered is not None and assessed < covered:
                    result.append((path.stem,gene,'REVIEW','Täckning skiljer sig från bedömbarhet',f'{int(covered)} täckta kodon; {assessed} bedömbara'))
    return result

def main():
    p=argparse.ArgumentParser(); p.add_argument('--reports',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
    args=p.parse_args(); findings=[]
    for f in sorted(args.reports.glob('*CMV.html')): findings.extend(inspect(f))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('w',newline='',encoding='utf-8') as h:
        w=csv.writer(h,delimiter='\t'); w.writerow(['sample','gene','severity','issue','detail']); w.writerows(findings)
    print(f'Inspected {len(list(args.reports.glob("*CMV.html")))} reports; {len(findings)} flags. Output: {args.out}')
if __name__=='__main__': main()
