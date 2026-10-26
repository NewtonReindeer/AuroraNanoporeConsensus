"""Read-only BAM audit without pysam. Inspects header/reference consistency and indel support.
BAM coordinate positions are 1-based in outputs. This does NOT call resistance.
"""
import argparse,csv,gzip,json,struct
from collections import Counter,defaultdict
from pathlib import Path
OPS='MIDNSHP=XB'

def records(path):
    with gzip.open(path,'rb') as f:
        if f.read(4)!=b'BAM\1': raise ValueError('Not a BAM')
        n=struct.unpack('<i',f.read(4))[0]; header=f.read(n).decode(errors='replace')
        refs=[]
        for _ in range(struct.unpack('<i',f.read(4))[0]):
            m=struct.unpack('<i',f.read(4))[0]; name=f.read(m).rstrip(b'\0').decode(); length=struct.unpack('<i',f.read(4))[0]; refs.append((name,length))
        yield ('header',header,refs)
        while True:
            b=f.read(4)
            if not b: break
            if len(b)!=4: raise ValueError('Truncated BAM')
            sz=struct.unpack('<i',b)[0]; raw=f.read(sz)
            if len(raw)!=sz: raise ValueError('Truncated BAM record')
            rid,pos,bin_mq_nl,flag_nc,lseq,_,_,_=struct.unpack_from('<iiIIiiii',raw)
            nl=bin_mq_nl&255; mq=(bin_mq_nl>>8)&255; flag=flag_nc>>16; nc=flag_nc&65535
            off=32+nl; cigar=struct.unpack_from('<'+'I'*nc,raw,off) if nc else ()
            yield ('read',rid,pos,mq,flag,[(OPS[x&15],x>>4) for x in cigar])

def audit(bam,minimum_mapq=20):
    refs=[]; header=''; events=defaultdict(Counter); coverage=Counter(); total=Counter()
    for rec in records(bam):
        if rec[0]=='header': header,refs=rec[1:];continue
        _,rid,pos,mq,flag,cigar=rec
        total['all']+=1
        if flag & (0x4|0x100|0x200|0x400|0x800) or mq<minimum_mapq: continue
        total['eligible']+=1
        cursor=pos
        for op,n in cigar:
            if op=='I': events[(rid,cursor,'insertion',n)]['support']+=1
            elif op=='D': events[(rid,cursor+1,'deletion',n)]['support']+=1
            if op in 'MDN=X': cursor+=n
        coverage[rid]+=1
    return {'bam':Path(bam).name,'reference':refs,'header_reference_consistent':all(str(x[0]) in header for x in refs),
            'counts':dict(total),'indels':[{'reference':refs[rid][0], 'position_1based':p,
            'type':typ,'length':n,'supporting_reads':counts['support']} for (rid,p,typ,n),counts in
            sorted(events.items(),key=lambda x:-x[1]['support'])]}

def main():
    p=argparse.ArgumentParser();p.add_argument('bams',nargs='+');p.add_argument('--out',default='bam_reference_audit.json');args=p.parse_args()
    results=[audit(x) for x in args.bams];Path(args.out).write_text(json.dumps(results,indent=2),encoding='utf8')
    for r in results:
        print(r['bam'],r['reference'],r['counts'])
        print('top indels',r['indels'][:8])
if __name__=='__main__':main()
