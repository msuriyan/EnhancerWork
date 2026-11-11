"""Generate synthetic test data (bigMaf, MAF, bigWig, bigBed, FASTA+fai, GFF3, Newick) with known ground truth.
Usage: python tests/make_test_data.py   ->  writes tests/data/
Requires: pip install pybigtools pyBigWig
"""
import os, pathlib, gzip, shutil
os.makedirs(pathlib.Path(__file__).resolve().parent / "data", exist_ok=True)
os.chdir(pathlib.Path(__file__).resolve().parent / "data")
import random, pybigtools, pyBigWig
random.seed(7)
B='ACGT'
chrom='NC_069100.1'; clen=3_000_000
sp=[('ilBany1.1',chrom,0.0),('ilArgPaph1.1','OX123456.1',0.05),('Danaus_plexippus','scaffold_3',0.12),('hg38','chr2',0.3),('speciesX','contig_7',0.2)]
ref=[random.choice(B) for _ in range(clen)]
open('ref.fa','w').write('>'+chrom+'\n'+'\n'.join(''.join(ref[i:i+60]) for i in range(0,clen,60))+'\n>chrUn\nACGTACGTNN\n')
# fai
import os
off=len('>'+chrom+'\n')

blocks=[]; pos=100_000; cnt={s[0]:1000 for s in sp}; truth={}  # truth[(sp,refpos)] = char
while pos<2_900_000:
    L=random.randint(200,3000)
    refseg=ref[pos:pos+L]
    cols=[]  # list of (refchar, {sp:char})
    for i,c in enumerate(refseg):
        if random.random()<0.004:  # insertion column(s) before this ref base
            for _ in range(random.randint(1,4)):
                cols.append(('-',{s[0]:(random.choice(B) if random.random()<0.5 else '-') for s in sp[1:]}))
        d={}
        for name,_,dv in sp[1:]:
            r=random.random()
            if r<dv*0.3: d[name]='-'
            elif r<dv: d[name]=random.choice(B)
            else: d[name]=c.lower() if random.random()<0.05 else c
        cols.append((c,d))
    present=[s for s in sp[1:] if random.random()>0.15]
    lines=['a score=1.0']
    rtxt=''.join(c for c,_ in cols)
    lines.append(f"s {sp[0][0]}.{chrom} {pos} {L} + {clen} {rtxt}")
    for name,ch,_ in present:
        txt=''.join(d[name] for _,d in cols)
        n=sum(1 for x in txt if x!='-')
        if n==0: continue
        strand='+' if random.random()<0.7 else '-'
        lines.append(f"s {name}.{ch} {cnt[name]} {n} {strand} 50000000 {txt}")
        cnt[name]+=n
        rp=pos
        for c,d in cols:
            if c!='-': truth[(name,rp)]=d[name]; rp+=1
    blocks.append((pos,pos+L,lines))
    pos+=L+random.randint(0,400)
with open('test.maf','w') as f:
    f.write('##maf version=1\n\n')
    for s,e,l in blocks: f.write('\n'.join(l)+'\n\n')
asql='''table bedMaf
"Bed3 with MAF block"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    lstring mafBlock;  "MAF block"
    )
'''
w=pybigtools.open('test.bigMaf.bb','w')
w.write({chrom:clen,'chrUn':10},[(chrom,s,e,';'.join(l)) for s,e,l in blocks],autosql=asql)
# truth sample for checks
import json
samp={f"{k[0]}|{k[1]}":v for k,v in random.sample(sorted(truth.items()),2000)}
json.dump(samp,open('truth.json','w'))
# bigwig
bw=pyBigWig.open('atac.bw','w'); bw.addHeader([(chrom,clen)],maxZooms=10)
starts=list(range(0,clen,25)); vals=[ (abs(random.gauss(0,1))+ (30 if (st//5000)%20==0 else 0)) for st in starts]
bw.addEntries([chrom]*len(starts),starts,ends=[s+25 for s in starts],values=vals); bw.close()
# narrowPeak + bigBed peaks
pk=[]
for i,st in enumerate(range(100000,2900000,100000)):
    pk.append((chrom,st,st+800,f'peak_{i}',500,'.',12.5,8.1,6.2,400))
open('peaks.narrowPeak','w').write(''.join('\t'.join(map(str,p))+'\n' for p in pk))
np_asql='''table narrowPeak
"BED6+4 Peaks"
(
string chrom; "chrom"
uint chromStart; "start"
uint chromEnd; "end"
string name; "name"
uint score; "score"
char[1] strand; "strand"
float signalValue; "sig"
float pValue; "p"
float qValue; "q"
int peak; "summit"
)
'''
w=pybigtools.open('peaks.bb','w'); w.write({chrom:clen},[(p[0],p[1],p[2],'\t'.join(map(str,p[3:]))) for p in pk],autosql=np_asql)
open('plain.bed','w').write(f"{chrom}\t150000\t151000\tmyRegion\t0\t+\n{chrom}\t152000\t153500\tother\t0\t-\n")
open('tree.nwk','w').write("((ilBany1.1:0.1,ilArgPaph1.1:0.1):0.1,(Danaus_plexippus:0.2,(hg38:0.3,speciesX:0.2):0.1):0.05);\n")
gff=["##gff-version 3",f"{chrom}\tt\tgene\t150001\t160000\t.\t+\t.\tID=g1;Name=engrailed",f"{chrom}\tt\tmRNA\t150001\t160000\t.\t+\t.\tID=m1;Parent=g1;Name=en-RA",
 f"{chrom}\tt\texon\t150001\t151000\t.\t+\t.\tParent=m1",f"{chrom}\tt\texon\t155001\t160000\t.\t+\t.\tParent=m1",f"{chrom}\tt\tCDS\t150501\t151000\t.\t+\t0\tParent=m1",f"{chrom}\tt\tCDS\t155001\t158000\t.\t+\t0\tParent=m1"]
open('genes.gff3','w').write('\n'.join(gff)+'\n')
print(len(blocks),'blocks')

# correct .fai
data=open('ref.fa','rb').read();recs=[];pos=0;cur=None
for line in data.splitlines(keepends=True):
    if line.startswith(b'>'):
        if cur: recs.append(cur)
        cur={'name':line[1:].strip().decode(),'len':0,'off':pos+len(line),'lb':None,'lw':None}
    else:
        if cur['lb'] is None: cur['lb']=len(line.rstrip(b'\n'));cur['lw']=len(line)
        cur['len']+=len(line.rstrip(b'\n'))
    pos+=len(line)
recs.append(cur)
open('ref.fa.fai','w').write(''.join(f"{r['name']}\t{r['len']}\t{r['off']}\t{r['lb']}\t{r['lw']}\n" for r in recs))
with open('test.maf','rb') as a, gzip.open('test.maf.gz','wb') as z: shutil.copyfileobj(a,z)
