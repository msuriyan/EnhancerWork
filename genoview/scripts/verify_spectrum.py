#!/usr/bin/env python3
"""Independent (pysam) re-calculation of the CRISPR deletion / indel size distribution,
used to check GenoView's numbers.  Rules (same definitions as GenoView, written separately):
  * a read = its primary + supplementary alignments (same chromosome and strand);
  * the read counts if its aligned span covers cut-W .. cut+W;
  * deletions: CIGAR D, plus gaps between split (supplementary) alignments;
    insertions: CIGAR I;  an event counts if it overlaps cut-W .. cut+W;
  * each distinct event is counted once per read;
  * % of deletions = deletions of that size / all counted deletions x 100 (sums to 100).
usage: verify_spectrum.py --cut 100031 --window 5 [--min-indel 1] [--kinds D|DI] a.bam b.bam ...
"""
import argparse, collections, json, pysam, sys
ap = argparse.ArgumentParser(); ap.add_argument('--cut', type=int, required=True, help='0-based base after which the cut falls (GenoView cut coordinate)')
ap.add_argument('--chrom', default='chr1'); ap.add_argument('--window', type=int, default=5); ap.add_argument('--flank', type=int, default=3000)
ap.add_argument('--min-indel', type=int, default=1); ap.add_argument('--kinds', default='D'); ap.add_argument('--mapq', type=int, default=20)
ap.add_argument('--json', action='store_true'); ap.add_argument('bams', nargs='+'); a = ap.parse_args()
A, B = a.cut - a.window, a.cut + a.window
out = {}
for path in a.bams:
    bam = pysam.AlignmentFile(path)
    segs = collections.defaultdict(list)
    for r in bam.fetch(a.chrom, max(0, a.cut - a.flank), a.cut + a.flank):
        if r.is_unmapped or r.is_secondary or r.is_duplicate or r.mapping_quality < a.mapq: continue
        segs[r.query_name].append(r)
    n = 0; dels = collections.Counter(); ins = collections.Counter(); reads_with = 0
    for name, rs in segs.items():
        prim = [r for r in rs if not r.is_supplementary]
        if not prim: continue
        p = prim[0]
        parts = sorted([r for r in rs if r.is_reverse == p.is_reverse], key=lambda r: r.reference_start)
        s0, s1 = parts[0].reference_start, max(r.reference_end for r in parts)
        if not (s0 <= A and s1 >= B): continue
        n += 1
        ev = set()
        for r in parts:
            pos = r.reference_start
            for op, l in r.cigartuples:
                if op in (0, 7, 8): pos += l
                elif op in (2, 3): ev.add(('D', pos, l)); pos += l
                elif op == 1: ev.add(('I', pos, l))
        for x, y in zip(parts, parts[1:]):
            if y.reference_start > x.reference_end: ev.add(('D', x.reference_end, y.reference_start - x.reference_end))
        touch = lambda k, pos, l: (pos <= a.cut + 1 and pos + l >= a.cut - 1) if k == 'D' else abs(pos - a.cut) <= 1
        kept = [(k, pos, l) for k, pos, l in ev
                if (k == 'D' and pos < B and pos + l > A) or (k == 'I' and A <= pos <= B)]
        kept = [e for e in kept if (a.kinds != 'D' or e[0] == 'D') and (e[2] >= a.min_indel or touch(*e))]
        for k, pos, l in kept: (dels if k == 'D' else ins)[l] += 1
        if kept: reads_with += 1
    nd = sum(dels.values()); ne = nd + sum(ins.values())
    out[path.split('/')[-1]] = {'reads': n, 'readsWith': reads_with, 'del': dict(sorted(dels.items())), 'ins': dict(sorted(ins.items())), 'deletions': nd, 'events': ne}
    if not a.json:
        print(f"{path}: {n} spanning reads, {reads_with} with {'a deletion' if a.kinds == 'D' else 'an indel'}, {nd} deletions")
        for l, c in sorted(dels.items()): print(f"  del {l:>6} bp  {c:>5}  {100*c/nd:6.2f}% of deletions  {100*c/n:6.2f}% of reads")
if a.json: print(json.dumps(out))
