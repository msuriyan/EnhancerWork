#!/usr/bin/env python3
"""
deletion_barplot.py — % of reads carrying a deletion at a CRISPR cut site, per BAM, as a bar plot.

For each BAM:
  spanning reads  = primary, mapped, MAPQ >= --mapq reads whose alignment (incl. split parts) covers cut ± window
  deletion reads  = spanning reads with >= 1 deletion (CIGAR 'D', or the gap between split alignments on the
                    same chromosome/strand) overlapping cut ± window, of length within --min-del .. --max-del
  % = deletion reads / spanning reads × 100

usage:
  python deletion_barplot.py --cut chr1:100031 *.bam -o deletion_pct
  python deletion_barplot.py --cut chr1:100031 --window 10 --min-del 3 --names names.tsv --groups groups.tsv *.bam
outputs: deletion_pct.tsv, deletion_pct.png (and .svg, .pdf)
requires: pysam, matplotlib
"""
import argparse, os, sys
import pysam
import logging
import matplotlib
matplotlib.use('Agg')
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
import matplotlib.pyplot as plt


def count(bam_path, chrom, cut, win, mapq, min_del, max_del, flank, skip_dups):
    lo, hi = cut - win, cut + win                     # 0-based window [lo, hi)
    span = dele = 0
    with pysam.AlignmentFile(bam_path) as bam:
        if chrom not in bam.references:
            alt = chrom[3:] if chrom.startswith('chr') else 'chr' + chrom
            if alt not in bam.references: sys.exit(f'{bam_path}: chromosome {chrom} not found')
            chrom = alt
        seen = set()
        for r in bam.fetch(chrom, max(0, cut - flank), cut + flank):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail: continue
            if skip_dups and r.is_duplicate: continue
            if r.mapping_quality < mapq or r.query_name in seen: continue
            seen.add(r.query_name)
            s0, s1 = r.reference_start, r.reference_end
            dels = []
            pos = r.reference_start
            for op, ln in r.cigartuples:
                if op in (0, 7, 8): pos += ln
                elif op in (2, 3): dels.append((pos, ln)); pos += ln
            # split read: supplementary on same chrom/strand separated by a reference gap = large deletion
            if r.has_tag('SA'):
                for sa in r.get_tag('SA').rstrip(';').split(';'):
                    c, p, strand, cig, mq = sa.split(',')[:5]
                    if c != chrom or int(mq) < mapq or (strand == '-') != r.is_reverse: continue
                    p = int(p) - 1
                    rl = sum(int(n) for n, o in __import__('re').findall(r'(\d+)([MDN=X])', cig))
                    if p >= s1: dels.append((s1, p - s1)); s1 = p + rl
                    elif p + rl <= s0: dels.append((p + rl, s0 - p - rl)); s0 = p
            if not (s0 <= lo and s1 >= hi): continue       # must span the window
            span += 1
            if any(min_del <= ln <= max_del and a < hi and a + ln > lo for a, ln in dels): dele += 1
    return span, dele


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('bams', nargs='+')
    ap.add_argument('--cut', required=True, help='chrom:pos — cut is between pos and pos+1 (1-based), e.g. chr1:100031')
    ap.add_argument('--window', type=int, default=5, help='bp either side of the cut (default 5; ~10 for long reads)')
    ap.add_argument('--mapq', type=int, default=20)
    ap.add_argument('--min-del', type=int, default=1, help='smallest deletion counted (e.g. 3 for noisy long reads)')
    ap.add_argument('--max-del', type=int, default=10**9)
    ap.add_argument('--flank', type=int, default=20000, help='region fetched around the cut (catches split reads)')
    ap.add_argument('--keep-dups', action='store_true')
    ap.add_argument('--names', help='TSV: bam_file_or_basename<TAB>label (x-axis names)')
    ap.add_argument('--groups', help='TSV: label<TAB>group (bars coloured by group)')
    ap.add_argument('--order', choices=['input', 'value', 'name'], default='input')
    ap.add_argument('--color', default='#4c72b0'); ap.add_argument('--palette', default='#4c72b0,#dd8452,#55a868,#c44e52,#8172b3,#937860')
    ap.add_argument('--title', default=''); ap.add_argument('--ylabel', default='Reads with deletion (%)'); ap.add_argument('--xlabel', default='')
    ap.add_argument('--font', default='Arial'); ap.add_argument('--fontsize', type=float, default=11)
    ap.add_argument('--width', type=float, default=0, help='inches (0 = auto)'); ap.add_argument('--height', type=float, default=3.5)
    ap.add_argument('--dpi', type=int, default=300); ap.add_argument('--ymax', type=float, default=0)
    ap.add_argument('--no-values', action='store_true', help='do not print % above bars')
    ap.add_argument('-o', '--out', default='deletion_pct')
    a = ap.parse_args()

    chrom, pos = a.cut.rsplit(':', 1); cut = int(pos.replace(',', ''))   # 0-based index of first base after the cut
    names = {}
    if a.names:
        for ln in open(a.names):
            if ln.strip() and not ln.startswith('#'): k, v = ln.rstrip('\n').split('\t')[:2]; names[k] = v
    rows = []
    for b in a.bams:
        span, dele = count(b, chrom, cut, a.window, a.mapq, a.min_del, a.max_del, a.flank, not a.keep_dups)
        base = os.path.basename(b)
        label = names.get(b, names.get(base, base.removesuffix('.bam')))
        pct = 100 * dele / span if span else float('nan')
        rows.append((label, b, span, dele, pct))
        print(f'{label:30s} spanning {span:8d}  with deletion {dele:8d}  {pct:6.2f}%')
    if a.order == 'value': rows.sort(key=lambda r: -r[4])
    elif a.order == 'name': rows.sort(key=lambda r: r[0])

    with open(a.out + '.tsv', 'w') as f:
        f.write('sample\tbam\tspanning_reads\treads_with_deletion\tpercent\n')
        for r in rows: f.write(f'{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.4f}\n')

    groups = {}
    if a.groups:
        for ln in open(a.groups):
            if ln.strip() and not ln.startswith('#'): k, v = ln.rstrip('\n').split('\t')[:2]; groups[k] = v
    pal = a.palette.split(','); gl = list(dict.fromkeys(groups.get(r[0], '') for r in rows))
    colors = [pal[gl.index(groups[r[0]]) % len(pal)] if r[0] in groups else a.color for r in rows]

    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': [a.font, 'Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'], 'font.size': a.fontsize, 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    w = a.width or max(3, 0.55 * len(rows) + 1.2)
    fig, ax = plt.subplots(figsize=(w, a.height))
    x = range(len(rows))
    bars = ax.bar(x, [r[4] for r in rows], color=colors, edgecolor='#333333', linewidth=0.6, width=0.75)
    if not a.no_values:
        for bar, r in zip(bars, rows): ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'{r[4]:.1f}', ha='center', va='bottom', fontsize=a.fontsize * 0.8)
    ax.set_xticks(list(x)); ax.set_xticklabels([r[0] for r in rows], rotation=45 if len(rows) > 4 else 0, ha='right' if len(rows) > 4 else 'center')
    ax.set_ylabel(a.ylabel); ax.set_xlabel(a.xlabel)
    ax.set_ylim(0, a.ymax or min(100, max([r[4] for r in rows if r[4] == r[4]] + [1]) * 1.15))
    if a.title: ax.set_title(a.title)
    ax.spines[['top', 'right']].set_visible(False)
    if groups:
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(facecolor=pal[i % len(pal)], edgecolor='#333333', label=g) for i, g in enumerate(gl) if g], frameon=False, fontsize=a.fontsize * 0.85)
    fig.tight_layout()
    for ext in ('png', 'svg', 'pdf'): fig.savefig(f'{a.out}.{ext}', dpi=a.dpi)
    print(f'wrote {a.out}.tsv/.png/.svg/.pdf')


if __name__ == '__main__':
    main()
