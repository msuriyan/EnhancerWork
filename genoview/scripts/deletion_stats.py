#!/usr/bin/env python3
"""
deletion_stats.py — statistics for CRISPR deletion-size distributions from GenoView exports.

Input : GenoView "Alleles TSV" files (CRISPR edits tab → Alleles TSV), one per sample/target,
        optionally a sample sheet (TSV: sample<TAB>group) to compare groups (e.g. guide, gene, condition).
Output: per_sample_summary.tsv, binned_distribution.tsv, group_tests.txt (in --out)

Unit of replication = the animal / sample (NOT the read): reads from one crispant are not independent,
so all group tests are done on per-sample values or with the sample as a cluster.

usage:
  python deletion_stats.py exports/*_alleles.tsv --groups samples.tsv --out results/
  python deletion_stats.py exports/*_alleles.tsv --out results/            # description only
requires: pandas, numpy, scipy (statsmodels optional, for the GEE models)
"""
import argparse, glob, os, re, sys
import numpy as np, pandas as pd
from scipy import stats

BIN_EDGES = [0, 20, 50, 100, 200, 500, 1000, 1500, 2000, np.inf]   # (lo, hi] in bp — same idea as GenoView's ranges
LARGE = 50                                                       # "large deletion" threshold (bp)


def read_alleles(paths):
    rows = []
    for p in paths:
        df = pd.read_csv(p, sep='\t', dtype={'allele': str})
        need = {'sample', 'allele', 'reads'}
        if not need <= set(df.columns):
            sys.exit(f'{p}: not a GenoView alleles TSV (needs columns {need})')
        df['file'] = os.path.basename(p)
        rows.append(df)
    df = pd.concat(rows, ignore_index=True)
    if 'target' not in df: df['target'] = 'all'
    df['sample'] = df['sample'].str.replace(r'\.bam$', '', regex=True)
    return df


def deletion_events(df):
    """one row per deletion event per allele, weighted by the allele's read count"""
    ev = []
    for r in df.itertuples(index=False):
        for m in re.finditer(r'D(\d+)@(-?\d+)', str(r.allele)):
            ev.append((r.sample, r.target, int(m.group(1)), int(m.group(2)), int(r.reads), '(split)' in str(getattr(r, 'description', ''))))
    return pd.DataFrame(ev, columns=['sample', 'target', 'length', 'pos', 'reads', 'split'])


def wquantile(x, w, q):
    o = np.argsort(x); x, w = np.asarray(x)[o], np.asarray(w, float)[o]
    c = np.cumsum(w) / w.sum()
    return float(x[min(len(x) - 1, np.searchsorted(c, q - 1e-12))])   # step quantile: an observed length


def per_sample(df, ev):
    out = []
    for (s, t), g in df.groupby(['sample', 'target']):
        n = g['reads'].sum()
        wt = g.loc[g['allele'] == 'WT', 'reads'].sum()
        has_del = g['allele'].str.contains(r'D\d+@', regex=True)
        fs = g.loc[g.get('frame', pd.Series('', index=g.index)) == 'frameshift', 'reads'].sum()
        e = ev[(ev['sample'] == s) & (ev['target'] == t)]
        row = dict(sample=s, target=t, reads=int(n), edited_pct=100 * (n - wt) / n, frameshift_pct=100 * fs / n,
                   reads_with_deletion_pct=100 * g.loc[has_del, 'reads'].sum() / n,
                   distinct_deletion_alleles_ge1pct=int(((g['reads'] / n >= 0.01) & has_del).sum()))
        if len(e):
            L, W = e['length'].to_numpy(), e['reads'].to_numpy()
            row.update(deletion_events=int(W.sum()), del_median_bp=wquantile(L, W, .5), del_q25_bp=wquantile(L, W, .25), del_q75_bp=wquantile(L, W, .75),
                       del_mean_bp=float(np.average(L, weights=W)), del_geomean_bp=float(np.exp(np.average(np.log(L), weights=W))),
                       del_max_bp=int(L.max()), large_del_pct_of_deletions=100 * W[L >= LARGE].sum() / W.sum(),
                       large_del_pct_of_reads=100 * W[L >= LARGE].sum() / n)
        out.append(row)
    return pd.DataFrame(out)


def binned(ev):
    lab = [f'{int(a) + 1}-{int(b)}' if np.isfinite(b) else f'>{int(a)}' for a, b in zip(BIN_EDGES[:-1], BIN_EDGES[1:])]
    ev = ev.copy(); ev['bin'] = pd.cut(ev['length'], BIN_EDGES, labels=lab, right=True)
    t = ev.pivot_table(index=['sample', 'target'], columns='bin', values='reads', aggfunc='sum', fill_value=0, observed=False)
    return (t.div(t.sum(axis=1), axis=0) * 100).round(3)


def sample_level_permutation(ev, groups, stat, n_perm=10000, seed=1):
    """permute group labels across SAMPLES (keeps each animal's reads together)"""
    rng = np.random.default_rng(seed)
    samples = sorted(groups); lab = np.array([groups[s] for s in samples])
    obs = stat(ev, dict(zip(samples, lab)))
    perm = [stat(ev, dict(zip(samples, rng.permutation(lab)))) for _ in range(n_perm)]
    return obs, (1 + np.sum(np.abs(perm) >= abs(obs))) / (n_perm + 1)


def group_tests(summ, ev, groups, out):
    lines = []
    summ = summ[summ['sample'].isin(groups)].copy(); summ['group'] = summ['sample'].map(groups)
    gnames = sorted(summ['group'].unique())
    lines.append(f'Groups: ' + ', '.join(f'{g} (n={sum(summ.group == g)} samples)' for g in gnames))
    if len(gnames) < 2: lines.append('Need ≥2 groups for comparisons.'); return lines
    metrics = ['edited_pct', 'frameshift_pct', 'reads_with_deletion_pct', 'del_median_bp', 'del_geomean_bp', 'large_del_pct_of_deletions', 'distinct_deletion_alleles_ge1pct']
    lines.append('\n1) Per-sample metrics (each animal = one value)')
    for mtr in metrics:
        if mtr not in summ: continue
        vals = [summ.loc[summ.group == g, mtr].dropna().to_numpy() for g in gnames]
        desc = ' | '.join(f'{g}: median {np.median(v):.2f} [{np.min(v):.2f}–{np.max(v):.2f}]' for g, v in zip(gnames, vals) if len(v))
        if any(len(v) < 2 for v in vals): lines.append(f'  {mtr:34s} {desc}  (too few samples for a test)'); continue
        if len(gnames) == 2: st, p = stats.mannwhitneyu(*vals, alternative='two-sided'); tn = 'Mann–Whitney U'
        else: st, p = stats.kruskal(*vals); tn = 'Kruskal–Wallis'
        lines.append(f'  {mtr:34s} {desc}  → {tn} p = {p:.3g}')
    # whole-distribution difference, each sample weighted equally, tested by permuting samples
    def w1(ev_, gmap):
        e = ev_[ev_['sample'].isin(gmap)].copy(); e['g'] = e['sample'].map(gmap)
        e['w'] = e['reads'] / e.groupby('sample')['reads'].transform('sum')       # each sample sums to 1
        a, b = [e[e.g == g] for g in gnames[:2]]
        if not len(a) or not len(b): return 0.0
        return stats.wasserstein_distance(np.log10(a['length']), np.log10(b['length']), a['w'], b['w'])
    def dmed(ev_, gmap):
        e = ev_[ev_['sample'].isin(gmap)]
        med = e.groupby('sample').apply(lambda d: wquantile(np.log2(d['length']), d['reads'], .5), include_groups=False)
        g = pd.Series(gmap); return med[g[g == gnames[0]].index.intersection(med.index)].mean() - med[g[g == gnames[1]].index.intersection(med.index)].mean()
    if len(gnames) == 2:
        lines.append(f'\n2) Whole-distribution comparison {gnames[0]} vs {gnames[1]} (sample-label permutation, 10,000 permutations)')
        o, p = sample_level_permutation(ev, groups, w1); lines.append(f'  Wasserstein distance of log10 deletion length = {o:.3f}  p = {p:.3g}')
        o, p = sample_level_permutation(ev, groups, dmed); lines.append(f'  Difference in mean per-sample median log2 length = {o:.3f} (×{2 ** o:.2f})  p = {p:.3g}')
        lines.append('  (with few animals per group the smallest attainable p is limited: e.g. 3 vs 3 → p ≥ 0.1)')
    # clustered regression (GEE): reads nested in samples
    try:
        import statsmodels.api as sm, statsmodels.formula.api as smf
        e = ev[ev['sample'].isin(groups)].copy(); e['group'] = pd.Categorical(e['sample'].map(groups), categories=gnames)
        rng = np.random.default_rng(2)
        # expand read counts, capped per sample so deep samples don't dominate
        parts = []
        for s, d in e.groupby('sample'):
            k = int(min(2000, d['reads'].sum())); idx = rng.choice(len(d), size=k, p=d['reads'] / d['reads'].sum()); parts.append(d.iloc[idx])
        x = pd.concat(parts); x['loglen'] = np.log2(x['length']); x['large'] = (x['length'] >= LARGE).astype(int)
        lines.append('\n3) GEE models, reads clustered within samples (exchangeable correlation, ≤2,000 reads per sample)')
        m1 = smf.gee('loglen ~ group', 'sample', x, cov_struct=sm.cov_struct.Exchangeable(), family=sm.families.Gaussian()).fit()
        for k in m1.params.index[1:]: lines.append(f'  log2(deletion length): {k} β = {m1.params[k]:.3f} (fold {2 ** m1.params[k]:.2f})  95% CI {2 ** m1.conf_int().loc[k, 0]:.2f}–{2 ** m1.conf_int().loc[k, 1]:.2f}  p = {m1.pvalues[k]:.3g}')
        if x['large'].nunique() > 1:
            m2 = smf.gee('large ~ group', 'sample', x, cov_struct=sm.cov_struct.Exchangeable(), family=sm.families.Binomial()).fit()
            for k in m2.params.index[1:]: lines.append(f'  large deletion (≥{LARGE} bp): {k} OR = {np.exp(m2.params[k]):.2f}  95% CI {np.exp(m2.conf_int().loc[k, 0]):.2f}–{np.exp(m2.conf_int().loc[k, 1]):.2f}  p = {m2.pvalues[k]:.3g}')
        lines.append('  GEE needs a reasonable number of samples (≥5–10 per group); with fewer, rely on 1) and 2).')
    except ImportError:
        lines.append('\n3) GEE models skipped (pip install statsmodels)')
    except Exception as ex:
        lines.append(f'\n3) GEE models could not be fitted: {ex}')
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('tsv', nargs='+', help='GenoView alleles TSV files (globs ok)')
    ap.add_argument('--groups', help='TSV: sample<TAB>group (sample = BAM name without .bam)')
    ap.add_argument('--target', help='only this target label (e.g. "cut 1")')
    ap.add_argument('--out', default='deletion_stats_out')
    a = ap.parse_args()
    paths = sorted({p for t in a.tsv for p in glob.glob(t)})
    df = read_alleles(paths)
    if a.target: df = df[df['target'] == a.target]
    os.makedirs(a.out, exist_ok=True)
    ev = deletion_events(df)
    summ = per_sample(df, ev); summ.to_csv(os.path.join(a.out, 'per_sample_summary.tsv'), sep='\t', index=False, float_format='%.3f')
    if len(ev): binned(ev).to_csv(os.path.join(a.out, 'binned_distribution.tsv'), sep='\t')
    rep = [f'Samples: {summ["sample"].nunique()}  ·  targets: {", ".join(map(str, summ["target"].unique()))}', '',
           summ.to_string(index=False, float_format=lambda v: f'{v:.2f}')]
    if a.groups:
        gm = pd.read_csv(a.groups, sep='\t', header=None, names=['sample', 'group'], comment='#')
        gm['sample'] = gm['sample'].str.replace(r'\.bam$', '', regex=True)
        rep += ['', *group_tests(summ, ev, dict(zip(gm['sample'], gm['group'].astype(str))), a.out)]
    txt = '\n'.join(rep); open(os.path.join(a.out, 'group_tests.txt'), 'w').write(txt + '\n'); print(txt)


if __name__ == '__main__':
    main()
