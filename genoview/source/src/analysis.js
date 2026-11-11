import { yieldNow } from './yield.js';
// GenoView analysis engine: depth, small indels, SVs (split reads, discordant pairs,
// soft clips, CIGAR), depth-based CNV segments and RNA-seq splice junctions.
// Works on any object exposing @gmod/bam's BamFile interface.

const F_PAIRED = 1, F_UNMAPPED = 4, F_MUNMAP = 8, F_REV = 16, F_MREV = 32,
  F_R1 = 64, F_SECONDARY = 256, F_QCFAIL = 512, F_DUP = 1024, F_SUPP = 2048;

export const DEFAULTS = {
  minMapq: 20,
  minIndel: 1,          // smallest CIGAR indel reported
  minSV: 50,            // CIGAR events >= this are SVs
  minSupport: 3,        // reads supporting a call (per call, pooled within a group)
  minAF: 0.1,           // small-indel allele fraction
  minClip: 20,          // soft-clip length counted as breakpoint evidence
  insertSize: 'auto',   // or {median, mad}
  insertSD: 5,          // MADs above median => discordant
  skipDuplicates: true,
  library: 'dna',       // 'dna' | 'rna'
  binSize: 'auto',
  cnvLossLog2: -0.6,
  cnvGainLog2: 0.45,
  cnvMinBins: 4,
  cnvBin: 250,
  windowSize: 250_000,
};

// ---------- helpers ----------
const sleep0 = () => yieldNow();

export function chromAlias(name, bamChroms) {
  if (bamChroms.has(name)) return name;
  const cands = [];
  if (name.startsWith('chr')) cands.push(name.slice(3)); else cands.push('chr' + name);
  if (/^(chr)?M(T)?$/.test(name)) cands.push('MT', 'chrM', 'M');
  for (const c of cands) if (bamChroms.has(c)) return c;
  return null;
}

function parseCigarString(s) {
  const out = []; const re = /(\d+)([MIDNSHP=X])/g; let m;
  while ((m = re.exec(s))) out.push([m[2], +m[1]]);
  return out;
}

function parseSA(sa) {
  return String(sa).split(';').filter(Boolean).map(x => {
    const [chr, pos, strand, cigar, mapq] = x.split(',');
    const ops = parseCigarString(cigar);
    let refLen = 0; for (const [op, l] of ops) if ('MDN=X'.includes(op)) refLen += l;
    const lc = ('SH'.includes(ops[0][0])) ? ops[0][1] : 0;
    const last = ops[ops.length - 1];
    const rc = ('SH'.includes(last[0])) ? last[1] : 0;
    return { chr, start: +pos - 1, end: +pos - 1 + refLen, rev: strand === '-', mapq: +mapq, lc, rc };
  });
}

function median(a) { if (!a.length) return 0; const s = Float64Array.from(a).sort(); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }

// ---------- sample model ----------
// sample = {id, name, group, bam (BamFile), isControl}
export async function openSampleHeader(bam) {
  const header = await bam.getHeader();
  const chroms = new Set((bam.indexToChr || []).map(c => c.refName));
  const sm = /@RG.*\tSM:([^\t\n]+)/.exec(header || '');
  return { chroms, refNames: bam.indexToChr || [], sampleName: sm ? sm[1] : null };
}

async function estimateInsert(sample, chr, start, end, opts) {
  const tl = [];
  const span = end - start;
  const probes = span <= 60_000 ? [[start, end]] :
    [0.15, 0.5, 0.85].map(f => { const c = start + Math.floor(span * f); return [c - 20_000, c + 20_000]; });
  for (const [s, e] of probes) {
    const recs = await sample.bam.getRecordsForRange(sample.bamChr, Math.max(0, s), e);
    for (const r of recs) {
      const f = r.flags;
      if (!(f & F_PAIRED) || (f & (F_UNMAPPED | F_MUNMAP | F_SECONDARY | F_SUPP))) continue;
      if (r.next_refid !== r.ref_id || (r.mq ?? 0) < opts.minMapq) continue;
      const t = r.template_length;
      if (t > 0 && !!(f & F_REV) !== !!(f & F_MREV)) tl.push(t);
      if (tl.length > 20000) break;
    }
  }
  if (tl.length < 100) return null;
  const med = median(tl);
  const mad = median(tl.map(x => Math.abs(x - med))) * 1.4826;
  return { median: med, mad: Math.max(mad, 10), n: tl.length };
}

// ---------- main scan ----------
export async function analyzeRegion({ samples, chr, start, end, options = {}, onProgress = () => {}, signal }) {
  const opts = { ...DEFAULTS, ...options };
  const span = end - start;
  const binSize = opts.binSize === 'auto' ? Math.max(1, Math.ceil(span / 1500)) : +opts.binSize;
  const nbins = Math.ceil(span / binSize);
  const groups = [...new Set(samples.map(s => s.group))];

  const perSample = [];
  const evidence = [];   // SV evidence
  const indels = new Map(); // key -> {type,chr,pos,len,seq,support:{group:n}}
  const junctions = new Map();
  const clips = new Map();  // `${pos}|${side}` -> {pos, side, support:{group:n}}

  let done = 0; const totalWork = samples.length * span;
  for (const s of samples) {
    if (signal?.aborted) throw new Error('cancelled');
    s.bamChr = chromAlias(chr, s.chroms);
    const depth = new Float64Array(nbins);
    const st = { id: s.id, name: s.name, group: s.group, depth, reads: 0, insert: null, bamChr: s.bamChr };
    perSample.push(st);
    if (!s.bamChr) { done += span; continue; }

    // insert size
    let ins = opts.insertSize !== 'auto' ? opts.insertSize : await estimateInsert(s, chr, start, end, opts);
    st.insert = ins;
    const maxInsert = ins ? ins.median + opts.insertSD * ins.mad : Infinity;
    const minInsert = ins ? Math.max(0, ins.median - opts.insertSD * ins.mad) : 0;

    for (let ws = start; ws < end; ws += opts.windowSize) {
      if (signal?.aborted) throw new Error('cancelled');
      const we = Math.min(end, ws + opts.windowSize);
      const recs = await s.bam.getRecordsForRange(s.bamChr, ws, we);
      for (const r of recs) {
        const rs = r.start;
        if (rs < ws && ws !== start) continue;   // counted in previous window
        const f = r.flags;
        if (f & (F_UNMAPPED | F_SECONDARY | F_QCFAIL)) continue;
        if (opts.skipDuplicates && (f & F_DUP)) continue;
        if ((r.mq ?? 0) < opts.minMapq) continue;
        st.reads++;
        const cig = r.NUMERIC_CIGAR;
        let pos = rs, qpos = 0;
        const nops = cig.length;
        for (let i = 0; i < nops; i++) {
          const op = cig[i] & 0xf, len = cig[i] >>> 4;
          if (op === 0 || op === 7 || op === 8) {        // M = X : depth
            addDepth(depth, pos, pos + len, start, end, binSize);
            pos += len; qpos += len;
          } else if (op === 2 || (op === 3 && opts.library === 'dna')) { // D (or N in DNA)
            if (len >= opts.minIndel && !(f & F_SUPP)) {
              if (len >= opts.minSV) evidence.push({ type: 'DEL', chr, pos1: pos, chr2: chr, pos2: pos + len, kind: 'cigar', group: s.group, precise: true });
              else addIndel(indels, 'DEL', chr, pos, len, '', s.group);
            }
            pos += len;
          } else if (op === 3) {                            // N in RNA = junction
            const key = pos + '-' + (pos + len);
            let j = junctions.get(key);
            if (!j) { j = { chr, start: pos, end: pos + len, strand: '.', support: {} }; junctions.set(key, j); }
            j.support[s.group] = (j.support[s.group] || 0) + 1;
            pos += len;
          } else if (op === 1) {                            // I
            if (len >= opts.minIndel && !(f & F_SUPP)) {
              if (len >= opts.minSV) evidence.push({ type: 'INS', chr, pos1: pos, chr2: chr, pos2: pos + 1, svlen: len, kind: 'cigar', group: s.group, precise: true });
              else addIndel(indels, 'INS', chr, pos, len, safeSeq(r, qpos, len), s.group);
            }
            qpos += len;
          } else if (op === 4) {                            // S
            if (len >= opts.minClip) {
              const side = i === 0 ? 'L' : 'R';
              const key = pos + '|' + side;
              let c = clips.get(key);
              if (!c) { c = { pos, side, support: {} }; clips.set(key, c); }
              c.support[s.group] = (c.support[s.group] || 0) + 1;
            }
            qpos += len;
          } else if (op === 5) { /* H */ if (len >= opts.minClip && !(f & F_SUPP)) { const side = i === 0 ? 'L' : 'R'; const key = pos + '|' + side; let c = clips.get(key); if (!c) { c = { pos, side, support: {} }; clips.set(key, c); } c.support[s.group] = (c.support[s.group] || 0) + 1; } }
        }
        const rend = pos;

        // split reads (primary only, to count each read once)
        if (!(f & F_SUPP)) {
          const sa = r.getTag ? r.getTag('SA') : undefined;
          if (sa) {
            const segs = parseSA(sa).filter(x => x.mapq >= Math.min(opts.minMapq, 1));
            if (segs.length) splitEvidence(evidence, chr, s, r, rs, rend, cig, segs[0], opts);
          }
        }

        // discordant pairs: count each pair once, from the leftmost mate
        if ((f & F_PAIRED) && !(f & (F_MUNMAP | F_SUPP))) {
          const mateChrIdx = r.next_refid;
          if (mateChrIdx !== r.ref_id) {
            if (mateChrIdx >= 0 && (f & F_R1)) {
              const mchr = s.refNames[mateChrIdx]?.refName;
              const rev = !!(f & F_REV);
              evidence.push({ type: 'BND', chr, pos1: rev ? rs : rend, chr2: mchr, pos2: r.next_pos, kind: 'disc', group: s.group, precise: false, strands: (rev ? '-' : '+') + (f & F_MREV ? '-' : '+') });
            }
          } else if (rs < r.next_pos || (rs === r.next_pos && (f & F_R1))) {
            const rev = !!(f & F_REV), mrev = !!(f & F_MREV);
            const tl = Math.abs(r.template_length) || (r.next_pos - rs + (rend - rs));
            let type = null;
            if (!rev && mrev) { if (tl > maxInsert && opts.library !== 'rna') type = 'DEL'; }
            else if (rev && !mrev) { if (tl > Math.max(maxInsert, (rend - rs) * 2)) type = 'DUP'; }
            else type = 'INV';
            if (type) evidence.push({ type, chr, pos1: rev ? rs : rend, chr2: chr, pos2: mrev ? r.next_pos + (rend - rs) : r.next_pos, kind: 'disc', group: s.group, precise: false, mateLeftFwd: !rev });
          }
        }
      }
      done += (we - ws);
      onProgress({ fraction: done / totalWork, sample: s.name, at: we });
      await sleep0();
    }
  }

  // ----- pool depth by group, normalise -----
  const groupDepth = {};
  for (const g of groups) groupDepth[g] = new Float64Array(nbins);
  for (const st of perSample) { const gd = groupDepth[st.group]; for (let i = 0; i < nbins; i++) gd[i] += st.depth[i]; }
  const pooled = new Float64Array(nbins);
  for (const g of groups) for (let i = 0; i < nbins; i++) pooled[i] += groupDepth[g][i];

  const depthAt = (g, pos) => { const b = Math.floor((pos - start) / binSize); const d = g === '*' ? pooled : groupDepth[g]; return (b >= 0 && b < nbins) ? d[b] : NaN; };
  const meanDepth = (g, a, b) => {
    const d = g === '*' ? pooled : groupDepth[g];
    let i0 = Math.max(0, Math.floor((a - start) / binSize)), i1 = Math.min(nbins - 1, Math.floor((b - 1 - start) / binSize));
    if (i1 < i0) return NaN; let s = 0; for (let i = i0; i <= i1; i++) s += d[i]; return s / (i1 - i0 + 1);
  };

  // ----- small indels -----
  const smallCalls = [];
  for (const v of indels.values()) {
    const tot = Object.values(v.support).reduce((a, b) => a + b, 0);
    if (tot < opts.minSupport) continue;
    const perGroup = {};
    let pass = false;
    for (const g of groups) {
      const n = v.support[g] || 0; const d = depthAt(g, v.pos);
      const af = d > 0 ? Math.min(1, n / Math.max(d, n)) : 0;
      perGroup[g] = { support: n, depth: d, af, gt: n === 0 ? '0/0' : af > 0.8 ? '1/1' : af >= opts.minAF ? '0/1' : './.' };
      if (n >= Math.min(opts.minSupport, 2) && af >= opts.minAF) pass = true;
    }
    if (!pass) continue;
    smallCalls.push({ type: v.type, chr, start: v.pos, end: v.type === 'DEL' ? v.pos + v.len : v.pos + 1, svlen: v.type === 'DEL' ? -v.len : v.len, seq: v.seq, support: tot, evidence: { cigar: tot }, perGroup, precise: true, small: true });
  }

  // ----- cluster SV evidence -----
  const tol = Math.max(...perSample.map(p => p.insert ? p.insert.median + 3 * p.insert.mad : 500), 500);
  const svCalls = clusterEvidence(evidence, tol).filter(c => c.support >= opts.minSupport);
  for (const c of svCalls) {
    c.perGroup = {};
    for (const g of groups) {
      const n = c.supportByGroup[g] || 0;
      const pg = { support: n };
      if (c.chr2 === c.chr && (c.type === 'DEL' || c.type === 'DUP' || c.type === 'INV') && c.end > c.start) {
        const fl = Math.max(binSize * 3, Math.min(5000, (c.end - c.start)));
        const inside = meanDepth(g, c.start, c.end);
        const flank = (meanDepth(g, c.start - fl, c.start) + meanDepth(g, c.end, c.end + fl)) / 2;
        pg.depthInside = inside; pg.depthFlank = flank; pg.ratio = flank > 0 ? inside / flank : NaN;
      }
      const bpd = depthAt(g, c.start);
      pg.af = bpd > 0 ? Math.min(1, n / (bpd + n * 0.5)) : 0;
      c.perGroup[g] = pg;
    }
    // overlapping soft-clip clusters confirm breakpoints
    c.clipSupport = clipNear(clips, c.start, 15) + (c.chr2 === c.chr ? clipNear(clips, c.end, 15) : 0);
  }

  let nid = 0; for (const c of [...svCalls, ...smallCalls]) c.id = 'call' + (++nid);
  // ----- depth-only CNV segments (per group vs median or control) -----
  const control = samples.find(s => s.isControl)?.group;
  const cnv = opts.library === 'rna' ? [] : depthSegments({ groupDepth, groups, control, nbins, binSize, start, chr, opts });
  // annotate CNVs that overlap a breakpoint call
  for (const seg of cnv) {
    seg.supportedBy = svCalls.find(c => c.chr2 === c.chr && c.type === seg.type && overlap(c.start, c.end, seg.start, seg.end) > 0.5 * Math.min(seg.end - seg.start, c.end - c.start))?.id || null;
  }

  // ----- junctions -----
  const juncs = [...junctions.values()].map(j => ({ ...j, total: Object.values(j.support).reduce((a, b) => a + b, 0) }))
    .filter(j => j.total >= Math.max(1, opts.minSupport - 1)).sort((a, b) => a.start - b.start);

  const calls = [...svCalls, ...smallCalls].sort((a, b) => a.start - b.start);
  for (const c of calls) c.score = scoreCall(c);
  return {
    chr, start, end, binSize, nbins, groups, control, perSample, groupDepth, pooled,
    calls, cnv, junctions: juncs, options: opts,
    medianDepth: Object.fromEntries(groups.map(g => [g, median(Array.from(groupDepth[g]).filter(x => x > 0))])),
  };
}

function overlap(a1, a2, b1, b2) { return Math.max(0, Math.min(a2, b2) - Math.max(a1, b1)); }

function safeSeq(r, qpos, len) { try { const s = r.seq; return s ? s.substr(qpos, len) : ''; } catch { return ''; } }

function addDepth(depth, a, b, start, end, bin) {
  a = Math.max(a, start); b = Math.min(b, end); if (b <= a) return;
  let i = Math.floor((a - start) / bin);
  while (a < b) {
    const binEnd = start + (i + 1) * bin;
    const e = Math.min(b, binEnd);
    depth[i] += (e - a) / bin;
    a = e; i++;
  }
}

function addIndel(map, type, chr, pos, len, seq, group) {
  const key = type + ':' + pos + ':' + len;
  let v = map.get(key);
  if (!v) { v = { type, chr, pos, len, seq, support: {} }; map.set(key, v); }
  if (!v.seq && seq) v.seq = seq;
  v.support[group] = (v.support[group] || 0) + 1;
}

function clipNear(clips, pos, w) {
  let n = 0;
  for (let p = pos - w; p <= pos + w; p++) for (const side of ['L', 'R']) { const c = clips.get(p + '|' + side); if (c) n += Object.values(c.support).reduce((a, b) => a + b, 0); }
  return n;
}

function splitEvidence(evidence, chr, s, r, rs, rend, cig, sa, opts) {
  const n = cig.length;
  const op0 = cig[0] & 0xf, opN = cig[n - 1] & 0xf;
  const lc = (op0 === 4 || op0 === 5) ? cig[0] >>> 4 : 0;
  const rc = (opN === 4 || opN === 5) ? cig[n - 1] >>> 4 : 0;
  const pSide = rc >= lc ? 'R' : 'L';
  const pBp = pSide === 'R' ? rend : rs;
  const rev = !!(r.flags & F_REV);
  const sSide = sa.rc > sa.lc ? 'R' : 'L';
  const sBp = sSide === 'R' ? sa.end : sa.start;
  const base = { kind: 'split', group: s.group, precise: true };
  if (sa.chr !== s.bamChr) {
    evidence.push({ ...base, type: 'BND', chr, pos1: pBp, chr2: sa.chr, pos2: sBp, strands: (pSide === 'R' ? '+' : '-') + (sSide === 'L' ? '+' : '-') });
    return;
  }
  const [a, b] = pBp <= sBp ? [pBp, sBp] : [sBp, pBp];
  if (sa.rev !== rev) { evidence.push({ ...base, type: 'INV', chr, pos1: a, chr2: chr, pos2: b }); return; }
  const leftSide = pBp <= sBp ? pSide : sSide;   // clip side of the leftmost segment
  if (b - a < opts.minSV) {
    if (pSide !== sSide) evidence.push({ ...base, type: 'INS', chr, pos1: a, chr2: chr, pos2: a + 1 });
    return;
  }
  evidence.push({ ...base, type: leftSide === 'R' ? 'DEL' : 'DUP', chr, pos1: a, chr2: chr, pos2: b });
}

function clusterEvidence(ev, tol) {
  const byKey = new Map();
  for (const e of ev) {
    const k = e.type + '|' + e.chr + '|' + e.chr2;
    if (!byKey.has(k)) byKey.set(k, []);
    byKey.get(k).push(e);
  }
  const calls = [];
  for (const list of byKey.values()) {
    list.sort((a, b) => a.pos1 - b.pos1);
    const clusters = [];
    for (const e of list) {
      let best = null;
      for (let i = clusters.length - 1; i >= 0; i--) {
        const c = clusters[i];
        if (e.pos1 - c.max1 > tol) break;
        if (Math.abs(e.pos2 - c.med2()) <= tol) { best = c; break; }
      }
      if (!best) { best = { items: [], max1: -Infinity, med2() { return median(this.items.map(x => x.pos2)); } }; clusters.push(best); }
      best.items.push(e); best.max1 = Math.max(best.max1, e.pos1);
    }
    for (const c of clusters) calls.push(summariseCluster(c.items));
  }
  return calls;
}

function mode(arr) { const m = new Map(); let best = arr[0], bc = 0; for (const x of arr) { const c = (m.get(x) || 0) + 1; m.set(x, c); if (c > bc) { bc = c; best = x; } } return best; }

function summariseCluster(items) {
  const e0 = items[0];
  const precise = items.filter(x => x.precise);
  const disc = items.filter(x => x.kind === 'disc');
  let start, end;
  if (precise.length) { start = mode(precise.map(x => x.pos1)); end = mode(precise.map(x => x.pos2)); }
  else {
    // discordant only: breakpoint lies after the furthest left-read end and before the nearest right-read start
    start = Math.max(...disc.map(x => x.pos1)); end = Math.min(...disc.map(x => x.pos2));
    if (end < start) [start, end] = [Math.min(...disc.map(x => x.pos1)), Math.max(...disc.map(x => x.pos2))];
  }
  const supportByGroup = {}, evidence = { split: 0, disc: 0, cigar: 0 };
  // a read pair may give split + disc; count templates once is not possible without names, so report both
  for (const x of items) { supportByGroup[x.group] = (supportByGroup[x.group] || 0) + 1; evidence[x.kind]++; }
  const support = items.length;
  const type = e0.type === 'BND' ? 'TRA' : e0.type;
  const svlen = e0.chr2 !== e0.chr ? null : type === 'INS' ? (e0.svlen || null) : type === 'DEL' ? -(end - start) : (end - start);
  return { type, chr: e0.chr, start: Math.min(start, e0.chr === e0.chr2 ? end : start), end: e0.chr === e0.chr2 ? Math.max(start, end) : start + 1, chr2: e0.chr2, pos2: end, svlen, support, supportByGroup, evidence, precise: precise.length > 0, strands: e0.strands || null };
}

function depthSegments({ groupDepth, groups, control, nbins, binSize, start, chr, opts }) {
  // re-bin to >= cnvBin bp so Poisson / read-boundary noise does not create micro-segments
  const f = Math.max(1, Math.ceil((opts.cnvBin || 250) / binSize));
  const nb = Math.floor(nbins / f), bs = binSize * f;
  if (nb < opts.cnvMinBins) return [];
  const rebin = d => { const o = new Float64Array(nb); for (let i = 0; i < nb; i++) { let a = 0; for (let k = 0; k < f; k++) a += d[i * f + k]; o[i] = a / f; } return o; };
  const segs = [];
  const ctrl = control ? rebin(groupDepth[control]) : null;
  const cmed = ctrl ? median(Array.from(ctrl).filter(x => x > 0)) : 0;
  for (const g of groups) {
    if (g === control) continue;
    const d = rebin(groupDepth[g]);
    const med = median(Array.from(d).filter(x => x > 0));
    if (!(med > 5)) continue;
    const lr = new Float64Array(nb);
    for (let i = 0; i < nb; i++) {
      const denom = ctrl ? (ctrl[i] / cmed) * med : med;
      lr[i] = denom > 0 ? Math.log2((d[i] + 0.5) / (denom + 0.5)) : 0;
    }
    const sm = new Float64Array(nb);
    for (let i = 0; i < nb; i++) { const w = [lr[Math.max(0, i - 1)], lr[i], lr[Math.min(nb - 1, i + 1)]].sort((a, b) => a - b); sm[i] = w[1]; }
    let i = 0;
    while (i < nb) {
      const t = sm[i] <= opts.cnvLossLog2 ? 'DEL' : sm[i] >= opts.cnvGainLog2 ? 'DUP' : null;
      if (!t) { i++; continue; }
      let j = i; let sum = 0;
      while (j < nb && (t === 'DEL' ? sm[j] <= opts.cnvLossLog2 * 0.7 : sm[j] >= opts.cnvGainLog2 * 0.7)) { sum += lr[j]; j++; }
      const mean = sum / (j - i);
      if (j - i >= opts.cnvMinBins && (t === 'DEL' ? mean <= opts.cnvLossLog2 : mean >= opts.cnvGainLog2))
        segs.push({ type: t, chr, group: g, start: start + i * bs, end: start + j * bs, log2: mean, bins: j - i, cn: +(2 * Math.pow(2, mean)).toFixed(2) });
      i = j;
    }
  }
  return segs;
}

function scoreCall(c) {
  const e = c.evidence || {};
  let q = (e.split || 0) * 3 + (e.disc || 0) * 2 + (e.cigar || 0) * 3 + (c.clipSupport || 0) * 0.5;
  for (const pg of Object.values(c.perGroup || {})) {
    if (pg.ratio != null && !isNaN(pg.ratio)) {
      if (c.type === 'DEL' && pg.ratio < 0.8) q += 10 * (1 - pg.ratio);
      if (c.type === 'DUP' && pg.ratio > 1.2) q += 10 * (pg.ratio - 1);
    }
  }
  return Math.round(Math.min(999, q));
}

// ---------- exports ----------
export function callsToVCF(res, sampleOrder) {
  const groups = sampleOrder || res.groups;
  const lines = ['##fileformat=VCFv4.2', `##source=GenoView`, `##region=${res.chr}:${res.start + 1}-${res.end}`,
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of variant">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
    '##INFO=<ID=CHR2,Number=1,Type=String,Description="Partner chromosome">',
    '##INFO=<ID=POS2,Number=1,Type=Integer,Description="Partner position">',
    '##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">',
    '##INFO=<ID=PE,Number=1,Type=Integer,Description="Discordant pair support">',
    '##INFO=<ID=CG,Number=1,Type=Integer,Description="CIGAR support">',
    '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise breakpoints">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (rough)">',
    '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Supporting reads">',
    '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Approx. allele fraction">',
    '##FORMAT=<ID=DR,Number=1,Type=Float,Description="Depth ratio inside/flank">',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + groups.join('\t')];
  for (const c of res.calls) {
    const info = [`SVTYPE=${c.type === 'TRA' ? 'BND' : c.type}`, `END=${c.end}`];
    if (c.svlen != null) info.push(`SVLEN=${c.svlen}`);
    if (c.type === 'TRA') info.push(`CHR2=${c.chr2}`, `POS2=${c.pos2 + 1}`);
    const e = c.evidence || {}; info.push(`SR=${e.split || 0}`, `PE=${e.disc || 0}`, `CG=${e.cigar || 0}`);
    if (!c.precise) info.push('IMPRECISE');
    const alt = c.small && c.type === 'INS' && c.seq ? 'N' + c.seq : `<${c.type === 'TRA' ? 'BND' : c.type}>`;
    const fmt = groups.map(g => { const p = c.perGroup[g] || {}; const af = p.af ?? 0; const gt = p.gt || (p.support ? (af > 0.8 ? '1/1' : '0/1') : '0/0'); return `${gt}:${p.support || 0}:${af.toFixed(2)}:${p.ratio != null && !isNaN(p.ratio) ? p.ratio.toFixed(2) : '.'}`; });
    lines.push([c.chr, c.start + 1, c.id, 'N', alt, c.score, 'PASS', info.join(';'), 'GT:AD:AF:DR', ...fmt].join('\t'));
  }
  return lines.join('\n') + '\n';
}

export function callsToBED(res) {
  return res.calls.filter(c => c.chr2 === c.chr || c.chr2 == null).map(c => `${c.chr}\t${c.start}\t${Math.max(c.end, c.start + 1)}\t${c.type}_${c.svlen ?? ''}_${c.id}\t${Math.min(1000, c.score)}`).join('\n') + '\n';
}

export function callsToBEDPE(res) {
  return res.calls.filter(c => !c.small).map(c => `${c.chr}\t${c.start}\t${c.start + 1}\t${c.chr2 || c.chr}\t${(c.type === 'TRA' ? c.pos2 : c.end)}\t${(c.type === 'TRA' ? c.pos2 : c.end) + 1}\t${c.type}_${c.id}\t${c.score}`).join('\n') + '\n';
}

export function callsToTSV(res) {
  const g = res.groups;
  const head = ['id', 'type', 'chr', 'start', 'end', 'chr2', 'pos2', 'svlen', 'score', 'precise', 'split', 'discordant', 'cigar', 'clips', ...g.flatMap(x => [`${x}_support`, `${x}_af`, `${x}_depth_ratio`])];
  const rows = res.calls.map(c => [c.id, c.type, c.chr, c.start + 1, c.end, c.chr2 || c.chr, c.pos2 != null ? c.pos2 + 1 : '', c.svlen ?? '', c.score, c.precise, c.evidence?.split || 0, c.evidence?.disc || 0, c.evidence?.cigar || 0, c.clipSupport || 0,
    ...g.flatMap(x => { const p = c.perGroup[x] || {}; return [p.support || 0, (p.af ?? 0).toFixed(3), p.ratio != null && !isNaN(p.ratio) ? p.ratio.toFixed(3) : '']; })]);
  return [head, ...rows].map(r => r.join('\t')).join('\n') + '\n';
}

export function depthToBedGraph(res, group) {
  const d = group === '*' ? res.pooled : res.groupDepth[group];
  const out = [`track type=bedGraph name="${group === '*' ? 'pooled' : group} depth"`];
  for (let i = 0; i < res.nbins; i++) { const s = res.start + i * res.binSize; out.push(`${res.chr}\t${s}\t${Math.min(res.end, s + res.binSize)}\t${d[i].toFixed(2)}`); }
  return out.join('\n') + '\n';
}

export function junctionsToBED(res) {
  return res.junctions.map((j, i) => `${j.chr}\t${j.start}\t${j.end}\tJ${i + 1}_${j.total}\t${j.total}\t${j.strand === '.' ? '+' : j.strand}`).join('\n') + '\n';
}
