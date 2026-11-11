// CRISPR editing outcome analysis for individual crispants (no control needed).
// Every read that spans the cut-site window is classified into an allele; indel sizes are
// summarised as an empirical distribution standardised to % of reads.

const F_UNMAPPED = 4, F_REV = 16, F_SECONDARY = 256, F_QCFAIL = 512, F_DUP = 1024, F_SUPP = 2048;
const COMP = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
export const revcomp = s => s.toUpperCase().split('').reverse().map(c => COMP[c] || 'N').join('');

export const NUCLEASES = {
  SpCas9: { pam: 'NGG', side: 3, cut: -3, label: 'SpCas9 (NGG)' },
  'SpCas9-NG': { pam: 'NG', side: 3, cut: -3, label: 'SpCas9-NG (NG)' },
  SpRY: { pam: 'NRN', side: 3, cut: -3, label: 'SpRY (NRN)' },
  SaCas9: { pam: 'NNGRRT', side: 3, cut: -3, label: 'SaCas9 (NNGRRT)' },
  Cas12a: { pam: 'TTTV', side: 5, cut: 18, label: 'Cas12a / Cpf1 (TTTV, 5′ PAM)' },
};
const IUPAC = { A: 'A', C: 'C', G: 'G', T: 'T', N: '[ACGT]', R: '[AG]', Y: '[CT]', S: '[GC]', W: '[AT]', K: '[GT]', M: '[AC]', B: '[CGT]', D: '[AGT]', H: '[ACT]', V: '[ACG]' };
const pamRe = p => new RegExp('^' + p.split('').map(c => IUPAC[c] || c).join('') + '$');

/** Find guide sites in seq (0-based offset `offset`). Returns [{strand, start, end, cut, pam}] */
export function findGuides(seq, offset, guide, nucKey = 'SpCas9') {
  const nuc = NUCLEASES[nucKey]; const g = guide.toUpperCase().replace(/[^ACGTN]/g, '');
  const S = seq.toUpperCase(); const out = []; const P = nuc.pam.length;
  const re = pamRe(nuc.pam), reRc = pamRe(revcomp(nuc.pam).replace(/[^ACGTN]/g, m => m));
  const rcPamPattern = new RegExp('^' + revcomp(nuc.pam.replace(/./g, c => c)).split('').map(c => IUPAC[c] || c).join('') + '$');
  // IUPAC revcomp
  const IU_RC = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N', R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K', B: 'V', V: 'B', D: 'H', H: 'D' };
  const pamRc = nuc.pam.split('').reverse().map(c => IU_RC[c]).join('');
  const reMinus = pamRe(pamRc);
  const L = g.length; if (L < 15) return out;
  const gRc = revcomp(g);
  let i = -1;
  while ((i = S.indexOf(g, i + 1)) >= 0) {           // + strand
    if (nuc.side === 3) { const pam = S.substr(i + L, P); if (re.test(pam)) out.push({ strand: '+', start: i + offset, end: i + L + offset, pam, cut: i + L + nuc.cut + offset }); }
    else { const pam = S.substr(i - P, P); if (i >= P && re.test(pam)) out.push({ strand: '+', start: i + offset, end: i + L + offset, pam, cut: i + nuc.cut + offset }); }
  }
  i = -1;
  while ((i = S.indexOf(gRc, i + 1)) >= 0) {         // − strand
    if (nuc.side === 3) { const pamF = S.substr(i - P, P); if (i >= P && reMinus.test(pamF)) out.push({ strand: '-', start: i + offset, end: i + L + offset, pam: revcomp(pamF), cut: i - nuc.cut + offset }); }
    else { const pamF = S.substr(i + L, P); if (reMinus.test(pamF)) out.push({ strand: '-', start: i + offset, end: i + L + offset, pam: revcomp(pamF), cut: i + L - nuc.cut + offset }); }
  }
  void reRc; void rcPamPattern;
  return out;
}

function parseSA(sa) {
  return String(sa).split(';').filter(Boolean).map(x => {
    const [chr, pos, strand, cigar, mapq] = x.split(',');
    const ops = [...cigar.matchAll(/(\d+)([MIDNSHP=X])/g)].map(m => [m[2], +m[1]]);
    let refLen = 0; for (const [op, l] of ops) if ('MDN=X'.includes(op)) refLen += l;
    return { chr, start: +pos - 1, end: +pos - 1 + refLen, rev: strand === '-', mapq: +mapq };
  });
}

/** Stage 1 — read the locus once. Keeps, per read, its aligned span (incl. split parts) and all indel events. */
export async function collectReads(sample, { regStart, regEnd, minMapq = 20, skipDup = false }) {
  const bamChr = sample.bamChr;
  const recs = await sample.bam.getRecordsForRange(bamChr, Math.max(0, regStart), regEnd);
  const raw = { name: sample.name, regStart: Math.max(0, regStart), regEnd, total: 0, lowMapq: 0, reads: [] };
  const seen = new Set();
  for (const r of recs) {
    const f = r.flags;
    if (f & (F_UNMAPPED | F_SECONDARY | F_SUPP | F_QCFAIL)) continue;
    if (skipDup && (f & F_DUP)) continue;
    const nm = r.name; if (seen.has(nm)) continue; seen.add(nm);
    raw.total++;
    if ((r.mq ?? 0) < minMapq) { raw.lowMapq++; continue; }
    const cig = r.NUMERIC_CIGAR; let pos = r.start, q = 0; const evs = [];
    let lclip = 0, rclip = 0; let seq = null;
    for (let i = 0; i < cig.length; i++) {
      const op = cig[i] & 0xf, len = cig[i] >>> 4;
      if (op === 0 || op === 7 || op === 8) { pos += len; q += len; }
      else if (op === 2 || op === 3) { evs.push({ kind: 'D', pos, len }); pos += len; }
      else if (op === 1) {
        let sq = '';
        if (pos >= regStart && pos <= regEnd) { try { seq = seq ?? r.seq; sq = seq ? seq.substr(q, len) : ''; } catch {} }
        evs.push({ kind: 'I', pos, len, seq: sq }); q += len;
      }
      else if (op === 4) { if (i === 0) lclip = len; else rclip = len; q += len; }
      else if (op === 5) { if (i === 0) lclip = len; else rclip = len; }
    }
    let s0 = r.start, s1 = pos, complex = null;
    const sa = r.getTag?.('SA');
    if (sa) {
      const rev = !!(f & F_REV);
      for (const x of parseSA(sa)) {
        if (x.mapq < minMapq) continue;
        if (x.chr !== bamChr) { complex = complex || 'translocation'; continue; }
        if (x.rev !== rev) { if (x.start < regEnd && x.end > regStart) complex = complex || 'inversion'; continue; }
        if (x.start >= s1 - 20 && x.start - s1 > 0) { evs.push({ kind: 'D', pos: s1, len: x.start - s1, split: true }); s1 = x.end; }
        else if (x.end <= s0 + 20 && s0 - x.end > 0) { evs.unshift({ kind: 'D', pos: x.end, len: s0 - x.end, split: true }); s0 = x.start; }
      }
    }
    // only keep events near the locus (memory)
    const keep = evs.filter(e => e.pos < regEnd && e.pos + (e.kind === 'D' ? e.len : 0) > regStart);
    raw.reads.push({ s0, s1, evs: keep, complex, lclip: lclip >= 15 ? r.start : null, rclip: rclip >= 15 ? pos : null });
  }
  return raw;
}

export const eventKey = e => `${e.kind}:${e.pos}:${e.len}${e.kind === 'I' && e.seq ? ':' + e.seq : ''}`;

/** parse "1-2, 9, >5000, <3" into a predicate on indel length */
export function parseSizeFilter(txt) {
  const parts = String(txt || '').split(/[\s,;]+/).filter(Boolean);
  const tests = [];
  for (const p of parts) {
    let m;
    if ((m = /^(\d+)\s*-\s*(\d+)$/.exec(p))) tests.push(n => n >= +m[1] && n <= +m[2]);
    else if ((m = /^(>=|<=|>|<)(\d+)$/.exec(p))) { const v = +m[2]; tests.push({ '>': n => n > v, '<': n => n < v, '>=': n => n >= v, '<=': n => n <= v }[m[1]]); }
    else if (/^\d+$/.test(p)) tests.push(n => n === +p);
  }
  return tests.length ? n => tests.some(t => t(n)) : null;
}

/**
 * Stage 2 — classify reads for one target. target = { id, label, cuts:[...], areas:[[s,e],...], origin }
 * o: { minIndel, hpFilter, mergeTol, largeDel, ignoreSize (fn|null), mask:Set, ref }
 */
export function classify(raw, target, o) {
  const { cuts, areas, origin } = target;
  const largeDel = o.largeDel || 50;
  const refAt = p => (o.ref && p >= o.ref.start && p < o.ref.start + o.ref.seq.length) ? o.ref.seq[p - o.ref.start].toUpperCase() : 'N';
  const runLen = (p, b) => { let a = p, z = p; while (refAt(a - 1) === b && p - a < 50) a--; while (refAt(z) === b && z - p < 50) z++; return z - a; };
  const isHp = ev => {
    if (!o.hpFilter || ev.len > 2 || !o.ref) return false;
    if (ev.kind === 'D') { const b = refAt(ev.pos); for (let k = 1; k < ev.len; k++) if (refAt(ev.pos + k) !== b) return false; return runLen(ev.pos, b) >= 4; }
    const b = ev.seq?.[0]; if (!b || [...ev.seq].some(c => c !== b)) return false;
    return Math.max(runLen(ev.pos, b), runLen(ev.pos - 1, b)) >= 3 && (refAt(ev.pos) === b || refAt(ev.pos - 1) === b);
  };
  const touchesCut = e => cuts.some(c => e.kind === 'D' ? (e.pos <= c + 1 && e.pos + e.len >= c - 1) : Math.abs(e.pos - c) <= 1);
  const inArea = e => areas.some(([a, b]) => e.kind === 'D' ? (e.pos < b && e.pos + e.len > a) : (e.pos >= a && e.pos <= b));
  const lo = Math.min(...areas.map(a => a[0])), hi = Math.max(...areas.map(a => a[1]));
  const regStart = raw.regStart, regLen = raw.regEnd - raw.regStart;
  const res = { name: raw.name, target, origin, cuts, areas, total: raw.total, lowMapq: raw.lowMapq, spanning: 0, notSpanning: 0, unresolved: 0,
    filteredHp: 0, filteredNoise: 0, filteredSize: 0, masked: 0, regStart, delCov: new Float64Array(regLen), insAt: new Float64Array(regLen) };
  const alleles = new Map();
  const add = (key, info) => { let a = alleles.get(key); if (!a) { a = { key, count: 0, ...info }; alleles.set(key, a); } a.count++; };
  for (const rd of raw.reads) {
    const spans = areas.every(([a, b]) => rd.s0 <= a && rd.s1 >= b);
    if (!spans) {
      const clip = (rd.rclip != null && cuts.some(c => Math.abs(rd.rclip - c) <= 10)) || (rd.lclip != null && cuts.some(c => Math.abs(rd.lclip - c) <= 10));
      if (clip && !rd.complex && rd.s0 <= hi && rd.s1 >= lo) { res.spanning++; res.unresolved++; add('unresolved', { label: 'Soft-clipped at cut (unaligned insert?)', net: null, cls: 'unresolved', events: [] }); }
      else res.notSpanning++;
      continue;
    }
    res.spanning++;
    const kept = [];
    for (const e of rd.evs) {
      if (!inArea(e)) continue;
      if (o.kinds === 'D' && e.kind !== 'D') { res.ignoredIns = (res.ignoredIns || 0) + 1; continue; }
      if ((o.mask && o.mask.has(eventKey(e))) || (o.ruleMask && o.ruleMask(e))) { res.masked++; continue; }
      if (o.ignoreSize && o.ignoreSize(e.len)) { res.filteredSize++; continue; }
      if (isHp(e)) { res.filteredHp++; continue; }
      if (e.len < (o.minIndel || 1) && !touchesCut(e)) { res.filteredNoise++; continue; }
      kept.push(e);
    }
    let net = 0, maxDel = 0;
    for (const e of kept) {
      if (e.kind === 'D') { net -= e.len; maxDel = Math.max(maxDel, e.len); const a = Math.max(0, e.pos - regStart), b = Math.min(regLen, e.pos + e.len - regStart); for (let p = a; p < b; p++) res.delCov[p]++; }
      else { net += e.len; const ix = e.pos - regStart; if (ix >= 0 && ix < regLen) res.insAt[ix]++; }
    }
    let cls;
    if (rd.complex) cls = 'complex';
    else if (!kept.length) cls = 'wt';
    else if (maxDel >= largeDel) cls = 'large';
    else if (kept.every(e => e.kind === 'D')) cls = 'del';
    else if (kept.every(e => e.kind === 'I')) cls = 'ins';
    else cls = 'mixed';
    const key = rd.complex ? rd.complex : kept.map(e => e.kind === 'D' ? `D${e.len}@${e.pos - origin}` : `I${e.len}@${e.pos - origin}:${e.seq}`).join('|') || 'WT';
    add(key, { net: rd.complex ? null : net, cls, events: kept.map(e => ({ ...e, rel: e.pos - origin })), complex: rd.complex });
  }
  res.alleles = alleles;
  finalize(res, o);
  return res;
}

/**
 * Whole-sample indel spectrum: every indel event in any counting area (all cuts / ROIs pooled),
 * from reads spanning at least one area. Each read contributes each distinct event once.
 */
export function sampleSpectrum(raw, areas, cuts, o) {
  const refAt = p => (o.ref && p >= o.ref.start && p < o.ref.start + o.ref.seq.length) ? o.ref.seq[p - o.ref.start].toUpperCase() : 'N';
  const runLen = (p, b) => { let a = p, z = p; while (refAt(a - 1) === b && p - a < 50) a--; while (refAt(z) === b && z - p < 50) z++; return z - a; };
  const isHp = ev => {
    if (!o.hpFilter || ev.len > 2 || !o.ref) return false;
    if (ev.kind === 'D') { const b = refAt(ev.pos); for (let k = 1; k < ev.len; k++) if (refAt(ev.pos + k) !== b) return false; return runLen(ev.pos, b) >= 4; }
    const b = ev.seq?.[0]; if (!b || [...ev.seq].some(c => c !== b)) return false;
    return Math.max(runLen(ev.pos, b), runLen(ev.pos - 1, b)) >= 3 && (refAt(ev.pos) === b || refAt(ev.pos - 1) === b);
  };
  const touchesCut = e => cuts.some(c => e.kind === 'D' ? (e.pos <= c + 1 && e.pos + e.len >= c - 1) : Math.abs(e.pos - c) <= 1);
  const hit = (e, [a, b]) => e.kind === 'D' ? (e.pos < b && e.pos + e.len > a) : (e.pos >= a && e.pos <= b);
  const del = new Map(), ins = new Map();
  let n = 0, edited = 0, events = 0;
  for (const rd of raw.reads) {
    const spanned = areas.filter(([a, b]) => rd.s0 <= a && rd.s1 >= b);
    if (!spanned.length) continue;
    n++;
    const seen = new Set(); let any = false;
    for (const e of rd.evs) {
      if (!spanned.some(ar => hit(e, ar))) continue;
      if (o.kinds === 'D' && e.kind !== 'D') continue;
      if (o.mask && o.mask.has(eventKey(e))) continue;
      if (o.ruleMask && o.ruleMask(e)) continue;
      if (o.ignoreSize && o.ignoreSize(e.len)) continue;
      if (isHp(e)) continue;
      if (e.len < (o.minIndel || 1) && !touchesCut(e)) continue;
      const k = eventKey(e); if (seen.has(k)) continue; seen.add(k);
      const m = e.kind === 'D' ? del : ins; m.set(e.len, (m.get(e.len) || 0) + 1);
      events++; any = true;
    }
    if (any) edited++;
  }
  const delOnly = o.kinds === 'D';
  return { name: raw.name, unitLabel: 'reads', total: n, reads: n, edited, events, cuts: cuts.length, areas: areas.length, delOnly,
    series: delOnly ? [{ key: 'del', label: 'Deletions', color: '#e34948', map: del }] : [{ key: 'del', label: 'Deletions', color: '#e34948', map: del }, { key: 'ins', label: 'Insertions', color: '#1baf7a', map: ins }] };
}

/** Backwards-compatible single-cut helper (used by tests). */
export async function analyzeEdits(sample, o) {
  const raw = await collectReads(sample, { regStart: o.cut - o.flank, regEnd: o.cut + o.flank, minMapq: o.minMapq, skipDup: o.skipDup });
  const W = o.window ?? 5;
  return classify(raw, { id: 's1', label: 'cut', cuts: [o.cut], areas: [[o.cut - W, o.cut + W]], origin: o.cut }, o);
}

function mergeNear(list, tol) {
  // fold minor single-event alleles into a major one with the same kind and ~same breakpoints (long-read jitter)
  const out = [];
  for (const a of list) {
    const e = a.events.length === 1 ? a.events[0] : null;
    const tgt = e && out.find(b => b.events.length === 1 && b.events[0].kind === e.kind && b.cls !== 'wt' &&
      Math.abs(b.events[0].len - e.len) <= Math.max(tol, Math.round(0.01 * e.len)) && Math.abs(b.events[0].pos - e.pos) <= Math.max(tol, Math.round(0.01 * e.len)) &&
      ((b.net % 3 === 0) === (a.net % 3 === 0) || e.len > 30) && b.count >= a.count * 3);
    if (tgt) { tgt.count += a.count; tgt.merged = (tgt.merged || 0) + a.count; } else out.push(a);
  }
  return out;
}

function finalize(res, o) {
  const n = res.spanning || 1;
  let list = [...res.alleles.values()].sort((a, b) => b.count - a.count);
  if (o.mergeTol > 0) list = mergeNear(list, o.mergeTol).sort((a, b) => b.count - a.count);
  for (const a of list) {
    a.pct = 100 * a.count / n;
    a.frame = a.net == null ? '?' : a.cls === 'wt' ? 'wt' : (a.net % 3 === 0 ? 'in-frame' : 'frameshift');
    a.label = a.label || describe(a, res.cuts);
    a.dropout = res.cuts && res.cuts.length > 1 && a.events.some(e => e.kind === 'D' && res.cuts.filter(c => e.pos <= c && e.pos + e.len >= c).length > 1);
  }
  res.alleleList = list;
  const sum = pred => list.filter(pred).reduce((s, a) => s + a.count, 0);
  const wt = sum(a => a.cls === 'wt');
  res.stats = {
    spanning: res.spanning, wtPct: 100 * wt / n, editedPct: 100 * (n - wt) / n,
    frameshiftPct: 100 * sum(a => a.frame === 'frameshift') / n,
    inframePct: 100 * sum(a => a.frame === 'in-frame') / n,
    largePct: 100 * sum(a => a.cls === 'large') / n,
    complexPct: 100 * sum(a => a.cls === 'complex' || a.cls === 'unresolved') / n,
    dropoutPct: 100 * sum(a => a.dropout) / n,
    readsWithDelPct: 100 * sum(a => (a.events || []).some(e => e.kind === 'D')) / n,
    insPct: 100 * sum(a => a.cls === 'ins') / n, delPct: 100 * sum(a => a.cls === 'del' || a.cls === 'large') / n,
    alleles1: list.filter(a => a.cls !== 'wt' && a.pct >= 1).length,
    alleles5: list.filter(a => a.cls !== 'wt' && a.pct >= 5).length,
  };
  const dist = new Map();
  for (const a of list) if (a.net != null) dist.set(a.net, (dist.get(a.net) || 0) + a.count);
  res.netDist = [...dist.entries()].sort((a, b) => a[0] - b[0]).map(([size, count]) => ({ size, count, pct: 100 * count / n, wt: size === 0 ? (list.find(a => a.cls === 'wt')?.count || 0) : 0 }));
  const ed = list.filter(a => a.cls !== 'wt' && a.net != null);
  const tot = ed.reduce((s, a) => s + a.count, 0);
  res.stats.meanNet = tot ? ed.reduce((s, a) => s + a.net * a.count, 0) / tot : 0;
  for (const arr of [res.delCov, res.insAt]) for (let i = 0; i < arr.length; i++) arr[i] = 100 * arr[i] / n;
  delete res.alleles;
}

function fmtSize(n) { n = Math.abs(n); return n >= 1000 ? (n / 1000).toFixed(n >= 1e4 ? 0 : 1) + ' kb' : n + ' bp'; }
function describe(a, cuts = []) {
  if (a.cls === 'wt') return 'Unmodified';
  if (cuts.length > 1) {
    const parts = a.events.map(e => {
      if (e.kind === 'D') { const hit = cuts.map((c, i) => (e.pos <= c && e.pos + e.len >= c) ? i + 1 : 0).filter(Boolean); if (hit.length > 1) return `dropout cut ${hit[0]}→${hit[hit.length - 1]} (−${fmtSize(e.len)})`; }
      return null;
    });
    if (parts.some(Boolean)) return a.events.map((e, i) => parts[i] || (e.kind === 'D' ? `−${fmtSize(e.len)} del` : `+${e.len} ins${e.seq && e.len <= 12 ? ' ' + e.seq : ''}`)).join(' & ');
  }
  if (a.complex) return a.complex === 'inversion' ? 'Inversion / inverted split read' : 'Inter-chromosomal split read';
  return a.events.map(e => e.kind === 'D' ? `−${fmtSize(e.len)} del${e.split ? ' (split)' : ''}` : `+${e.len} ins${e.seq && e.len <= 12 ? ' ' + e.seq : ''}`).join(' & ');
}

/** Build aligned reference/allele strings for display around each cut (joined by ' ··· '). */
export function alleleAlignment(a, ref, cuts, half = 25) {
  cuts = Array.isArray(cuts) ? cuts : [cuts];
  // merge overlapping windows
  const wins = [];
  for (const c of [...cuts].sort((x, y) => x - y)) { const w = [c - half, c + half]; const last = wins[wins.length - 1]; if (last && w[0] <= last[1]) last[1] = w[1]; else wins.push(w); }
  const base = p => (ref && p >= ref.start && p < ref.start + ref.seq.length) ? ref.seq[p - ref.start].toUpperCase() : 'N';
  const ins = new Map(); for (const ev of (a.events || [])) if (ev.kind === 'I') ins.set(ev.pos, ev.seq || 'N'.repeat(ev.len));
  const del = p => (a.events || []).some(ev => ev.kind === 'D' && p >= ev.pos && p < ev.pos + ev.len);
  const cutSet = new Set(cuts);
  let top = '', bot = '', marks = '';
  wins.forEach(([s, e], wi) => {
    if (wi) { top += ' ··· '; bot += ' ··· '; marks += '     '; }
    for (let p = s; p < e; p++) {
      if (cutSet.has(p)) { top += '|'; bot += '|'; marks += 'c'; }
      if (ins.has(p)) { const x = ins.get(p); const show = x.length > 12 ? x.slice(0, 5) + '…' + x.slice(-5) : x; top += '-'.repeat(show.length); bot += show.toLowerCase(); marks += 'i'.repeat(show.length); }
      const b = base(p);
      top += b; if (del(p)) { bot += '-'; marks += 'd'; } else { bot += b; marks += ' '; }
    }
  });
  return { top, bot, marks };
}

export function editsToTSV(res) {
  const rows = [['sample', 'target', 'net_indel_bp', 'reads', 'percent_of_spanning_reads']];
  for (const d of res.netDist) rows.push([res.name, res.target?.label || '', d.size, d.count, d.pct.toFixed(4)]);
  return rows.map(r => r.join('\t')).join('\n') + '\n';
}
export function allelesToTSV(res) {
  const rows = [['sample', 'target', 'rank', 'allele', 'description', 'class', 'frame', 'net_indel_bp', 'reads', 'percent']];
  res.alleleList.forEach((a, i) => rows.push([res.name, res.target?.label || '', i + 1, a.key, a.label, a.cls, a.frame, a.net ?? '', a.count, a.pct.toFixed(4)]));
  return rows.map(r => r.join('\t')).join('\n') + '\n';
}
