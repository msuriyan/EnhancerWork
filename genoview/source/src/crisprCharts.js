// Single-sample charts for CRISPR outcomes: standardised indel-size histogram / ECDF, deletion footprint.
const INK = '#0b0b0b', INK2 = '#52514e', INK3 = '#8a8984', GRID = '#ebeae5', SURF = '#fcfcfb';
export const FRAME_COLORS = { wt: '#b5b4ae', 'in-frame': '#2a78d6', frameshift: '#eb6834' };
const FONT = 'font-family="system-ui,-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif" font-size="11"';

function niceMax(v) { if (v <= 0) return 1; const m = Math.pow(10, Math.floor(Math.log10(v))); for (const k of [1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]) if (k * m >= v) return k * m; return 10 * m; }
function ticks(max, n = 4) { const st = niceMax(max / n); const out = []; for (let v = 0; v <= max + 1e-9; v += st) out.push(+v.toFixed(6)); return out; }

function attachTip(el, svg, W, items) {
  const tip = el.querySelector('.tip');
  svg.querySelectorAll('[data-tip]').forEach(n => {
    n.addEventListener('mousemove', e => {
      tip.innerHTML = items[+n.dataset.tip]; tip.style.display = 'block';
      const r = el.getBoundingClientRect();
      let x = e.clientX - r.left + 12; if (x + tip.offsetWidth > r.width) x = e.clientX - r.left - tip.offsetWidth - 12;
      tip.style.left = x + 'px'; tip.style.top = (e.clientY - r.top - 10) + 'px';
    });
    n.addEventListener('mouseleave', () => { tip.style.display = 'none'; });
  });
}

/**
 * opts: range (bp), bin (bp), basis 'all'|'edited', showWT, cumulative
 */
export function indelHistogram(el, res, opts) {
  const { range = 30, bin = 1, basis = 'all', showWT = true, cumulative = false } = opts;
  const W = Math.max(380, el.clientWidth || 640), H = 290, m = { l: 48, r: 16, t: 26, b: 42 };
  const iw = W - m.l - m.r, ih = H - m.t - m.b;
  const n = res.stats.spanning || 1;
  const wtCount = res.alleleList.find(a => a.cls === 'wt')?.count || 0;
  const denom = basis === 'edited' ? Math.max(1, n - wtCount) : n;
  // per-read net sizes with frame classification
  const bins = new Map(); // key -> {lo, hi, wt, inf, fs, label}
  const keyOf = s => s < -range ? 'lo' : s > range ? 'hi' : Math.floor((s + (bin > 1 ? Math.floor(bin / 2) : 0)) / bin);
  const binLo = k => k === 'lo' ? -Infinity : k === 'hi' ? range + 1 : k * bin - (bin > 1 ? Math.floor(bin / 2) : 0);
  for (const a of res.alleleList) {
    if (a.net == null) continue;
    if (a.cls === 'wt' && (!showWT || basis === 'edited')) continue;
    const k = keyOf(a.net);
    let b = bins.get(k); if (!b) { b = { k, wt: 0, inf: 0, fs: 0, alleles: [] }; bins.set(k, b); }
    if (a.cls === 'wt') b.wt += a.count; else if (a.net % 3 === 0) b.inf += a.count; else b.fs += a.count;
    b.alleles.push(a);
  }
  const nb = Math.ceil((2 * range + 1) / bin);
  const kMin = keyOf(-range), kMax = keyOf(range);
  const slots = ['lo']; for (let k = kMin; k <= kMax; k++) slots.push(k); slots.push('hi');
  const slotW = iw / slots.length;
  const xOfSlot = i => m.l + i * slotW;
  const unres = res.alleleList.filter(a => a.net == null).reduce((s, a) => s + a.count, 0);
  let s = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" ${FONT}><rect width="${W}" height="${H}" fill="${SURF}"/>`;
  const tips = [];

  if (cumulative) {
    // empirical CDF of net indel size over spanning reads (or edited reads)
    const pts = []; let acc = 0;
    const sizes = res.netDist
      .map(d => ({ size: d.size, c: basis === 'edited' && d.size === 0 ? d.count - d.wt : d.count })).filter(d => d.c > 0);
    const clampS = v => Math.max(-range - 1, Math.min(range + 1, v));
    const x = v => m.l + (clampS(v) + range + 1) / (2 * range + 2) * iw;
    const y = v => m.t + ih - v / 100 * ih;
    pts.push([x(-range - 1), y(0)]);
    for (const d of sizes) { const yb = acc; acc += 100 * d.c / denom; pts.push([x(d.size), y(yb)], [x(d.size), y(acc)]); tips.push(`net <b>${d.size > 0 ? '+' : ''}${d.size}</b> bp → cumulative <b>${acc.toFixed(1)}%</b>`); }
    pts.push([x(range + 1), y(acc)]);
    for (const v of [0, 25, 50, 75, 100]) s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${y(v)}" y2="${y(v)}" stroke="${GRID}"/><text x="${m.l - 6}" y="${y(v) + 3.5}" text-anchor="end" fill="${INK2}">${v}</text>`;
    s += `<line x1="${x(0)}" x2="${x(0)}" y1="${m.t}" y2="${m.t + ih}" stroke="${INK3}" stroke-dasharray="3 3"/>`;
    s += `<path d="M${pts.map(p => p[0].toFixed(1) + ' ' + p[1].toFixed(1)).join('L')}" fill="none" stroke="#2a78d6" stroke-width="2" stroke-linejoin="round"/>`;
    const xt = []; const st = range >= 50 ? 25 : range >= 20 ? 10 : 5; for (let v = -Math.floor(range / st) * st; v <= range; v += st) if (Math.abs(v) < range - st / 2) xt.push(v);
    for (const v of xt) s += `<text x="${x(v)}" y="${m.t + ih + 16}" text-anchor="middle" fill="${INK2}">${v > 0 ? '+' + v : v}</text>`;
    s += `<text x="${x(-range - 1)}" y="${m.t + ih + 16}" text-anchor="start" fill="${INK3}" font-size="10">≤−${range + 1}</text><text x="${x(range + 1)}" y="${m.t + ih + 16}" text-anchor="end" fill="${INK3}" font-size="10">≥+${range + 1}</text>`;
  } else {
    let ymax = 0; for (const b of bins.values()) ymax = Math.max(ymax, 100 * (b.wt + b.inf + b.fs) / denom);
    ymax = niceMax(ymax * 1.08);
    const y = v => m.t + ih - v / ymax * ih;
    for (const v of ticks(ymax)) s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${y(v)}" y2="${y(v)}" stroke="${GRID}"/><text x="${m.l - 6}" y="${y(v) + 3.5}" text-anchor="end" fill="${INK2}">${v}</text>`;
    const i0 = slots.indexOf(keyOf(0));
    s += `<line x1="${xOfSlot(i0) + slotW / 2}" x2="${xOfSlot(i0) + slotW / 2}" y1="${m.t}" y2="${m.t + ih}" stroke="${INK3}" stroke-dasharray="3 3"/>`;
    const gap = slotW > 6 ? 1 : 0, bw = Math.max(1, slotW - 2 * gap);
    const labelled = [...bins.values()].sort((a, b) => (b.wt + b.inf + b.fs) - (a.wt + a.inf + a.fs)).slice(0, 4).map(b => b.k);
    slots.forEach((k, i) => {
      const b = bins.get(k); if (!b) return;
      const x0 = xOfSlot(i) + gap;
      let acc = 0;
      for (const [part, col] of [['wt', FRAME_COLORS.wt], ['inf', FRAME_COLORS['in-frame']], ['fs', FRAME_COLORS.frameshift]]) {
        if (!b[part]) continue;
        const v = 100 * b[part] / denom, ya = y(acc + v), yb = y(acc);
        const top = acc + v === 100 * (b.wt + b.inf + b.fs) / denom;
        s += top && yb - ya > 4 ? `<path d="M${x0} ${yb}V${ya + 2}q0 -2 2 -2h${bw - 4}q2 0 2 2V${yb}z" fill="${col}"/>` : `<rect x="${x0}" y="${ya}" width="${bw}" height="${Math.max(0.5, yb - ya - (acc ? 1 : 0))}" fill="${col}"/>`;
        acc += v;
      }
      const lbl = k === 'lo' ? `≤ −${range + 1} bp` : k === 'hi' ? `≥ +${range + 1} bp` : bin === 1 ? `${k > 0 ? '+' : ''}${k} bp` : `${binLo(k)} … ${binLo(k) + bin - 1} bp`;
      const tot = b.wt + b.inf + b.fs;
      tips.push(`<b>${lbl}</b> · ${(100 * tot / denom).toFixed(2)}% (${tot} reads)<br>` + b.alleles.slice(0, 6).map(a => `${a.pct.toFixed(1)}% · ${a.label}`).join('<br>') + (b.alleles.length > 6 ? `<br>+${b.alleles.length - 6} more` : ''));
      s += `<rect data-tip="${tips.length - 1}" x="${xOfSlot(i)}" y="${m.t}" width="${slotW}" height="${ih}" fill="transparent"/>`;
      if (labelled.includes(k) && 100 * tot / denom >= 0.5) s += `<text x="${x0 + bw / 2}" y="${y(100 * tot / denom) - 4}" text-anchor="middle" fill="${INK}" font-weight="600" font-size="10.5">${(100 * tot / denom).toFixed(1)}%</text>`;
    });
    // x labels
    const step = Math.max(1, Math.ceil(28 / slotW));
    const lblEvery = [1, 2, 5, 10, 20, 25, 50].find(v => v >= step * bin) || 50;
    slots.forEach((k, i) => {
      const cx = xOfSlot(i) + slotW / 2;
      if (k === 'lo' || k === 'hi') { s += `<text x="${cx}" y="${m.t + ih + 15}" text-anchor="middle" fill="${INK3}" font-size="10">${k === 'lo' ? '≤−' + (range + 1) : '≥+' + (range + 1)}</text>`; return; }
      const v = binLo(k) + (bin > 1 ? Math.floor(bin / 2) : 0);
      if (v % lblEvery === 0 && Math.abs(v) <= range - Math.max(bin, lblEvery / 2)) s += `<line x1="${cx}" x2="${cx}" y1="${m.t + ih}" y2="${m.t + ih + 4}" stroke="${INK3}"/><text x="${cx}" y="${m.t + ih + 15}" text-anchor="middle" fill="${INK2}">${v > 0 ? '+' + v : v}</text>`;
    });
  }
  s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${m.t + ih}" y2="${m.t + ih}" stroke="${INK3}"/>`;
  s += `<text x="${m.l + iw / 2}" y="${H - 6}" text-anchor="middle" fill="${INK2}">net indel size per read (bp; − deletion, + insertion)</text>`;
  if (res.target && res.target.label) s += `<text x="${m.l + iw / 2}" y="14" text-anchor="middle" fill="${INK2}">${esc(res.target.label)}</text>`;
  s += `<text transform="translate(12 ${m.t + ih / 2}) rotate(-90)" text-anchor="middle" fill="${INK2}">${cumulative ? 'cumulative ' : ''}% of ${basis === 'edited' ? 'edited' : 'spanning'} reads</text>`;
  s += `<text x="${m.l}" y="14" fill="${INK}" font-weight="600" font-size="12">${esc(res.name)}</text><text x="${m.l + iw}" y="14" text-anchor="end" fill="${INK2}">n = ${res.stats.spanning.toLocaleString('en-US')} reads${basis === 'edited' ? ` · ${(n - wtCount).toLocaleString('en-US')} edited` : ''}${unres ? ` · ${unres} unresolved not plotted` : ''}</text>`;
  s += `</svg>`;
  const legend = cumulative ? '' : `<div class="legend">${showWT && basis !== 'edited' ? `<span><span class="swatch" style="background:${FRAME_COLORS.wt}"></span>unmodified</span>` : ''}<span><span class="swatch" style="background:${FRAME_COLORS['in-frame']}"></span>in-frame (net ≡ 0 mod 3)</span><span><span class="swatch" style="background:${FRAME_COLORS.frameshift}"></span>frameshift</span></div>`;
  el.innerHTML = legend + s + '<div class="tip"></div>';
  attachTip(el, el.querySelector('svg'), W, tips);
}

const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

/** % of spanning reads with each base deleted, around the cut(s); insertion sites as ticks. */
export function footprintChart(el, res, sites, half) {
  const W = Math.max(360, el.clientWidth || 520), H = 230, m = { l: 44, r: 14, t: 24, b: 36 };
  const iw = W - m.l - m.r, ih = H - m.t - m.b;
  const o = res.origin, cuts = res.cuts.length ? res.cuts : [o];
  const lo = Math.max(res.regStart, Math.min(...cuts) - half), hi = Math.min(res.regStart + res.delCov.length, Math.max(...cuts) + half);
  const x = p => m.l + (p - lo) / Math.max(1, hi - lo) * iw;
  const at = (arr, p) => arr[p - res.regStart] || 0;
  let ymax = 0; for (let p = lo; p < hi; p++) ymax = Math.max(ymax, at(res.delCov, p), at(res.insAt, p));
  ymax = niceMax(Math.max(1, ymax * 1.1));
  const y = v => m.t + ih - v / ymax * ih;
  let s = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" ${FONT}><rect width="${W}" height="${H}" fill="${SURF}"/>`;
  for (const a of res.areas) if (a[1] > lo && a[0] < hi) s += `<rect x="${x(Math.max(lo, a[0]))}" y="${m.t}" width="${Math.max(1, x(Math.min(hi, a[1])) - x(Math.max(lo, a[0])))}" height="${ih}" fill="#eda100" fill-opacity="0.10"/>`;
  for (const st of sites || []) {
    if (st.start == null || st.end < lo || st.start > hi) continue;
    s += `<rect x="${x(Math.max(lo, st.start))}" y="${m.t}" width="${x(Math.min(hi, st.end)) - x(Math.max(lo, st.start))}" height="6" fill="#2a78d6" fill-opacity="0.55"><title>guide ${st.strand} ${st.guide || ''}</title></rect>`;
  }
  for (const v of ticks(ymax)) s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${y(v)}" y2="${y(v)}" stroke="${GRID}"/><text x="${m.l - 6}" y="${y(v) + 3.5}" text-anchor="end" fill="${INK2}">${v}</text>`;
  let d = `M${x(lo)} ${y(0)}`;
  const step = Math.max(1, Math.floor((hi - lo) / iw));
  for (let p = lo; p < hi; p += step) { const v = at(res.delCov, p); d += `L${x(p).toFixed(1)} ${y(v).toFixed(1)}L${x(Math.min(hi, p + step)).toFixed(1)} ${y(v).toFixed(1)}`; }
  d += `L${x(hi)} ${y(0)}Z`;
  s += `<path d="${d}" fill="#e34948" fill-opacity="0.22" stroke="#e34948" stroke-width="1.5" stroke-linejoin="round"/>`;
  for (let p = lo; p < hi; p++) { const v = at(res.insAt, p); if (v > 0.05) s += `<rect x="${x(p) - 1.5}" y="${y(v)}" width="3" height="${y(0) - y(v)}" fill="#1baf7a" rx="1"><title>${v.toFixed(1)}% reads with insertion at ${p - o >= 0 ? '+' : ''}${p - o}</title></rect>`; }
  cuts.forEach((c, i) => { if (c < lo || c > hi) return; s += `<line x1="${x(c)}" x2="${x(c)}" y1="${m.t}" y2="${m.t + ih}" stroke="${INK}" stroke-width="1" stroke-dasharray="4 3"/><text x="${x(c)}" y="${m.t - 6}" text-anchor="middle" fill="${INK}" font-weight="600" font-size="10.5">${cuts.length > 1 ? 'cut ' + (i + 1) : 'cut'}</text>`; });
  s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${m.t + ih}" y2="${m.t + ih}" stroke="${INK3}"/>`;
  const span = hi - lo;
  const st = [5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000].find(v => span / v <= 8) || 20000;
  for (let v = Math.ceil((lo - o) / st) * st; v <= hi - o; v += st) s += `<line x1="${x(o + v)}" x2="${x(o + v)}" y1="${m.t + ih}" y2="${m.t + ih + 4}" stroke="${INK3}"/><text x="${x(o + v)}" y="${m.t + ih + 15}" text-anchor="middle" fill="${INK2}">${v > 0 ? '+' + v : v}</text>`;
  s += `<text x="${m.l + iw / 2}" y="${H - 4}" text-anchor="middle" fill="${INK2}">position relative to ${cuts.length > 1 ? 'cut 1' : 'cut site'} (bp)</text>`;
  s += `<text transform="translate(12 ${m.t + ih / 2}) rotate(-90)" text-anchor="middle" fill="${INK2}">% of reads</text></svg>`;
  el.innerHTML = `<div class="legend"><span><span class="swatch" style="background:#e34948;opacity:.6"></span>base deleted</span><span><span class="swatch" style="background:#1baf7a"></span>insertion site</span><span><span class="swatch" style="background:#eda100;opacity:.35"></span>counting window</span><span><span class="swatch" style="background:#2a78d6;opacity:.6"></span>guide</span></div>` + s;
}

/** Per-sample bar plot (one bar per sample), returns an SVG string. rows: [{name, value, n, group}] */
export function sampleBarSVG(rows, o) {
  const W = o.width || Math.max(360, 70 * rows.length + 140), H = o.height || 340;
  const font = o.font || 'Arial, Helvetica, sans-serif', fs = o.fontSize || 12;
  const pal = ['#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3', '#937860', '#da8bc3', '#8c8c8c'];
  const groups = [...new Set(rows.map(r => r.group))];
  const colorOf = r => o.byGroup && groups.length > 1 ? pal[groups.indexOf(r.group) % pal.length] : o.color || '#4c72b0';
  const maxLab = Math.max(...rows.map(r => r.name.length));
  const rot = maxLab * fs * 0.6 > (W - 120) / rows.length;
  const m = { l: 62, r: 16, t: o.title ? 34 : 14, b: (rot ? Math.min(150, maxLab * fs * 0.5 + 18) : fs + 18) + (o.byGroup && groups.length > 1 ? fs + 14 : 0) };
  const iw = W - m.l - m.r, ih = H - m.t - m.b;
  const vmax = Math.max(1, ...rows.map(r => r.value));
  const top = o.ymax || Math.min(100, Math.ceil(vmax * 1.15 / 10) * 10 || 10);
  const Y = v => m.t + ih - v / top * ih;
  const slot = iw / rows.length, bw = slot * 0.72;
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" font-family="${esc(font)}" font-size="${fs}"><rect width="${W}" height="${H}" fill="#ffffff"/>`;
  if (o.title) s += `<text x="${m.l}" y="20" font-size="${fs + 3}" font-weight="700" fill="#1a1a1a">${esc(o.title)}</text>`;
  const step = top <= 10 ? 2 : top <= 50 ? 10 : 20;
  for (let v = 0; v <= top + 1e-9; v += step) s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${Y(v)}" y2="${Y(v)}" stroke="#e8e8e8"/><text x="${m.l - 6}" y="${Y(v) + fs * 0.35}" text-anchor="end" fill="#333">${v}</text>`;
  rows.forEach((r, i) => {
    const cx = m.l + slot * (i + 0.5), y = Y(Math.max(0, r.value));
    s += `<rect x="${(cx - bw / 2).toFixed(2)}" y="${y.toFixed(2)}" width="${bw.toFixed(2)}" height="${Math.max(0, m.t + ih - y).toFixed(2)}" fill="${colorOf(r)}" stroke="#333" stroke-width="0.7"><title>${esc(r.name)}: ${r.value.toFixed(2)}% of ${r.n.toLocaleString('en-US')} reads</title></rect>`;
    if (o.showValues !== false) s += `<text x="${cx}" y="${y - 4}" text-anchor="middle" font-size="${fs * 0.85}" fill="#1a1a1a">${r.value.toFixed(1)}</text>`;
    const ly = m.t + ih + fs + 4;
    s += rot ? `<text transform="translate(${cx} ${m.t + ih + 8}) rotate(-45)" text-anchor="end" fill="#1a1a1a">${esc(r.name)}</text>` : `<text x="${cx}" y="${ly}" text-anchor="middle" fill="#1a1a1a">${esc(r.name)}</text>`;
  });
  s += `<line x1="${m.l}" x2="${m.l}" y1="${m.t}" y2="${m.t + ih}" stroke="#333"/><line x1="${m.l}" x2="${m.l + iw}" y1="${m.t + ih}" y2="${m.t + ih}" stroke="#333"/>`;
  s += `<text transform="translate(16 ${m.t + ih / 2}) rotate(-90)" text-anchor="middle" fill="#1a1a1a">${esc(o.ylabel || '% of reads')}</text>`;
  if (o.byGroup && groups.length > 1) { let lx = m.l; const ly = H - 8; groups.forEach((g, k) => { s += `<rect x="${lx}" y="${ly - fs + 2}" width="${fs}" height="${fs}" fill="${pal[k % pal.length]}" stroke="#333" stroke-width="0.7"/><text x="${lx + fs + 5}" y="${ly}" fill="#1a1a1a">${esc(g)}</text>`; lx += fs + 14 + String(g).length * fs * 0.58; }); }
  return s + '</svg>';
}
