// SVG charts for GenoView: multi-sample depth profile with call shading + hover, evidence bars.
export const SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'];
export const TYPE_COLORS = { DEL: '#e34948', LOSS: '#e34948', DUP: '#2a78d6', GAIN: '#2a78d6', INV: '#6250d6', TRA: '#eb6834', INS: '#1baf7a', CNV: '#eda100' };
const INK = '#0b0b0b', INK2 = '#52514e', INK3 = '#8a8984', GRID = '#ebeae5', POOLED = '#52514e';

const fmtBp = (v, span) => {
  if (span >= 2e6) return (v / 1e6).toFixed(span >= 2e7 ? 0 : 2).replace(/\.?0+$/, '') + ' Mb';
  if (span >= 5e3) return (v / 1e3).toFixed(span >= 5e4 ? 0 : 1).replace(/\.0$/, '') + ' kb';
  return Math.round(v).toLocaleString('en-US');
};
function niceTicks(a, b, n = 6) {
  const span = b - a; if (span <= 0) return [a];
  const step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => span / s <= n) || 10 * mag;
  const out = []; for (let v = Math.ceil(a / step) * step; v <= b + 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}
function quantile(arr, q) { const s = Float64Array.from(arr.filter(Number.isFinite)).sort(); if (!s.length) return 0; return s[Math.min(s.length - 1, Math.floor(q * (s.length - 1)))]; }
function median(arr) { return quantile(arr.filter(x => x > 0), 0.5); }
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

function seriesFor(res, mode, pooled, medians) {
  const out = [];
  const med = g => (medians && medians[g]) || median(Array.from(g === '*' ? res.pooled : res.groupDepth[g]));
  const ctrl = res.control;
  for (const g of res.groups) {
    const d = res.groupDepth[g]; let v;
    if (mode === 'raw') v = Array.from(d);
    else if (mode === 'norm') { const m = med(g) || 1; v = Array.from(d, x => 2 * x / m); }
    else { if (g === ctrl || !ctrl) continue; const m = med(g) || 1, mc = med(ctrl) || 1, c = res.groupDepth[ctrl]; v = Array.from(d, (x, i) => (c[i] + 0.5) > 0 ? Math.log2(((x + 0.5) / m) / ((c[i] + 0.5) / mc)) : 0); }
    out.push({ key: g, name: g + (g === ctrl && mode !== 'log2' ? ' (control)' : ''), color: res.colors[g], values: v, dash: null });
  }
  if (pooled && mode !== 'log2') {
    const d = res.pooled; const m = med('*') || 1;
    out.push({ key: '*', name: 'pooled (all BAMs)', color: POOLED, values: mode === 'raw' ? Array.from(d) : Array.from(d, x => 2 * x / m), dash: '5 3' });
  }
  return out;
}

export function depthChart(el, res, opts) {
  const { mode = 'raw', pooled = true, shade = [], selected = null, compact = false, onClick, medians } = opts;
  const W = Math.max(360, el.clientWidth || 900), H = compact ? 230 : 280;
  const m = { l: 50, r: 14, t: 24, b: 30 };
  const iw = W - m.l - m.r, ih = H - m.t - m.b;
  const series = seriesFor(res, mode, pooled, medians);
  const x = p => m.l + (p - res.start) / (res.end - res.start) * iw;
  let y0, y1;
  if (mode === 'log2') { y0 = -2.5; y1 = 2; }
  else { y0 = 0; const all = series.flatMap(s => s.values); y1 = Math.max(mode === 'norm' ? 4 : 5, quantile(all, 0.995) * 1.15); }
  const y = v => m.t + ih - (Math.max(y0, Math.min(y1, v)) - y0) / (y1 - y0) * ih;
  const bx = i => x(res.start + (i + 0.5) * res.binSize);

  let s = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" font-family="system-ui,-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif" font-size="11">`;
  s += `<rect x="0" y="0" width="${W}" height="${H}" fill="#fcfcfb"/>`;
  // shading (behind everything)
  for (const c of shade) {
    if (c.end < res.start || c.start > res.end) continue;
    const x0 = Math.max(m.l, x(c.start)), x1 = Math.min(m.l + iw, x(Math.max(c.end, c.start + 1)));
    const col = TYPE_COLORS[c.type] || TYPE_COLORS.CNV;
    const isSel = selected && selected.id === c.id;
    s += `<rect x="${x0}" y="${m.t}" width="${Math.max(1.5, x1 - x0)}" height="${ih}" fill="${col}" fill-opacity="${isSel ? 0.16 : 0.09}"${isSel ? ` stroke="${col}" stroke-width="1"` : ''}/>`;
    if (x1 - x0 > 26 || isSel) s += `<text x="${(x0 + x1) / 2}" y="${m.t - 6}" text-anchor="middle" fill="${INK2}" font-size="10.5" font-weight="600">${c.source === 'depth' ? 'CNV ' + c.type.toLowerCase() : c.type}</text>`;
    if (compact) for (const bp of [c.start, c.end]) if (bp >= res.start && bp <= res.end) s += `<line x1="${x(bp)}" x2="${x(bp)}" y1="${m.t}" y2="${m.t + ih}" stroke="${col}" stroke-width="1" stroke-dasharray="3 2"/>`;
  }
  // grid + y axis
  const yt = niceTicks(y0, y1, compact ? 4 : 5);
  for (const v of yt) s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${y(v)}" y2="${y(v)}" stroke="${GRID}"/><text x="${m.l - 6}" y="${y(v) + 3.5}" text-anchor="end" fill="${INK2}">${+v.toFixed(2)}</text>`;
  const refLine = mode === 'norm' ? 2 : mode === 'log2' ? 0 : null;
  if (refLine != null) s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${y(refLine)}" y2="${y(refLine)}" stroke="${INK3}" stroke-dasharray="4 3"/>`;
  const ylab = mode === 'raw' ? 'depth (×)' : mode === 'norm' ? 'copy number' : 'log2 ratio vs control';
  s += `<text transform="translate(12 ${m.t + ih / 2}) rotate(-90)" text-anchor="middle" fill="${INK2}">${ylab}</text>`;
  // x axis
  const xt = niceTicks(res.start, res.end, Math.max(3, Math.floor(iw / 110)));
  s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${m.t + ih}" y2="${m.t + ih}" stroke="${INK3}"/>`;
  for (const v of xt) s += `<line x1="${x(v)}" x2="${x(v)}" y1="${m.t + ih}" y2="${m.t + ih + 4}" stroke="${INK3}"/><text x="${x(v)}" y="${m.t + ih + 16}" text-anchor="middle" fill="${INK2}">${fmtBp(v, res.end - res.start)}</text>`;
  s += `<text x="${m.l + iw}" y="${H - 2}" text-anchor="end" fill="${INK3}" font-size="10">${esc(res.chr)}</text>`;
  // lines
  const step = Math.max(1, Math.floor(res.nbins / (iw * 1.5)));
  for (const se of series) {
    let d = '';
    for (let i = 0; i < res.nbins; i += step) {
      let v = se.values[i];
      if (step > 1) { let a = 0; for (let k = i; k < Math.min(res.nbins, i + step); k++) a += se.values[k]; v = a / Math.min(step, res.nbins - i); }
      d += (i ? 'L' : 'M') + bx(Math.min(res.nbins - 1, i + (step - 1) / 2)).toFixed(1) + ' ' + y(v).toFixed(1);
    }
    s += `<path d="${d}" fill="none" stroke="${se.color}" stroke-width="${se.key === '*' ? 1.5 : 2}" stroke-linejoin="round" stroke-linecap="round"${se.dash ? ` stroke-dasharray="${se.dash}"` : ''}/>`;
  }
  // hover layer
  s += `<g class="hover" style="display:none"><line class="hx" y1="${m.t}" y2="${m.t + ih}" stroke="${INK2}" stroke-width="1"/>${series.map(se => `<circle r="4" fill="${se.color}" stroke="#fcfcfb" stroke-width="2"/>`).join('')}</g>`;
  s += `<rect class="hit" x="${m.l}" y="${m.t}" width="${iw}" height="${ih}" fill="transparent" style="cursor:${onClick ? 'pointer' : 'crosshair'}"/>`;
  s += `</svg>`;

  const legend = `<div class="legend">${series.map(se => `<span><i class="ln" style="border-top:2px ${se.dash ? 'dashed' : 'solid'} ${se.color}"></i>${esc(se.name)}</span>`).join('')}${refLine != null ? `<span><i class="ln" style="border-top:1px dashed ${INK3}"></i>${mode === 'norm' ? 'diploid (2)' : 'no change'}</span>` : ''}</div>`;
  el.innerHTML = legend + s + '<div class="tip"></div>';
  const svg = el.querySelector('svg'), hover = svg.querySelector('.hover'), tip = el.querySelector('.tip'), hit = svg.querySelector('.hit');
  const circles = [...hover.querySelectorAll('circle')], hx = hover.querySelector('.hx');
  const toLocal = e => { const r = svg.getBoundingClientRect(); return (e.clientX - r.left) * (W / r.width); };
  hit.addEventListener('mousemove', e => {
    const px = toLocal(e); const pos = res.start + (px - m.l) / iw * (res.end - res.start);
    const i = Math.max(0, Math.min(res.nbins - 1, Math.floor((pos - res.start) / res.binSize)));
    const cx = bx(i); hover.style.display = ''; hx.setAttribute('x1', cx); hx.setAttribute('x2', cx);
    series.forEach((se, k) => { circles[k].setAttribute('cx', cx); circles[k].setAttribute('cy', y(se.values[i])); });
    const binS = res.start + i * res.binSize;
    const inCalls = shade.filter(c => binS < c.end && binS + res.binSize > c.start).map(c => c.source === 'depth' ? 'CNV ' + c.type.toLowerCase() : c.type);
    tip.innerHTML = `<b>${esc(res.chr)}:${(binS + 1).toLocaleString('en-US')}-${Math.min(res.end, binS + res.binSize).toLocaleString('en-US')}</b>${inCalls.length ? ' · ' + inCalls.join(', ') : ''}<br>` +
      series.map(se => `<span class="swatch" style="background:${se.color}"></span> ${esc(se.name)}: <b>${se.values[i].toFixed(mode === 'raw' ? 1 : 2)}</b>`).join('<br>');
    tip.style.display = 'block';
    const r = el.getBoundingClientRect(), sr = svg.getBoundingClientRect();
    const left = (cx / W) * sr.width + (sr.left - r.left);
    tip.style.left = (left + 14 + tip.offsetWidth > r.width ? left - tip.offsetWidth - 14 : left + 14) + 'px';
    tip.style.top = (sr.top - r.top + 30) + 'px';
  });
  hit.addEventListener('mouseleave', () => { hover.style.display = 'none'; tip.style.display = 'none'; });
  if (onClick) hit.addEventListener('click', e => { const px = toLocal(e); onClick(res.start + (px - m.l) / iw * (res.end - res.start)); });
}

const EV = [['split', 'split reads', '#184f95'], ['disc', 'discordant pairs', '#3987e5'], ['cigar', 'CIGAR indel', '#86b6ef']];
export function evidenceChart(el, call, res) {
  // per-group support is pooled across evidence kinds in the engine; split by kind proportionally for display only when single group
  const groups = res.groups;
  const W = Math.max(260, el.clientWidth || 320), rowH = 24, m = { l: 92, r: 40, t: 6, b: 18 };
  const H = m.t + m.b + rowH * groups.length;
  const tot = call.evidence ? (call.evidence.split || 0) + (call.evidence.disc || 0) + (call.evidence.cigar || 0) : 0;
  const max = Math.max(1, ...groups.map(g => call.perGroup[g]?.support || 0));
  const iw = W - m.l - m.r;
  let s = `<svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" font-family="system-ui,-apple-system,Segoe UI,Roboto,sans-serif" font-size="11">`;
  groups.forEach((g, k) => {
    const n = call.perGroup[g]?.support || 0; const yy = m.t + k * rowH;
    s += `<text x="${m.l - 8}" y="${yy + rowH / 2 + 4}" text-anchor="end" fill="${INK2}">${esc(g.length > 14 ? g.slice(0, 13) + '…' : g)}</text>`;
    let x0 = m.l;
    const w = n / max * iw;
    // apportion this group's reads over evidence kinds by the call-wide mix
    EV.forEach(([key, , col], j) => {
      const share = tot ? (call.evidence[key] || 0) / tot : 0; const ww = w * share;
      if (ww <= 0) return;
      s += `<rect x="${x0}" y="${yy + 5}" width="${Math.max(0, ww - 2)}" height="${rowH - 10}" rx="2" fill="${col}"><title>${esc(g)} · ${EV[j][1]} ≈ ${Math.round(n * share)}</title></rect>`;
      x0 += ww;
    });
    s += `<text x="${m.l + w + 5}" y="${yy + rowH / 2 + 4}" fill="${INK}" font-weight="600">${n}</text>`;
    s += `<rect x="${m.l - 86}" y="${yy + rowH / 2 - 4}" width="8" height="8" rx="2" fill="${res.colors[g]}"/>`;
  });
  s += `<line x1="${m.l}" x2="${m.l}" y1="${m.t}" y2="${H - m.b}" stroke="${INK3}"/></svg>`;
  el.innerHTML = `<div class="legend">${EV.map(([k, l, c]) => `<span><span class="swatch" style="background:${c}"></span>${l} <b>${call.evidence?.[k] || 0}</b></span>`).join('')}<span>soft-clips near breakpoints <b>${call.clipSupport ?? 0}</b></span></div>` + s;
}
