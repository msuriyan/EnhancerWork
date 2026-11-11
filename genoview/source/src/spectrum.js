// Publication figure: indel / SV size spectrum (small → large), any number of series, fully styleable.
// spec = { name, unitLabel:'reads'|'calls', total, edited, events, series:[{key,label,color,map:Map(len->count)}],
//          perSample?: [spec, ...] (for mean ± SD), pooled?: 'reads'|'mean'|null, subtitle?: string }

export const BIN_PRESETS = {
  equal: { label: 'Equal-width bins (set bin width and maximum below)', edges: null, equal: true },
  ranges: { label: 'Size ranges: 1–20, 21–50, 51–100, 101–200, 201–500, 501–1000, 1001–1500, 1501–2000, >2000', edges: [20, 50, 100, 200, 500, 1000, 1500, 2000] },
  log: { label: '1–5 bp single, then 6–10, 11–20, 21–50 … >10 kb', edges: [1, 2, 3, 4, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000] },
  fine: { label: '1–20 bp single, then log-style', edges: [...Array.from({ length: 20 }, (_, i) => i + 1), 50, 100, 200, 500, 1000, 2000, 5000, 10000] },
  sv: { label: 'SV classes (50–100, 100–500, 500–1k, 1–5k, 5–10k, 10–100k, 100k–1M, >1M)', edges: [49, 100, 500, 1000, 5000, 10000, 100000, 1000000] },
  decade: { label: 'Decades (1–9, 10–99, 100–999, 1–9.9 kb, …)', edges: [9, 99, 999, 9999, 99999, 999999] },
  exact: { label: 'Every length (1-bp bins up to the largest)', edges: null },
};

export const DEFAULT_STYLE = {
  title: '', subtitle: 'auto', xLabel: 'Indel length (bp)', yLabel: 'auto',
  font: 'Arial, Helvetica, sans-serif', titleSize: 16, subtitleSize: 11, axisTitleSize: 13, tickSize: 11, valueSize: 10, legendSize: 11,
  textColor: '#1a1a1a', axisColor: '#555555', gridColor: '#e6e6e6', background: '#ffffff', transparent: false,
  layout: 'histogram',          // histogram | grouped | stacked | mirror | cumulative
  measure: 'pct',               // pct (of reads / calls) | count | pctEvents
  bins: 'equal', customEdges: '20, 50, 100, 200, 500, 1000, 1500, 2000', trimEmpty: false,
  showValues: false, valueDecimals: 1, showGrid: true, legend: 'top-right', errorBars: true,
  width: 760, height: 430, barWidth: 0.78, rotateTicks: 'auto', lineWidth: 1, dpi: 300,
  series: {},                   // key -> { label, color, show }
  binWidth: 100, binMax: 0, overflowBin: true, mergeSeries: false,
  barStroke: '#333333', barStrokeWidth: 1, fillOpacity: 0.85, xTickStep: 0,
};

const fmtLen = n => n >= 1e6 ? +(n / 1e6).toFixed(1) + ' Mb' : n >= 10000 ? +(n / 1000).toFixed(1) + ' kb' : n + '';
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
function niceMax(v) { if (v <= 0) return 1; const m = Math.pow(10, Math.floor(Math.log10(v))); for (const k of [1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]) if (k * m >= v * 1.0001) return k * m; return 10 * m; }
function ticks(max, n = 5) { const raw = max / n; const m = Math.pow(10, Math.floor(Math.log10(raw))); const st = [1, 2, 2.5, 5, 10].map(k => k * m).find(k => k >= raw) || 10 * m; const out = []; for (let v = 0; v <= max + st * 1e-6; v += st) out.push(+v.toFixed(8)); return out; }

export function seriesStyle(S, se) {
  const o = S.series?.[se.key] || {};
  return { label: o.label ?? se.label, color: o.color ?? se.color, show: o.show ?? true };
}

function withMerge(spec, S) {
  if (!S.mergeSeries || spec.series.length < 2) return spec;
  const map = new Map(); for (const se of spec.series) for (const [l, c] of se.map) map.set(l, (map.get(l) || 0) + c);
  const merged = { ...spec, series: [{ key: 'all', label: 'All indels', color: '#c9c9c9', map }] };
  if (spec.perSample) merged.perSample = spec.perSample.map(sp => withMerge(sp, S));
  return merged;
}
/** in deletions-only mode, default titles say "Deletion" instead of "Indel" (custom titles are kept) */
const delWord = (spec, t) => spec.delOnly && /^indel size (spectrum|distribution)$/i.test(t) ? t.replace(/^indel/i, 'Deletion') : t;
function maxLength(spec) { return Math.max(1, ...spec.series.flatMap(se => [...se.map.keys()])); }
function binEdges(spec, S) {
  let edges;
  if (S.bins === 'equal') {
    const w = Math.max(1, +S.binWidth || 100);
    const max = +S.binMax > 0 ? +S.binMax : Math.ceil(maxLength(spec) / w) * w;
    edges = []; for (let e = w; e <= max + 1e-9 && edges.length < 2000; e += w) edges.push(e);
    const bins = []; let lo = 1;
    for (const e of edges) { bins.push({ lo, hi: e, x0: lo - 1, x1: e }); lo = e + 1; }
    if (S.overflowBin) bins.push({ lo, hi: Infinity, x0: lo - 1, x1: lo - 1 + w, overflow: true });
    for (const b of bins) b.label = b.hi === Infinity ? `>${fmtLen(b.lo - 1)}` : `${fmtLen(b.lo)}–${fmtLen(b.hi)}`;
    return bins;
  }
  if (S.bins === 'custom') edges = String(S.customEdges).split(/[\s,;]+/).map(Number).filter(x => x > 0).sort((a, b) => a - b);
  else edges = BIN_PRESETS[S.bins]?.edges;
  if (!edges || !edges.length) {
    const maxLen = Math.max(1, ...spec.series.flatMap(se => [...se.map.keys()]));
    edges = Array.from({ length: Math.min(maxLen, 2000) }, (_, i) => i + 1);
  }
  const bins = []; let lo = 1;
  for (const e of edges) { if (e < lo) continue; bins.push({ lo, hi: e }); lo = e + 1; }
  bins.push({ lo, hi: Infinity });
  for (const b of bins) b.label = b.hi === Infinity ? `>${fmtLen(b.lo - 1)}` : b.lo === b.hi ? fmtLen(b.lo) : `${fmtLen(b.lo)}–${fmtLen(b.hi)}`;
  return bins;
}
function binCounts(spec, bins) {
  return spec.series.map(se => { const a = new Float64Array(bins.length); for (const [len, c] of se.map) { const i = bins.findIndex(b => len >= b.lo && len <= b.hi); if (i >= 0) a[i] += c; } return a; });
}

/** values per series per bin for the chosen measure, plus SD across samples (mean pooling) */
export function computeTable(spec, S) {
  spec = withMerge(spec, S);
  let bins = binEdges(spec, S);
  const denom = sp => S.measure === 'count' ? 1 : S.measure === 'pctEvents' ? Math.max(1, sp.events) / 100 : Math.max(1, sp.total) / 100;
  let vals, sds = null;
  if (spec.pooled === 'mean' && spec.perSample?.length && S.measure !== 'count') {
    // % of events: a sample with no events (e.g. uninjected) has no distribution, so it is left out of the mean
    const use = S.measure === 'pctEvents' ? spec.perSample.filter(sp => sp.events > 0) : spec.perSample;
    const per = (use.length ? use : spec.perSample).map(sp => binCounts(sp, bins).map(a => a.map(c => c / denom(sp))));
    vals = spec.series.map((_, k) => bins.map((_, i) => per.reduce((s, p) => s + p[k][i], 0) / per.length));
    sds = spec.series.map((_, k) => bins.map((_, i) => { const m = vals[k][i]; const v = per.reduce((s, p) => s + (p[k][i] - m) ** 2, 0) / Math.max(1, per.length - 1); return Math.sqrt(v); }));
  } else {
    const raw = binCounts(spec, bins); const d = denom(spec);
    vals = raw.map(a => Array.from(a, c => c / d));
  }
  if (S.trimEmpty) {
    let last = 0; bins.forEach((b, i) => { if (vals.some(v => v[i] > 0)) last = i; });
    bins = bins.slice(0, last + 1); vals = vals.map(v => v.slice(0, last + 1)); if (sds) sds = sds.map(v => v.slice(0, last + 1));
  }
  return { bins, vals, sds, spec };
}

/** returns { svg } sized S.width × S.height */
export function renderSpectrum(spec0, S) {
  if (S.layout === 'histogram') return renderHistogram(spec0, S);
  const { bins, vals, sds, spec } = computeTable(spec0, S);
  const ser = spec.series.map((se, k) => ({ ...se, ...seriesStyle(S, se), v: vals[k], sd: sds?.[k] })).filter(x => x.show);
  const layout = ser.length < 2 && (S.layout === 'mirror' || S.layout === 'stacked') ? 'grouped' : S.layout;
  const W = S.width, H = S.height, F = `font-family="${esc(S.font)}"`;
  const unit = spec.unitLabel || 'reads';
  const yLabel = S.yLabel === 'auto' ? (layout === 'cumulative' ? `Cumulative % of ${S.measure === 'pctEvents' ? (spec.delOnly ? 'deletions' : 'events') : unit}` : S.measure === 'count' ? (unit[0].toUpperCase() + unit.slice(1)) : S.measure === 'pctEvents' ? (spec.delOnly ? '% of deletions' : '% of all indel events') : `% of ${unit}`) + (spec.pooled === 'mean' && S.measure !== 'count' && layout !== 'cumulative' ? ' (mean ± SD)' : '') : S.yLabel;
  const subtitle = S.subtitle === 'auto' ? (spec.subtitle || spec.name) : S.subtitle;
  const title = delWord(spec, S.title || 'Indel size spectrum');
  const nb = bins.length;
  const legendH = S.legend === 'bottom' ? S.legendSize + 14 : 0;
  const maxLabelChars = Math.max(...bins.map(b => b.label.length));
  const rotate = S.rotateTicks === 'auto' ? (maxLabelChars * S.tickSize * 0.6 > (W - 90) / nb * 0.95) : S.rotateTicks === 'yes';
  const tickBlock = rotate ? Math.min(100, maxLabelChars * S.tickSize * 0.55 + 10) : S.tickSize + 8;
  const top = 12 + (title ? S.titleSize + 6 : 0) + (subtitle ? S.subtitleSize + 6 : 0);
  const m = { l: 18 + S.axisTitleSize + 8 + S.tickSize * 3.4, r: 16, t: top + 10, b: tickBlock + S.axisTitleSize + 16 + legendH };
  const iw = Math.max(50, W - m.l - m.r), ih = Math.max(50, H - m.t - m.b);
  const slot = iw / nb, bw = slot * Math.min(0.98, Math.max(0.1, S.barWidth));
  const fmtV = v => S.measure === 'count' && layout !== 'cumulative' && spec.pooled !== 'mean' ? Math.round(v).toLocaleString('en-US') : v.toFixed(S.valueDecimals);
  // cumulative transform
  let C = ser;
  if (layout === 'cumulative') C = ser.map(x => { const tot = x.v.reduce((a, b) => a + b, 0) || 1; let acc = 0; return { ...x, v: x.v.map(v => (acc += v) / tot * 100), sd: null }; });
  // y-range
  let ymaxPos = 0, ymaxNeg = 0;
  const hi = (x, i) => x.v[i] + (S.errorBars && x.sd ? x.sd[i] : 0);
  if (layout === 'cumulative') ymaxPos = 100;
  else if (layout === 'stacked') ymaxPos = Math.max(...bins.map((_, i) => C.reduce((s, x) => s + x.v[i], 0)));
  else if (layout === 'mirror') { ymaxNeg = Math.max(...bins.map((_, i) => hi(C[0], i))); ymaxPos = Math.max(...bins.map((_, i) => C.slice(1).reduce((s, x) => s + x.v[i], 0) + (C.length === 2 && S.errorBars && C[1].sd ? C[1].sd[i] : 0))); }
  else ymaxPos = Math.max(...C.flatMap(x => bins.map((_, i) => hi(x, i))));
  const room = S.showValues ? 1.14 : 1.04;
  if (layout !== 'cumulative') { ymaxPos = ymaxPos > 0 ? niceMax(ymaxPos * room) : (layout === 'mirror' ? 0 : 1); ymaxNeg = ymaxNeg > 0 ? niceMax(ymaxNeg * room) : 0; }
  if (!ymaxPos && !ymaxNeg) ymaxPos = 1;
  const span = ymaxPos + ymaxNeg;
  const y = v => m.t + (ymaxPos - v) / span * ih;
  const y0 = y(0);
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" ${F}>`;
  if (!S.transparent) s += `<rect width="${W}" height="${H}" fill="${S.background}"/>`;
  let ty = 12;
  if (title) { ty += S.titleSize; s += `<text x="${m.l}" y="${ty}" font-size="${S.titleSize}" font-weight="700" fill="${S.textColor}">${esc(title)}</text>`; ty += 6; }
  if (subtitle) { ty += S.subtitleSize; s += `<text x="${m.l}" y="${ty}" font-size="${S.subtitleSize}" fill="${S.axisColor}">${esc(subtitle)}</text>`; }
  const step = (() => { const t = ticks(Math.max(ymaxPos, ymaxNeg) || 1e-4, ymaxNeg && ymaxPos ? 4 : 5); return t[1] - t[0]; })();
  const tks = []; for (let v = 0; v <= ymaxPos + step * 1e-6; v += step) tks.push(+v.toFixed(8)); for (let v = step; v <= ymaxNeg + step * 1e-6; v += step) tks.push(-(+v.toFixed(8)));
  if (ymaxNeg) { ymaxNeg = Math.max(ymaxNeg, -Math.min(...tks)); }
  for (const v of tks) {
    if (S.showGrid) s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${y(v)}" y2="${y(v)}" stroke="${S.gridColor}" stroke-width="${S.lineWidth}"/>`;
    s += `<line x1="${m.l - 4}" x2="${m.l}" y1="${y(v)}" y2="${y(v)}" stroke="${S.axisColor}" stroke-width="${S.lineWidth}"/><text x="${m.l - 7}" y="${y(v) + S.tickSize * 0.35}" text-anchor="end" font-size="${S.tickSize}" fill="${S.textColor}">${S.measure === 'count' && layout !== 'cumulative' && spec.pooled !== 'mean' ? Math.abs(v).toLocaleString('en-US') : +Math.abs(v).toFixed(3)}</text>`;
  }
  const rect = (x, w, a, b, col, tip) => { const t = Math.min(y(a), y(b)), h = Math.abs(y(a) - y(b)); return h > 0 ? `<rect x="${x.toFixed(2)}" y="${t.toFixed(2)}" width="${Math.max(0.5, w).toFixed(2)}" height="${Math.max(0.5, h).toFixed(2)}" fill="${col}"><title>${esc(tip)}</title></rect>` : ''; };
  const errbar = (cx, v, sd, dir = 1) => { if (!(S.errorBars && sd > 0 && v >= 0)) return ''; const a = y(dir * v), b = y(dir * (v + sd)), c = y(dir * Math.max(0, v - sd)); const w = Math.min(6, bw / 6); return `<path d="M${cx} ${c}V${b}M${cx - w} ${b}H${cx + w}" stroke="${S.textColor}" stroke-width="${S.lineWidth}" fill="none"/>`; };
  const valTxt = (cx, v, top, dir = 1) => S.showValues && v >= 0.5 * Math.pow(10, -S.valueDecimals) ? `<text x="${cx.toFixed(2)}" y="${(dir > 0 ? y(top) - 4 : y(-top) + S.valueSize + 2).toFixed(2)}" text-anchor="middle" font-size="${S.valueSize}" fill="${S.textColor}">${fmtV(v)}</text>` : '';
  bins.forEach((b, i) => {
    const cx = m.l + slot * (i + 0.5);
    if (layout === 'grouped') {
      const w = bw / C.length;
      C.forEach((x, k) => { const xx = cx - bw / 2 + k * w; const sd = x.sd?.[i] || 0; s += rect(xx + 0.5, w - 1, 0, x.v[i], x.color, `${x.label} ${b.label}: ${fmtV(x.v[i])}`) + errbar(xx + w / 2, x.v[i], sd) + valTxt(xx + w / 2, x.v[i], x.v[i] + (S.errorBars ? sd : 0)); });
    } else if (layout === 'stacked') {
      let acc = 0; C.forEach(x => { s += rect(cx - bw / 2, bw, acc, acc + x.v[i], x.color, `${x.label} ${b.label}: ${fmtV(x.v[i])}`); acc += x.v[i]; });
      s += valTxt(cx, acc, acc);
    } else if (layout === 'mirror') {
      const d = C[0]; const sd0 = d.sd?.[i] || 0;
      s += rect(cx - bw / 2, bw, 0, -d.v[i], d.color, `${d.label} ${b.label}: ${fmtV(d.v[i])}`) + errbar(cx, d.v[i], sd0, -1) + valTxt(cx, d.v[i], d.v[i] + (S.errorBars ? sd0 : 0), -1);
      let acc = 0; C.slice(1).forEach(x => { s += rect(cx - bw / 2, bw, acc, acc + x.v[i], x.color, `${x.label} ${b.label}: ${fmtV(x.v[i])}`); acc += x.v[i]; });
      if (C.length === 2 && C[1].sd) s += errbar(cx, acc, C[1].sd[i]);
      s += valTxt(cx, acc, acc + (C.length === 2 && S.errorBars && C[1].sd ? C[1].sd[i] : 0));
    }
    const ly = m.t + ih + 6;
    s += rotate ? `<text transform="translate(${cx.toFixed(2)} ${ly + 4}) rotate(-45)" text-anchor="end" font-size="${S.tickSize}" fill="${S.textColor}">${esc(b.label)}</text>`
      : `<text x="${cx.toFixed(2)}" y="${ly + S.tickSize}" text-anchor="middle" font-size="${S.tickSize}" fill="${S.textColor}">${esc(b.label)}</text>`;
  });
  if (layout === 'cumulative') C.forEach(x => {
    const pts = x.v.map((v, i) => [m.l + slot * (i + 0.5), y(v)]);
    s += `<path d="M${m.l} ${y0}L${pts.map(p => p[0].toFixed(1) + ' ' + p[1].toFixed(1)).join('L')}" fill="none" stroke="${x.color}" stroke-width="${Math.max(2, S.lineWidth * 2)}" stroke-linejoin="round"/>`;
    pts.forEach((p, i) => { s += `<circle cx="${p[0]}" cy="${p[1]}" r="3.5" fill="${x.color}"><title>${esc(x.label)} ≤ ${esc(bins[i].label)}: ${x.v[i].toFixed(1)}%</title></circle>`; });
  });
  s += `<line x1="${m.l}" x2="${m.l}" y1="${m.t}" y2="${m.t + ih}" stroke="${S.axisColor}" stroke-width="${S.lineWidth}"/><line x1="${m.l}" x2="${m.l + iw}" y1="${y0}" y2="${y0}" stroke="${S.axisColor}" stroke-width="${S.lineWidth}"/>`;
  if (layout === 'mirror') s += `<text x="${m.l + 6}" y="${m.t + ih - 6}" font-size="${S.legendSize}" fill="${S.axisColor}">▼ ${esc(C[0].label)}</text>`;
  s += `<text x="${m.l + iw / 2}" y="${H - 8 - legendH}" text-anchor="middle" font-size="${S.axisTitleSize}" fill="${S.textColor}">${esc(spec.delOnly && /^indel length/i.test(S.xLabel) ? S.xLabel.replace(/^indel/i, 'Deletion') : S.xLabel)}</text>`;
  s += `<text transform="translate(${16 + S.axisTitleSize * 0.4} ${m.t + ih / 2}) rotate(-90)" text-anchor="middle" font-size="${S.axisTitleSize}" fill="${S.textColor}">${esc(yLabel)}</text>`;
  if (S.legend !== 'none' && C.length) {
    const sw = S.legendSize, gap = 14;
    const widths = C.map(x => sw + 6 + x.label.length * S.legendSize * 0.56);
    const total = widths.reduce((a, b) => a + b, 0) + gap * (C.length - 1);
    let lx = S.legend === 'bottom' ? m.l + (iw - total) / 2 : S.legend === 'top-left' ? m.l + 10 : m.l + iw - total - 6;
    const ly = S.legend === 'bottom' ? H - 6 : m.t + sw + 4;
    C.forEach((x, k) => { s += `<rect x="${lx}" y="${ly - sw + 1}" width="${sw}" height="${sw}" rx="2" fill="${x.color}"/><text x="${lx + sw + 6}" y="${ly}" font-size="${S.legendSize}" fill="${S.textColor}">${esc(x.label)}</text>`; lx += widths[k] + gap; });
  }
  return { svg: s + '</svg>' };
}

export function spectrumTSV(spec0, S) {
  const { bins, vals, sds, spec } = computeTable(spec0, { ...S, trimEmpty: false, layout: 'grouped' });
  const ser = spec.series;
  const meas = S.measure === 'count' ? (spec.pooled === 'mean' ? 'value' : 'count') : S.measure === 'pctEvents' ? 'pct_events' : `pct_${spec.unitLabel || 'reads'}`;
  const head = ['source', 'bin', 'min_bp', 'max_bp', ...ser.flatMap(se => [`${se.key}_${meas}`, ...(sds ? [`${se.key}_sd`] : [])])];
  const rows = [head];
  bins.forEach((b, i) => rows.push([spec.name, b.label, b.lo, b.hi === Infinity ? '' : b.hi, ...ser.flatMap((se, k) => [vals[k][i].toFixed(4), ...(sds ? [sds[k][i].toFixed(4)] : [])])]));
  rows.push([], ['# exact lengths (raw counts' + (spec.pooled === 'mean' ? ', summed over samples' : '') + ')'], ['source', 'type', 'length_bp', 'count']);
  for (const se of ser) for (const [len, c] of [...se.map].sort((a, b) => a[0] - b[0])) rows.push([spec.name, se.key, len, +c.toFixed ? +c.toFixed(4) : c]);
  return rows.map(r => r.join('\t')).join('\n') + '\n';
}

/** pool specs from several samples: 'reads' sums counts; 'mean' keeps per-sample specs for mean ± SD */
export function poolSpecs(specs, mode, name) {
  const series = specs[0].series.map(se => ({ ...se, map: new Map() }));
  for (const sp of specs) sp.series.forEach((se, k) => { for (const [len, c] of se.map) series[k].map.set(len, (series[k].map.get(len) || 0) + c); });
  const sum = f => specs.reduce((s, sp) => s + sp[f], 0);
  return { name, unitLabel: specs[0].unitLabel, total: sum('total'), edited: sum('edited'), events: sum('events'), series, delOnly: specs.every(x => x.delOnly), pooled: mode, perSample: mode === 'mean' ? specs : null };
}

/** Continuous histogram: bars touch, outlined, numeric size axis; several series are stacked. */
function renderHistogram(spec0, S) {
  const { bins, vals, sds, spec } = computeTable(spec0, S);
  const ser = spec.series.map((se, k) => ({ ...se, ...seriesStyle(S, se), v: vals[k], sd: sds?.[k] })).filter(x => x.show);
  const W = S.width, H = S.height, F = `font-family="${esc(S.font)}"`;
  const unit = spec.unitLabel || 'reads';
  const yLabel = S.yLabel === 'auto' ? (S.measure === 'count' ? (unit[0].toUpperCase() + unit.slice(1)) : S.measure === 'pctEvents' ? (spec.delOnly ? '% of deletions' : '% of all indel events') : `% of ${unit}`) + (spec.pooled === 'mean' && S.measure !== 'count' ? ' (mean ± SD)' : '') : S.yLabel;
  const subtitle = S.subtitle === 'auto' ? (spec.subtitle || spec.name) : S.subtitle;
  const title = delWord(spec, S.title || 'Indel size distribution');
  const legendH = S.legend === 'bottom' ? S.legendSize + 14 : 0;
  const top = 12 + (title ? S.titleSize + 6 : 0) + (subtitle ? S.subtitleSize + 6 : 0);
  const m = { l: 18 + S.axisTitleSize + 8 + S.tickSize * 3.4, r: 18, t: top + 10, b: S.tickSize + S.axisTitleSize + 26 + legendH };
  const iw = Math.max(50, W - m.l - m.r), ih = Math.max(50, H - m.t - m.b);
  // numeric x: equal bins use real bp; other binnings fall back to equal slots
  const numeric = bins.every(b => b.x0 !== undefined);
  const xmin = numeric ? bins[0].x0 : 0, xmax = numeric ? bins[bins.length - 1].x1 : bins.length;
  const X = v => m.l + (v - xmin) / (xmax - xmin) * iw;
  const bx0 = (b, i) => numeric ? X(b.x0) : X(i), bx1 = (b, i) => numeric ? X(b.x1) : X(i + 1);
  const stackTop = bins.map((_, i) => ser.reduce((a, x) => a + x.v[i], 0));
  let ymax = Math.max(...bins.map((_, i) => stackTop[i] + (ser.length === 1 && S.errorBars && ser[0].sd ? ser[0].sd[i] : 0)));
  ymax = ymax > 0 ? niceMax(ymax * (S.showValues ? 1.12 : 1.04)) : 1;
  const Y = v => m.t + ih - v / ymax * ih;
  const pct = S.measure !== 'count' || spec.pooled === 'mean';
  const fmtV = v => pct ? v.toFixed(S.valueDecimals) : Math.round(v).toLocaleString('en-US');
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" ${F}>`;
  if (!S.transparent) s += `<rect width="${W}" height="${H}" fill="${S.background}"/>`;
  let ty = 12;
  if (title) { ty += S.titleSize; s += `<text x="${m.l}" y="${ty}" font-size="${S.titleSize}" font-weight="700" fill="${S.textColor}">${esc(title)}</text>`; ty += 6; }
  if (subtitle) { ty += S.subtitleSize; s += `<text x="${m.l}" y="${ty}" font-size="${S.subtitleSize}" fill="${S.axisColor}">${esc(subtitle)}</text>`; }
  for (const v of ticks(ymax)) {
    if (S.showGrid) s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${Y(v)}" y2="${Y(v)}" stroke="${S.gridColor}" stroke-width="${S.lineWidth}"/>`;
    s += `<line x1="${m.l - 4}" x2="${m.l}" y1="${Y(v)}" y2="${Y(v)}" stroke="${S.axisColor}" stroke-width="${S.lineWidth}"/><text x="${m.l - 7}" y="${Y(v) + S.tickSize * 0.35}" text-anchor="end" font-size="${S.tickSize}" fill="${S.textColor}">${pct ? +v.toFixed(3) + (S.measure === 'count' ? '' : '%') : v.toLocaleString('en-US')}</text>`;
  }
  // bars (touching), stacked by series, with outline
  bins.forEach((b, i) => {
    const x0 = bx0(b, i), x1 = bx1(b, i); let acc = 0;
    for (const x of ser) {
      const v = x.v[i]; if (!(v > 0)) continue;
      s += `<rect x="${x0.toFixed(2)}" y="${Y(acc + v).toFixed(2)}" width="${Math.max(0.5, x1 - x0).toFixed(2)}" height="${Math.max(0.5, Y(acc) - Y(acc + v)).toFixed(2)}" fill="${x.color}" fill-opacity="${S.fillOpacity}" stroke="${S.barStrokeWidth > 0 ? S.barStroke : 'none'}" stroke-width="${S.barStrokeWidth}"><title>${esc(x.label)} ${esc(b.label)} bp: ${fmtV(v)}</title></rect>`;
      acc += v;
    }
    if (ser.length === 1 && S.errorBars && ser[0].sd && ser[0].sd[i] > 0 && acc > 0) { const cx = (x0 + x1) / 2, w = Math.min(5, (x1 - x0) / 4); s += `<path d="M${cx} ${Y(Math.max(0, acc - ser[0].sd[i]))}V${Y(acc + ser[0].sd[i])}M${cx - w} ${Y(acc + ser[0].sd[i])}H${cx + w}" stroke="${S.textColor}" stroke-width="${S.lineWidth}" fill="none"/>`; }
    if (S.showValues && acc >= 0.5 * Math.pow(10, -S.valueDecimals)) s += `<text x="${((x0 + x1) / 2).toFixed(2)}" y="${(Y(acc + (ser.length === 1 && S.errorBars && ser[0].sd ? ser[0].sd[i] : 0)) - 4).toFixed(2)}" text-anchor="middle" font-size="${S.valueSize}" fill="${S.textColor}">${fmtV(acc)}</text>`;
  });
  // axes
  s += `<line x1="${m.l}" x2="${m.l}" y1="${m.t}" y2="${m.t + ih}" stroke="${S.axisColor}" stroke-width="${S.lineWidth}"/><line x1="${m.l}" x2="${m.l + iw}" y1="${m.t + ih}" y2="${m.t + ih}" stroke="${S.axisColor}" stroke-width="${S.lineWidth}"/>`;
  // x ticks: numeric bp for equal bins, otherwise bin labels
  const ly = m.t + ih + 6 + S.tickSize;
  if (numeric) {
    const span = xmax - xmin;
    let step = +S.xTickStep > 0 ? +S.xTickStep : 0;
    if (!step) { const raw = span / Math.max(3, Math.floor(iw / (S.tickSize * 5.5))); const mg = Math.pow(10, Math.floor(Math.log10(raw))); step = [1, 2, 2.5, 5, 10].map(k => k * mg).find(k => k >= raw) || 10 * mg; const bw = bins[0].x1 - bins[0].x0; step = Math.max(bw, Math.round(step / bw) * bw); }
    const lastReal = bins.filter(b => !b.overflow).slice(-1)[0];
    const ovb = bins.find(b => b.overflow);
    for (let v = Math.ceil(xmin / step) * step; v <= (lastReal ? lastReal.x1 : xmax) + 1e-9; v += step) if (!(ovb && Math.abs(X(v) - (X(ovb.x0) + X(ovb.x1)) / 2) < S.tickSize * 2.2)) s += `<line x1="${X(v)}" x2="${X(v)}" y1="${m.t + ih}" y2="${m.t + ih + 4}" stroke="${S.axisColor}" stroke-width="${S.lineWidth}"/><text x="${X(v)}" y="${ly}" text-anchor="middle" font-size="${S.tickSize}" fill="${S.textColor}">${fmtLen(v)}</text>`;
    const ov = bins.find(b => b.overflow); if (ov) s += `<text x="${(X(ov.x0) + X(ov.x1)) / 2}" y="${ly}" text-anchor="middle" font-size="${S.tickSize}" fill="${S.axisColor}">${esc(ov.label)}</text>`;
  } else {
    const every = Math.max(1, Math.ceil(bins.length / Math.max(3, Math.floor(iw / (S.tickSize * 6)))));
    bins.forEach((b, i) => { if (i % every === 0) s += `<text x="${(bx0(b, i) + bx1(b, i)) / 2}" y="${ly}" text-anchor="middle" font-size="${S.tickSize}" fill="${S.textColor}">${esc(b.label)}</text>`; });
  }
  s += `<text x="${m.l + iw / 2}" y="${H - 8 - legendH}" text-anchor="middle" font-size="${S.axisTitleSize}" fill="${S.textColor}">${esc(spec.delOnly && /^indel length/i.test(S.xLabel) ? S.xLabel.replace(/^indel/i, 'Deletion') : S.xLabel)}</text>`;
  s += `<text transform="translate(${16 + S.axisTitleSize * 0.4} ${m.t + ih / 2}) rotate(-90)" text-anchor="middle" font-size="${S.axisTitleSize}" fill="${S.textColor}">${esc(yLabel)}</text>`;
  if (S.legend !== 'none' && ser.length > 1) {
    const sw = S.legendSize, gap = 14; const widths = ser.map(x => sw + 6 + x.label.length * S.legendSize * 0.56); const total = widths.reduce((a, b) => a + b, 0) + gap * (ser.length - 1);
    let lx = S.legend === 'bottom' ? m.l + (iw - total) / 2 : S.legend === 'top-left' ? m.l + 10 : m.l + iw - total - 6; const lyy = S.legend === 'bottom' ? H - 6 : m.t + sw + 4;
    ser.forEach((x, k) => { s += `<rect x="${lx}" y="${lyy - sw + 1}" width="${sw}" height="${sw}" fill="${x.color}" fill-opacity="${S.fillOpacity}" stroke="${S.barStroke}" stroke-width="${Math.min(1, S.barStrokeWidth)}"/><text x="${lx + sw + 6}" y="${lyy}" font-size="${S.legendSize}" fill="${S.textColor}">${esc(x.label)}</text>`; lx += widths[k] + gap; });
  }
  return { svg: s + '</svg>' };
}
