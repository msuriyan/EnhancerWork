// Per-sample bar plot: renderer + full editor (data filter, labels, colours, fonts, sizes, export).
import { svgToPngBlob } from './figureEditor.js';
import { openAnnotator } from './annotate.js';

const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const FONTS = ['Arial, Helvetica, sans-serif', 'Helvetica, Arial, sans-serif', '"Times New Roman", Times, serif', 'Georgia, serif', 'Calibri, Carlito, sans-serif', 'Verdana, sans-serif', '"Courier New", monospace'];
export const BAR_PALETTE = ['#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3', '#937860', '#da8bc3', '#8c8c8c', '#ccb974', '#64b5cd'];

export const BAR_DEFAULT = {
  title: '', ylabel: 'auto', xlabel: '',
  font: 'Arial, Helvetica, sans-serif', titleSize: 15, axisSize: 13, tickSize: 11, valueSize: 10, legendSize: 11,
  textColor: '#1a1a1a', axisColor: '#333333', gridColor: '#e8e8e8', background: '#ffffff', transparent: false,
  colorMode: 'single', color: '#4c72b0', groupColors: {}, sampleColors: {}, labels: {}, hidden: {},
  showNames: true, rotate: 'auto', showN: false,
  barWidth: 0.72, barStroke: '#333333', barStrokeWidth: 0.7, fillOpacity: 1,
  showValues: true, valueDecimals: 1, showGrid: true, ymax: 0, ystep: 0,
  legend: 'top-right', width: 0, height: 360, dpi: 300,
};

/** rows: [{id, name, value, n, group}] */
export function renderBars(rows0, S, meta = {}) {
  const rows = rows0.filter(r => !S.hidden[r.id]);
  const W = S.width || Math.max(360, 70 * rows.length + 150), H = S.height;
  const groups = [...new Set(rows.map(r => r.group || ''))];
  const gcol = g => S.groupColors[g] || BAR_PALETTE[groups.indexOf(g) % BAR_PALETTE.length];
  const colorOf = r => S.sampleColors[r.id] || (S.colorMode === 'group' ? gcol(r.group || '') : S.color);
  const label = r => S.labels[r.id] ?? r.name;
  const names = rows.map(label);
  const maxLab = Math.max(1, ...names.map(n => n.length));
  const slot0 = (W - 120) / Math.max(1, rows.length);
  const rot = !S.showNames ? 0 : S.rotate === 'auto' ? (maxLab * S.tickSize * 0.6 > slot0 ? 45 : 0) : +S.rotate;
  const nameBlock = !S.showNames ? 6 : rot === 90 ? maxLab * S.tickSize * 0.62 + 10 : rot === 45 ? maxLab * S.tickSize * 0.46 + 14 : S.tickSize + 10;
  const nBlock = S.showN ? S.tickSize + 4 : 0;
  const showLegend = S.legend !== 'none' && S.colorMode === 'group' && groups.filter(Boolean).length > 1;
  const legBottom = showLegend && S.legend === 'bottom' ? S.legendSize + 14 : 0;
  const m = { l: 22 + S.axisSize + S.tickSize * 3.2, r: 16, t: (S.title ? S.titleSize + 18 : 12), b: nameBlock + nBlock + (S.xlabel ? S.axisSize + 12 : 4) + legBottom + 6 };
  const iw = Math.max(40, W - m.l - m.r), ih = Math.max(40, H - m.t - m.b);
  const vmax = Math.max(0.0001, ...rows.map(r => r.value));
  let top = S.ymax > 0 ? S.ymax : (() => { const v = vmax * (S.showValues ? 1.15 : 1.05); const mg = Math.pow(10, Math.floor(Math.log10(v))); return Math.min(100, [1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10].map(k => k * mg).find(k => k >= v) || 10 * mg); })();
  const step = S.ystep > 0 ? S.ystep : (() => { const raw = top / 5; const mg = Math.pow(10, Math.floor(Math.log10(raw))); return [1, 2, 2.5, 5, 10].map(k => k * mg).find(k => k >= raw) || 10 * mg; })();
  const Y = v => m.t + ih - Math.min(v, top) / top * ih;
  const slot = iw / Math.max(1, rows.length), bw = slot * Math.min(0.98, Math.max(0.1, S.barWidth));
  const F = `font-family="${esc(S.font)}"`;
  let s = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" ${F}>`;
  if (!S.transparent) s += `<rect width="${W}" height="${H}" fill="${S.background}"/>`;
  if (S.title) s += `<text x="${m.l}" y="${S.titleSize + 4}" font-size="${S.titleSize}" font-weight="700" fill="${S.textColor}">${esc(S.title)}</text>`;
  for (let v = 0; v <= top + step * 1e-6; v += step) {
    if (S.showGrid) s += `<line x1="${m.l}" x2="${m.l + iw}" y1="${Y(v)}" y2="${Y(v)}" stroke="${S.gridColor}"/>`;
    s += `<line x1="${m.l - 4}" x2="${m.l}" y1="${Y(v)}" y2="${Y(v)}" stroke="${S.axisColor}"/><text x="${m.l - 7}" y="${Y(v) + S.tickSize * 0.35}" text-anchor="end" font-size="${S.tickSize}" fill="${S.textColor}">${+v.toFixed(4)}</text>`;
  }
  rows.forEach((r, i) => {
    const cx = m.l + slot * (i + 0.5), y = Y(Math.max(0, r.value));
    s += `<rect x="${(cx - bw / 2).toFixed(2)}" y="${y.toFixed(2)}" width="${bw.toFixed(2)}" height="${Math.max(0, m.t + ih - y).toFixed(2)}" fill="${colorOf(r)}" fill-opacity="${S.fillOpacity}" stroke="${S.barStrokeWidth > 0 ? S.barStroke : 'none'}" stroke-width="${S.barStrokeWidth}"><title>${esc(label(r))}: ${r.value.toFixed(2)}% (n = ${r.n.toLocaleString('en-US')} reads)</title></rect>`;
    if (S.showValues) s += `<text x="${cx.toFixed(2)}" y="${(y - 4).toFixed(2)}" text-anchor="middle" font-size="${S.valueSize}" fill="${S.textColor}">${r.value.toFixed(S.valueDecimals)}</text>`;
    if (S.showNames) {
      const ly = m.t + ih + 8;
      s += rot ? `<text transform="translate(${cx.toFixed(2)} ${ly}) rotate(-${rot})" text-anchor="end" dominant-baseline="${rot === 90 ? 'middle' : 'auto'}" font-size="${S.tickSize}" fill="${S.textColor}">${esc(label(r))}</text>`
        : `<text x="${cx.toFixed(2)}" y="${ly + S.tickSize}" text-anchor="middle" font-size="${S.tickSize}" fill="${S.textColor}">${esc(label(r))}</text>`;
    }
    if (S.showN) s += `<text x="${cx.toFixed(2)}" y="${m.t + ih + nameBlock + S.tickSize}" text-anchor="middle" font-size="${S.tickSize * 0.85}" fill="${S.axisColor}">n=${r.n.toLocaleString('en-US')}</text>`;
  });
  s += `<line x1="${m.l}" x2="${m.l}" y1="${m.t}" y2="${m.t + ih}" stroke="${S.axisColor}"/><line x1="${m.l}" x2="${m.l + iw}" y1="${m.t + ih}" y2="${m.t + ih}" stroke="${S.axisColor}"/>`;
  const yl = S.ylabel === 'auto' ? (meta.ylabel || '% of reads') : S.ylabel;
  if (yl) s += `<text transform="translate(${14 + S.axisSize * 0.35} ${m.t + ih / 2}) rotate(-90)" text-anchor="middle" font-size="${S.axisSize}" fill="${S.textColor}">${esc(yl)}</text>`;
  if (S.xlabel) s += `<text x="${m.l + iw / 2}" y="${H - 6 - legBottom}" text-anchor="middle" font-size="${S.axisSize}" fill="${S.textColor}">${esc(S.xlabel)}</text>`;
  if (showLegend) {
    const gs = groups.filter(Boolean); const sw = S.legendSize, gap = 14;
    const widths = gs.map(g => sw + 6 + g.length * S.legendSize * 0.58); const total = widths.reduce((a, b) => a + b, 0) + gap * (gs.length - 1);
    let lx = S.legend === 'bottom' ? m.l + (iw - total) / 2 : S.legend === 'top-left' ? m.l + 10 : m.l + iw - total - 6;
    const ly = S.legend === 'bottom' ? H - 6 : m.t + sw + 4;
    gs.forEach((g, k) => { s += `<rect x="${lx}" y="${ly - sw + 1}" width="${sw}" height="${sw}" fill="${gcol(g)}" fill-opacity="${S.fillOpacity}" stroke="${S.barStroke}" stroke-width="${Math.min(1, S.barStrokeWidth)}"/><text x="${lx + sw + 6}" y="${ly}" font-size="${S.legendSize}" fill="${S.textColor}">${esc(g)}</text>`; lx += widths[k] + gap; });
  }
  return s + '</svg>';
}

function loadStyle() { try { return { ...BAR_DEFAULT, ...(JSON.parse(localStorage.getItem('genoview.fig.bars') || 'null') || {}) }; } catch { return { ...BAR_DEFAULT }; } }
function saveStyle(S) { try { localStorage.setItem('genoview.fig.bars', JSON.stringify(S)); } catch {} }
export function barStyle() { return loadStyle(); }

/**
 * opts: { metrics:[[key,label,axis]], data:{metric,minDel,maxDel}, getRows(data) -> rows, onData(data), filename }
 */
export function openBarEditor(opts) {
  document.getElementById('barEditor')?.remove();
  const S = loadStyle(); const D = { ...opts.data };
  const wrap = document.createElement('div'); wrap.id = 'barEditor'; wrap.className = 'figmodal';
  wrap.innerHTML = `<div class="figbox"><div class="fighead"><b>Per-sample bar plot — figure editor</b><span class="muted" style="margin-left:8px">the preview is exactly what is exported</span><span style="flex:1"></span><button id="beReset">Reset style</button><button id="beClose" class="primary">Close</button></div>
    <div class="figbody"><div class="figctl" id="beCtl"></div><div class="figprev"><div id="bePrev" style="flex:1;overflow:auto;padding:16px;background:repeating-conic-gradient(#f3f3f0 0% 25%,#fff 0% 50%) 50%/20px 20px;display:flex;align-items:flex-start;justify-content:center"></div>
    <div class="figexp"><label class="status">PNG <select id="beDpi">${[150, 300, 600, 1200].map(d => `<option ${d === S.dpi ? 'selected' : ''} value="${d}">${d} dpi</option>`).join('')}</select></label>
      <button id="bePng" class="primary">Download PNG</button><button id="beSvg">SVG</button><button id="bePdf">PDF</button><button id="beTsv">Data TSV</button><button id="beAnn" title="Add your own text boxes / arrows, then export">Annotate (text boxes)…</button><span class="status" id="beInfo"></span></div></div></div></div>`;
  document.body.appendChild(wrap);
  const $ = q => wrap.querySelector(q);
  let rows = opts.getRows(D);
  const metricAxis = () => (opts.metrics.find(m => m[0] === D.metric) || [])[2] || '% of reads';
  const yAuto = () => D.metric === 'readsWithDelPct' && (D.minDel > 1 || D.maxDel) ? `Reads with deletion ${D.maxDel ? `${D.minDel || 1}–${D.maxDel}` : '≥ ' + D.minDel} bp (%)` : metricAxis();
  const F = (l, h, t = '') => `<label class="fe" title="${esc(t)}"><span>${l}</span>${h}</label>`;
  const num = (k, a, b, st = 1) => `<input type="number" data-k="${k}" value="${S[k]}" min="${a}" max="${b}" step="${st}">`;
  const col = k => `<input type="color" data-k="${k}" value="${S[k]}">`;
  const chk = k => `<input type="checkbox" data-k="${k}" ${S[k] ? 'checked' : ''}>`;
  const opt = (k, list) => `<select data-k="${k}">${list.map(([v, l]) => `<option value="${esc(v)}" ${String(S[k]) === String(v) ? 'selected' : ''}>${esc(l)}</option>`).join('')}</select>`;
  function build() {
    const groups = [...new Set(rows.map(r => r.group || ''))].filter(Boolean);
    $('#beCtl').innerHTML = `
      <fieldset><legend>Data</legend>
        ${F('Value', `<select id="beMetric">${opts.metrics.map(([v, l]) => `<option value="${v}" ${v === D.metric ? 'selected' : ''}>${esc(l)}</option>`).join('')}</select>`)}
        ${D.metric === 'readsWithDelPct' ? `<div class="fe-grid">${F('Deletion ≥ bp', `<input type="number" id="beMin" value="${D.minDel || 1}" min="1" step="1">`, 'Count only reads whose deletion is at least this long')}${F('and ≤ bp (0 = any)', `<input type="number" id="beMax" value="${D.maxDel || 0}" min="0" step="1">`)}</div>
        <div class="fe-presets">${[[1, 0, 'all'], [2, 0, '≥2'], [10, 0, '≥10'], [50, 0, '≥50'], [100, 0, '≥100'], [500, 0, '≥500'], [1, 49, '1–49']].map(([a, b, l]) => `<button data-pre="${a},${b}">${l}</button>`).join('')}</div>` : ''}
        ${F('Order', `<select id="beOrder">${[['input', 'as loaded'], ['value', 'high → low'], ['valueAsc', 'low → high'], ['name', 'by name'], ['group', 'by group']].map(([v, l]) => `<option value="${v}" ${v === (D.order || 'input') ? 'selected' : ''}>${l}</option>`).join('')}</select>`)}
      </fieldset>
      <fieldset><legend>Samples (x-axis)</legend>
        ${F('Show sample names', chk('showNames'), 'Untick to hide the x-axis names (add your own with Annotate)')}
        ${F('Show n (reads) under bars', chk('showN'))}
        ${F('Rotate names', opt('rotate', [['auto', 'auto'], ['0', '0°'], ['45', '45°'], ['90', '90°']]))}
        <div class="hint" style="margin:4px 0">Per bar: show · colour · label</div>
        ${rows.map(r => `<div class="fe-series"><input type="checkbox" data-hide="${esc(r.id)}" ${S.hidden[r.id] ? '' : 'checked'} title="Show this bar"><input type="color" data-scol="${esc(r.id)}" value="${S.sampleColors[r.id] || (S.colorMode === 'group' ? (S.groupColors[r.group] || BAR_PALETTE[[...new Set(rows.map(x => x.group || ''))].indexOf(r.group || '') % BAR_PALETTE.length]) : S.color)}" title="Colour of this bar (overrides)"><input type="text" data-lab="${esc(r.id)}" value="${esc(S.labels[r.id] ?? r.name)}" title="Label on the x-axis"></div>`).join('')}
        <div class="rulebar" style="margin:4px 0"><button id="beClrS">Reset bar colours</button><button id="beClrL">Reset labels</button></div>
      </fieldset>
      <fieldset><legend>Colours</legend>
        ${F('Colour bars', opt('colorMode', [['single', 'one colour'], ['group', 'by group (track-list group)']]))}
        ${S.colorMode === 'single' ? F('Bar colour', col('color')) : groups.map((g, k) => `<div class="fe-series"><input type="color" data-gcol="${esc(g)}" value="${S.groupColors[g] || BAR_PALETTE[k % BAR_PALETTE.length]}"><span style="font-size:12px">${esc(g)}</span></div>`).join('') || '<div class="hint">Set groups in the track list (BAM → group).</div>'}
        <div class="fe-grid">${F('Outline', col('barStroke'))}${F('Outline w', num('barStrokeWidth', 0, 5, 0.25))}${F('Fill opacity', num('fillOpacity', 0.05, 1, 0.05))}${F('Bar width', num('barWidth', 0.1, 0.98, 0.02))}${F('Text', col('textColor'))}${F('Axes', col('axisColor'))}${F('Grid', col('gridColor'))}${F('Background', col('background'))}</div>
        ${F('Transparent background', chk('transparent'))}${F('Grid lines', chk('showGrid'))}
        ${F('Legend', opt('legend', [['top-right', 'top right'], ['top-left', 'top left'], ['bottom', 'bottom'], ['none', 'none']]))}
      </fieldset>
      <fieldset><legend>Text &amp; fonts</legend>
        ${F('Title', `<input type="text" data-k="title" value="${esc(S.title)}">`)}
        ${F('Y-axis title', `<input type="text" data-k="ylabel" value="${esc(S.ylabel)}" placeholder="auto">`, 'Type auto to follow the value / deletion filter')}
        ${F('X-axis title', `<input type="text" data-k="xlabel" value="${esc(S.xlabel)}">`)}
        ${F('Font', `<select data-k="font">${[...new Set([S.font, ...FONTS])].map(f => `<option value='${esc(f)}' ${f === S.font ? 'selected' : ''}>${esc(f.split(',')[0].replace(/"/g, ''))}</option>`).join('')}</select>`)}
        ${F('…or installed font', `<input type="text" id="beFontFree" placeholder="e.g. Myriad Pro">`)}
        <div class="fe-grid">${F('Title', num('titleSize', 4, 72))}${F('Axis titles', num('axisSize', 4, 72))}${F('Tick labels', num('tickSize', 4, 72))}${F('Bar values', num('valueSize', 4, 72))}${F('Legend', num('legendSize', 4, 72))}${F('Decimals', num('valueDecimals', 0, 4))}</div>
        ${F('Values on bars', chk('showValues'))}
      </fieldset>
      <fieldset><legend>Axis &amp; size</legend>
        <div class="fe-grid">${F('Y max (0=auto)', num('ymax', 0, 100, 5))}${F('Y step (0=auto)', num('ystep', 0, 100, 1))}${F('Width px (0=auto)', num('width', 0, 4000, 10))}${F('Height px', num('height', 150, 4000, 10))}</div>
        <div class="fe-presets">${[['Single column', 340, 280], ['1.5 column', 520, 340], ['Double column', 720, 400]].map(([l, w, h]) => `<button data-w="${w}" data-h="${h}">${l}</button>`).join('')}</div>
      </fieldset>`;
    const ctl = $('#beCtl');
    ctl.querySelectorAll('[data-k]').forEach(el => el.addEventListener(el.type === 'text' || el.type === 'number' || el.type === 'color' ? 'input' : 'change', () => {
      const k = el.dataset.k; S[k] = el.type === 'checkbox' ? el.checked : el.type === 'number' ? +el.value : el.value;
      if (k === 'colorMode') build(); draw();
    }));
    ctl.querySelectorAll('[data-hide]').forEach(el => el.onchange = () => { S.hidden[el.dataset.hide] = !el.checked; draw(); });
    ctl.querySelectorAll('[data-scol]').forEach(el => el.oninput = () => { S.sampleColors[el.dataset.scol] = el.value; draw(); });
    ctl.querySelectorAll('[data-lab]').forEach(el => el.oninput = () => { S.labels[el.dataset.lab] = el.value; draw(); });
    ctl.querySelectorAll('[data-gcol]').forEach(el => el.oninput = () => { S.groupColors[el.dataset.gcol] = el.value; draw(); });
    ctl.querySelectorAll('[data-w]').forEach(b => b.onclick = () => { S.width = +b.dataset.w; S.height = +b.dataset.h; build(); draw(); });
    $('#beClrS').onclick = () => { S.sampleColors = {}; build(); draw(); };
    $('#beClrL').onclick = () => { S.labels = {}; build(); draw(); };
    const ff = $('#beFontFree'); ff.onchange = () => { if (ff.value.trim()) { S.font = `"${ff.value.trim()}", Arial, sans-serif`; build(); draw(); } };
    const reData = () => { rows = opts.getRows(D); opts.onData?.(D); build(); draw(); };
    $('#beMetric').onchange = e => { D.metric = e.target.value; reData(); };
    $('#beOrder').onchange = e => { D.order = e.target.value; reData(); };
    if ($('#beMin')) { $('#beMin').onchange = e => { D.minDel = Math.max(1, +e.target.value || 1); reData(); }; $('#beMax').onchange = e => { D.maxDel = Math.max(0, +e.target.value || 0); reData(); }; }
    ctl.querySelectorAll('[data-pre]').forEach(b => b.onclick = () => { const [a, c] = b.dataset.pre.split(',').map(Number); D.minDel = a; D.maxDel = c; reData(); });
  }
  const svg = () => renderBars(rows, S, { ylabel: yAuto() });
  function draw() {
    const out = svg(); $('#bePrev').innerHTML = out;
    const w = +(/width="([\d.]+)"/.exec(out)?.[1] || 600), h = +(/height="([\d.]+)"/.exec(out)?.[1] || 400); S.dpi = +$('#beDpi').value;
    $('#beInfo').textContent = `${Math.round(w)}×${Math.round(h)} px → PNG ${Math.round(w * S.dpi / 96)}×${Math.round(h * S.dpi / 96)} px at ${S.dpi} dpi`;
    saveStyle(S);
  }
  const fname = ext => `${opts.filename || 'per_sample_bars'}.${ext}`;
  const dl = (n, b) => { const a = document.createElement('a'); a.href = URL.createObjectURL(b); a.download = n; document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500); };
  const dims = s => [+(/width="([\d.]+)"/.exec(s)[1]), +(/height="([\d.]+)"/.exec(s)[1])];
  $('#beDpi').onchange = draw;
  $('#bePng').onclick = async () => { const s = svg(); const [w, h] = dims(s); dl(fname('png'), await svgToPngBlob(s, w, h, S.dpi, S.transparent)); };
  $('#beSvg').onclick = () => dl(fname('svg'), new Blob(['<?xml version="1.0" encoding="UTF-8"?>\n' + svg()], { type: 'image/svg+xml' }));
  $('#bePdf').onclick = () => { const s = svg(); const [w, h] = dims(s); const win = window.open('', '_blank'); if (!win) return alert('Allow pop-ups to export PDF.'); const mw = w / 96 * 25.4, mh = h / 96 * 25.4; win.document.write(`<!DOCTYPE html><html><head><title>${esc(fname('pdf'))}</title><style>@page{size:${mw}mm ${mh}mm;margin:0}html,body{margin:0}svg{display:block;width:${mw}mm;height:${mh}mm}</style></head><body>${s}<script>setTimeout(()=>print(),300)<\/script></body></html>`); win.document.close(); };
  $('#beTsv').onclick = () => dl(fname('tsv'), new Blob(['sample\tlabel\tgroup\tspanning_reads\t' + yAuto().replace(/\s+/g, '_') + '\n' + rows.filter(r => !S.hidden[r.id]).map(r => `${r.name}\t${S.labels[r.id] ?? r.name}\t${r.group || ''}\t${r.n}\t${r.value.toFixed(4)}`).join('\n') + '\n'], { type: 'text/tab-separated-values' }));
  $('#beAnn').onclick = () => openAnnotator({ svg: svg(), key: 'bars', title: 'Annotate · per-sample bar plot', filename: opts.filename || 'per_sample_bars' });
  $('#beClose').onclick = () => wrap.remove();
  $('#beReset').onclick = () => { Object.keys(S).forEach(k => delete S[k]); Object.assign(S, JSON.parse(JSON.stringify(BAR_DEFAULT))); build(); draw(); };
  build(); draw();
}
