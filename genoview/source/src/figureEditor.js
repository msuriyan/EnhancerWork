// Modal figure editor for the indel/SV size spectrum: live preview, every text/font/colour/size editable,
// export PNG (with DPI metadata), SVG, PDF (via print), TSV. Style is remembered per figure type.
import { renderSpectrum, spectrumTSV, DEFAULT_STYLE, BIN_PRESETS } from './spectrum.js';
import { openAnnotator } from './annotate.js';

const FONTS = ['Arial, Helvetica, sans-serif', 'Helvetica, Arial, sans-serif', '"Times New Roman", Times, serif', 'Georgia, serif', 'Calibri, Carlito, sans-serif', 'Verdana, sans-serif', '"Trebuchet MS", sans-serif', '"Courier New", monospace', 'system-ui, sans-serif'];
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

function loadStyle(key, defaults) {
  try { const s = JSON.parse(localStorage.getItem('genoview.fig.' + key) || 'null'); if (s && !('binWidth' in s) && !defaults.layout) { s.layout = 'histogram'; s.bins = 'equal'; s.showValues = false; } if (s) return { ...DEFAULT_STYLE, ...defaults, ...s, series: { ...(defaults.series || {}), ...(s.series || {}) } }; } catch {}
  return { ...DEFAULT_STYLE, ...defaults };
}
function saveStyle(key, S) { try { localStorage.setItem('genoview.fig.' + key, JSON.stringify(S)); } catch {} }

/**
 * opts: { key, title, selects: [{ id, label, options:[{value,label}], value }], getSpec(selValues) -> spec, defaults, filename(selValues) }
 */
export function openFigureEditor(opts) {
  document.getElementById('figEditor')?.remove();
  const S = loadStyle(opts.key, opts.defaults || {});
  const sel = Object.fromEntries(opts.selects.map(s => [s.id, s.value ?? s.options[0]?.value]));
  const wrap = document.createElement('div');
  wrap.id = 'figEditor'; wrap.className = 'figmodal';
  wrap.innerHTML = `<div class="figbox">
    <div class="fighead"><b>${esc(opts.title)}</b><span class="muted" style="margin-left:8px">edit anything — the preview is exactly what is exported</span><span style="flex:1"></span>
      <button id="feReset" title="Restore default style">Reset style</button><button id="feClose" class="primary">Close</button></div>
    <div class="figbody"><div class="figctl" id="feCtl"></div><div class="figprev"><div id="fePrev"></div>
      <div class="figexp">
        <label class="status">PNG <select id="feDpi">${[150, 300, 600, 1200].map(d => `<option ${d === S.dpi ? 'selected' : ''} value="${d}">${d} dpi</option>`).join('')}</select></label>
        <button id="fePng" class="primary">Download PNG</button><button id="feSvg">SVG (vector)</button><button id="fePdf" title="Opens the print dialog — choose “Save as PDF”">PDF</button><button id="feTsv">Data TSV</button><button id="feAnn" title="Add free text boxes, arrows and boxes, then export at any dpi">Annotate…</button><button id="feCopy" title="Copy the SVG markup">Copy SVG</button>
        <span class="status" id="feSize"></span></div></div></div></div>`;
  document.body.appendChild(wrap);
  const $ = q => wrap.querySelector(q);

  let spec = opts.getSpec(sel);
  const ctl = $('#feCtl');
  const field = (label, html, title = '') => `<label class="fe" title="${esc(title)}"><span>${label}</span>${html}</label>`;
  const txt = (k, ph = '') => `<input type="text" data-k="${k}" value="${esc(S[k] ?? '')}" placeholder="${esc(ph)}">`;
  const num = (k, min, max, step = 1) => `<input type="number" data-k="${k}" value="${S[k]}" min="${min}" max="${max}" step="${step}">`;
  const col = k => `<input type="color" data-k="${k}" value="${S[k]}">`;
  const chk = k => `<input type="checkbox" data-k="${k}" ${S[k] ? 'checked' : ''}>`;
  const opt = (k, list) => `<select data-k="${k}">${list.map(([v, l]) => `<option value="${esc(v)}" ${String(S[k]) === String(v) ? 'selected' : ''}>${esc(l)}</option>`).join('')}</select>`;

  function buildControls() {
    const shownSeries = S.mergeSeries && spec.series.length > 1 ? [{ key: 'all', label: 'All indels', color: '#c9c9c9' }] : spec.series;
    const seriesRows = shownSeries.map(se => {
      const st = { label: S.series[se.key]?.label ?? se.label, color: S.series[se.key]?.color ?? se.color, show: S.series[se.key]?.show ?? true };
      return `<div class="fe-series"><input type="checkbox" data-s="${se.key}" data-p="show" ${st.show ? 'checked' : ''} title="Show"><input type="color" data-s="${se.key}" data-p="color" value="${st.color}"><input type="text" data-s="${se.key}" data-p="label" value="${esc(st.label)}"></div>`;
    }).join('');
    ctl.innerHTML = `
      <fieldset><legend>Data</legend>
        ${opts.selects.map(s => field(s.label, `<select data-sel="${s.id}">${s.options.map(o => `<option value="${esc(o.value)}" ${o.value === sel[s.id] ? 'selected' : ''}>${esc(o.label)}</option>`).join('')}</select>`)).join('')}
        ${field('Y value', opt('measure', [['pct', `% of ${spec.unitLabel || 'reads'}`], ['count', spec.pooled === 'mean' ? 'mean count' : 'count'], ['pctEvents', spec.delOnly ? '% of deletions (sums to 100%)' : '% of all indel events (sums to 100%)']]))}
        ${field('Chart', opt('layout', [['histogram', 'Histogram (continuous, bars touch)'], ['grouped', 'Grouped bars'], ['stacked', 'Stacked bars'], ['mirror', 'Mirrored (first series down)'], ['cumulative', 'Cumulative (CDF)']]))}
        ${field('Combine deletions + insertions', chk('mergeSeries'), 'One series of all indel sizes (series name/colour under Series → all)')}
        ${field('Size bins', opt('bins', [...Object.entries(BIN_PRESETS).map(([k, v]) => [k, v.label]), ['custom', 'Custom cut-points…']]))}
        ${S.bins === 'equal' ? `<div class="fe-grid">${field('Bin width bp', num('binWidth', 1, 100000))}${field('Max bp (0=auto)', num('binMax', 0, 10000000, 50))}</div>${field('Extra bin for longer', chk('overflowBin'))}` : ''}
        ${field('Cut-points (bp)', `<input type="text" data-k="customEdges" value="${esc(S.customEdges)}" placeholder="20, 50, 100, 500, 1000, 2000">`, 'Upper bound of each bin; a final >last bin is added. Used when Size bins = Custom.')}
        ${field('Drop empty bins at the end', chk('trimEmpty'))}
        ${spec.pooled === 'mean' ? field('Error bars (± SD)', chk('errorBars')) : ''}
      </fieldset>
      <fieldset><legend>Series</legend>${seriesRows}</fieldset>
      <fieldset><legend>Text</legend>
        ${field('Title', txt('title', 'Indel size spectrum'))}
        ${field('Subtitle', txt('subtitle', 'auto'), 'Type “auto” for sample / n, or leave empty for none')}
        ${field('X-axis title', txt('xLabel'))}
        ${field('Y-axis title', txt('yLabel', 'auto'), 'Type “auto” to follow the Y value')}
      </fieldset>
      <fieldset><legend>Fonts</legend>
        ${field('Font', `<select data-k="font">${[...new Set([S.font, ...FONTS])].map(f => `<option value='${esc(f)}' ${f === S.font ? 'selected' : ''}>${esc(f.split(',')[0].replace(/"/g, ''))}</option>`).join('')}</select>`)}
        ${field('…or any installed font', `<input type="text" id="feFontFree" placeholder="e.g. Myriad Pro">`)}
        <div class="fe-grid">${[['titleSize', 'Title'], ['subtitleSize', 'Subtitle'], ['axisTitleSize', 'Axis titles'], ['tickSize', 'Tick labels'], ['valueSize', 'Bar values'], ['legendSize', 'Legend']].map(([k, l]) => field(l, num(k, 4, 72))).join('')}</div>
      </fieldset>
      <fieldset><legend>Appearance</legend>
        <div class="fe-grid">${field('Text', col('textColor'))}${field('Axes', col('axisColor'))}${field('Grid', col('gridColor'))}${field('Background', col('background'))}</div>
        ${field('Transparent background', chk('transparent'))}
        ${field('Grid lines', chk('showGrid'))}
        ${field('Values on bars', chk('showValues'))}${field('Decimals', num('valueDecimals', 0, 4))}
        ${field('Legend', opt('legend', [['top-right', 'Top right'], ['top-left', 'Top left'], ['bottom', 'Bottom'], ['none', 'None']]))}
        ${S.layout === 'histogram' ? `<div class="fe-grid">${field('Bar outline', col('barStroke'))}${field('Outline w', num('barStrokeWidth', 0, 5, 0.25))}${field('Fill opacity', num('fillOpacity', 0.05, 1, 0.05))}${field('X tick every bp (0=auto)', num('xTickStep', 0, 1000000, 50))}</div>` : ''}
        ${field('Rotate x labels', opt('rotateTicks', [['auto', 'Auto'], ['yes', 'Yes (45°)'], ['no', 'No']]))}
        ${field('Bar width', num('barWidth', 0.1, 0.98, 0.02))}${field('Line width', num('lineWidth', 0.5, 4, 0.25))}
      </fieldset>
      <fieldset><legend>Size</legend>
        <div class="fe-grid">${field('Width px', num('width', 240, 4000, 10))}${field('Height px', num('height', 180, 4000, 10))}</div>
        <div class="fe-presets">${[['Single column', 340, 260], ['1.5 column', 520, 340], ['Double column', 720, 420], ['Slide 16:9', 1280, 720]].map(([l, w, h]) => `<button data-w="${w}" data-h="${h}">${l}</button>`).join('')}</div>
      </fieldset>`;
    ctl.querySelectorAll('[data-k]').forEach(el => el.addEventListener(el.type === 'text' || el.type === 'number' ? 'input' : 'change', () => {
      const k = el.dataset.k; let v = el.type === 'checkbox' ? el.checked : el.type === 'number' ? +el.value : el.value;
      S[k] = v; if (k === 'customEdges' && S.bins !== 'custom') { S.bins = 'custom'; ctl.querySelector('[data-k="bins"]').value = 'custom'; }
      if (k === 'layout' && v === 'histogram' && S.bins !== 'equal') { S.bins = 'equal'; }
      if (['layout', 'bins', 'mergeSeries'].includes(k)) { buildControls(); }
      draw();
    }));
    ctl.querySelectorAll('[data-s]').forEach(el => el.addEventListener(el.type === 'text' ? 'input' : 'change', () => {
      const k = el.dataset.s; S.series[k] = { ...(S.series[k] || {}), [el.dataset.p]: el.type === 'checkbox' ? el.checked : el.value }; draw();
    }));
    ctl.querySelectorAll('[data-sel]').forEach(el => el.addEventListener('change', () => { sel[el.dataset.sel] = el.value; spec = opts.getSpec(sel); buildControls(); draw(); }));
    ctl.querySelectorAll('[data-w]').forEach(b => b.onclick = () => { S.width = +b.dataset.w; S.height = +b.dataset.h; buildControls(); draw(); });
    const ff = ctl.querySelector('#feFontFree'); ff.onchange = () => { if (ff.value.trim()) { S.font = `"${ff.value.trim()}", Arial, sans-serif`; buildControls(); draw(); } };
  }
  function draw() {
    S.dpi = +$('#feDpi').value;
    const { svg } = renderSpectrum(spec, S);
    $('#fePrev').innerHTML = svg;
    const px = Math.round(S.width * S.dpi / 96), py = Math.round(S.height * S.dpi / 96);
    $('#feSize').textContent = `${S.width}×${S.height} px canvas · PNG ${px}×${py} px = ${(S.width / 96 * 25.4).toFixed(0)}×${(S.height / 96 * 25.4).toFixed(0)} mm at ${S.dpi} dpi`;
    saveStyle(opts.key, S);
  }
  const fname = ext => (opts.filename ? opts.filename(sel) : 'indel_spectrum') + '.' + ext;
  $('#feDpi').onchange = draw;
  $('#feAnn').onclick = () => openAnnotator({ svg: renderSpectrum(spec, S).svg, key: 'spectrum-' + opts.key, title: 'Annotate & export · size spectrum', filename: opts.filename ? opts.filename(sel) : 'indel_spectrum' });
  $('#feClose').onclick = () => wrap.remove();
  wrap.addEventListener('mousedown', e => { if (e.target === wrap) wrap.remove(); });
  $('#feReset').onclick = () => { Object.keys(S).forEach(k => delete S[k]); Object.assign(S, { ...DEFAULT_STYLE, ...(opts.defaults || {}), series: { ...((opts.defaults || {}).series || {}) } }); buildControls(); draw(); };
  $('#feSvg').onclick = () => dl(fname('svg'), new Blob(['<?xml version="1.0" encoding="UTF-8"?>\n' + renderSpectrum(spec, S).svg], { type: 'image/svg+xml' }));
  $('#feTsv').onclick = () => dl(fname('tsv'), new Blob([spectrumTSV(spec, S)], { type: 'text/tab-separated-values' }));
  $('#feCopy').onclick = async () => { try { await navigator.clipboard.writeText(renderSpectrum(spec, S).svg); $('#feCopy').textContent = 'Copied ✓'; setTimeout(() => $('#feCopy').textContent = 'Copy SVG', 1500); } catch { alert('Clipboard not available here — use SVG download.'); } };
  $('#fePng').onclick = () => svgToPngBlob(renderSpectrum(spec, S).svg, S.width, S.height, S.dpi, S.transparent).then(b => dl(fname('png'), b));
  $('#fePdf').onclick = () => {
    const w = window.open('', '_blank'); if (!w) return alert('Allow pop-ups to export PDF.');
    const mmW = S.width / 96 * 25.4, mmH = S.height / 96 * 25.4;
    w.document.write(`<!DOCTYPE html><html><head><title>${esc(fname('pdf'))}</title><style>@page{size:${mmW}mm ${mmH}mm;margin:0}html,body{margin:0}svg{display:block;width:${mmW}mm;height:${mmH}mm}</style></head><body>${renderSpectrum(spec, S).svg}<script>setTimeout(()=>{print();},300)<\/script></body></html>`);
    w.document.close();
  };
  buildControls(); draw();
}

function dl(name, blob) { const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500); }

// PNG with a pHYs chunk so journals / Word / Illustrator read the intended DPI
export function svgToPngBlob(svg, w, h, dpi, transparent) {
  return new Promise((resolve, reject) => {
    const scale = dpi / 96; const img = new Image();
    img.onload = () => {
      const cv = document.createElement('canvas'); cv.width = Math.round(w * scale); cv.height = Math.round(h * scale);
      const ctx = cv.getContext('2d'); if (!transparent) { ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, cv.width, cv.height); }
      ctx.drawImage(img, 0, 0, cv.width, cv.height);
      cv.toBlob(async b => { try { resolve(await setPngDpi(b, dpi)); } catch { resolve(b); } }, 'image/png');
    };
    img.onerror = reject;
    img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg);
  });
}
const CRC = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
function crc32(bytes) { let c = 0xffffffff; for (const b of bytes) c = CRC[(c ^ b) & 0xff] ^ (c >>> 8); return (c ^ 0xffffffff) >>> 0; }
async function setPngDpi(blob, dpi) {
  const src = new Uint8Array(await blob.arrayBuffer());
  const ppm = Math.round(dpi / 0.0254);
  const chunk = new Uint8Array(21); const dv = new DataView(chunk.buffer);
  dv.setUint32(0, 9); chunk.set([0x70, 0x48, 0x59, 0x73], 4); dv.setUint32(8, ppm); dv.setUint32(12, ppm); chunk[16] = 1;
  dv.setUint32(17, crc32(chunk.subarray(4, 17)));
  const ihdrEnd = 8 + 25; // signature + IHDR chunk
  const out = new Uint8Array(src.length + 21); out.set(src.subarray(0, ihdrEnd)); out.set(chunk, ihdrEnd); out.set(src.subarray(ihdrEnd), ihdrEnd + 21);
  return new Blob([out], { type: 'image/png' });
}
