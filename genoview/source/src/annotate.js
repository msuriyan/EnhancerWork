// Annotation editor: put any figure (SVG) on a stage, add text boxes / arrows / boxes, export at any DPI.
import { svgToPngBlob } from './figureEditor.js';

const FONTS = ['Arial, Helvetica, sans-serif', 'Helvetica, Arial, sans-serif', '"Times New Roman", Times, serif', 'Georgia, serif', 'Calibri, Carlito, sans-serif', 'Verdana, sans-serif', '"Courier New", monospace', 'system-ui, sans-serif'];
const NS = 'http://www.w3.org/2000/svg';
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const store = {};   // per-key annotations, kept for the session
let lastTextStyle = null;
export function getAnnotationStore() { return JSON.parse(JSON.stringify(store, (k, v) => k === '_bb' ? undefined : v)); }
export function setAnnotationStore(obj) { for (const k of Object.keys(store)) delete store[k]; Object.assign(store, obj || {}); }

const TEXT_DEFAULT = { type: 'text', text: 'Label', font: 'Arial, Helvetica, sans-serif', size: 16, color: '#111111', bold: false, italic: false, underline: false, align: 'left', bg: '#ffffff', bgOpacity: 0, border: '#111111', borderWidth: 0, pad: 4, rotate: 0, lineHeight: 1.2 };
const ARROW_DEFAULT = { type: 'arrow', color: '#e34948', width: 2.5, head: 12, dashed: false, headStart: false };
const BOX_DEFAULT = { type: 'box', color: '#e34948', width: 2, fill: '#e34948', fillOpacity: 0, dashed: true, radius: 3 };

function parseSize(svg) {
  const d = new DOMParser().parseFromString(svg, 'image/svg+xml').documentElement;
  const vb = (d.getAttribute('viewBox') || '').split(/[\s,]+/).map(Number);
  const w = parseFloat(d.getAttribute('width')) || vb[2] || 800, h = parseFloat(d.getAttribute('height')) || vb[3] || 600;
  return { w, h };
}

/** opts: { svg, key, title, filename } */
export function openAnnotator(opts) {
  document.getElementById('annEditor')?.remove();
  const { w: BW, h: BH } = parseSize(opts.svg);
  const st = store[opts.key] ||= { items: [], margin: { t: 0, r: 0, b: 0, l: 0 }, bg: '#ffffff', transparent: false, dpi: 300, scale: 1 };
  st.ov ||= {};
  st.crop ||= { t: 0, r: 0, b: 0, l: 0 };   // overlay set name -> { show, labels, opacity }
  const overlays = opts.overlays || [];
  overlays.forEach(o => { st.ov[o.name] ||= { show: true, labels: true, lines: true }; });
  let sel = null;
  const wrap = document.createElement('div'); wrap.id = 'annEditor'; wrap.className = 'figmodal';
  wrap.innerHTML = `<div class="figbox">
    <div class="fighead"><b>${esc(opts.title || 'Annotate & export')}</b>
      <span class="muted" style="margin-left:8px">drag items to move · double-click text to edit · Delete removes · arrow keys nudge</span><span style="flex:1"></span>
      <button id="anText" class="primary">+ Text box</button><button id="anArrow">+ Arrow</button><button id="anBox">+ Box</button>
      <span class="sep" style="width:1px;height:20px;background:var(--line);margin:0 4px"></span><button id="anClose">Close</button></div>
    <div class="figbody"><div class="figctl" id="anCtl"></div>
      <div class="figprev"><div id="anStageWrap" style="flex:1;overflow:auto;padding:16px;background:repeating-conic-gradient(#f3f3f0 0% 25%,#fff 0% 50%) 50%/20px 20px"><div id="anStage" style="width:max-content;margin:auto"></div></div>
      <div class="figexp">
        <label class="status">Zoom <select id="anZoom">${[0.5, 0.75, 1, 1.5, 2].map(z => `<option value="${z}" ${z === 1 ? 'selected' : ''}>${z * 100}%</option>`).join('')}</select></label>
        <span class="sep" style="width:1px;height:20px;background:var(--line)"></span>
        <label class="status">Resolution <select id="anDpi">${[96, 150, 300, 600, 1200].map(d => `<option value="${d}" ${d === st.dpi ? 'selected' : ''}>${d} dpi</option>`).join('')}<option value="custom">custom…</option></select></label>
        <input type="number" id="anDpiCustom" min="50" max="2400" value="${st.dpi}" style="width:70px;display:none">
        <button id="anPng" class="primary">PNG</button><button id="anJpg">JPEG</button><button id="anSvg">SVG (vector)</button><button id="anPdf" title="Print dialog → Save as PDF">PDF</button>
        <span class="status" id="anInfo"></span></div></div></div></div>`;
  document.body.appendChild(wrap);
  const $ = q => wrap.querySelector(q);
  const stage = $('#anStage');

  // build the SVG stage
  let svg, itemsG, baseG;
  function buildStage() {
    const m = st.margin, c = st.crop, CW = Math.max(20, BW - c.l - c.r), CH = Math.max(20, BH - c.t - c.b);
    const W = CW + m.l + m.r, H = CH + m.t + m.b;
    const z = +$('#anZoom').value;
    stage.innerHTML = `<svg xmlns="${NS}" id="anSvgRoot" width="${W * z}" height="${H * z}" viewBox="0 0 ${W} ${H}" style="background:${st.transparent ? 'transparent' : st.bg};box-shadow:0 1px 6px rgba(0,0,0,.15);display:block;user-select:none">
      <defs><marker id="anHead" viewBox="0 0 10 10" refX="9" refY="5" markerUnits="userSpaceOnUse" orient="auto-start-reverse" markerWidth="12" markerHeight="12"><path d="M0,0L10,5L0,10z"/></marker></defs>
      ${st.transparent ? '' : `<rect class="an-bg" width="${W}" height="${H}" fill="${st.bg}"/>`}
      <clipPath id="anCrop"><rect x="${m.l}" y="${m.t}" width="${CW}" height="${CH}"/></clipPath>
      <g clip-path="url(#anCrop)"><g class="an-base" transform="translate(${m.l - c.l} ${m.t - c.t})"></g></g><g class="an-items"></g></svg>`;
    svg = stage.querySelector('svg'); baseG = svg.querySelector('.an-base'); itemsG = svg.querySelector('.an-items');
    // nest the original figure as an inner <svg> so its own coordinate system is kept
    const inner = new DOMParser().parseFromString(opts.svg, 'image/svg+xml').documentElement;
    inner.setAttribute('x', 0); inner.setAttribute('y', 0); inner.setAttribute('width', BW); inner.setAttribute('height', BH);
    inner.querySelectorAll('.hover, .tip').forEach(n => n.remove());
    baseG.appendChild(document.importNode(inner, true));
    drawOverlays();
    drawItems();
    $('#anInfo').textContent = `${Math.round(W)}×${Math.round(H)} px figure`;
    updateInfo();
  }
  const dpi = () => $('#anDpi').value === 'custom' ? Math.max(50, +$('#anDpiCustom').value || 300) : +$('#anDpi').value;
  function updateInfo() {
    const m = st.margin, c = st.crop, W = Math.max(20, BW - c.l - c.r) + m.l + m.r, H = Math.max(20, BH - c.t - c.b) + m.t + m.b, d = dpi();
    $('#anInfo').textContent = `${Math.round(W)}×${Math.round(H)} px figure → PNG ${Math.round(W * d / 96)}×${Math.round(H * d / 96)} px · ${(W / 96 * 25.4).toFixed(0)}×${(H / 96 * 25.4).toFixed(0)} mm at ${d} dpi`;
  }

  // ---------- highlighted regions from the browser (guides, cut sites, calls) ----------
  function solid(c) { const m = /rgba?\(([^)]+)\)/.exec(c); if (!m) return c; const p = m[1].split(',').map(x => x.trim()); return `rgb(${p[0]},${p[1]},${p[2]})`; }
  function alpha(c) { const m = /rgba\(([^)]+)\)/.exec(c); return m ? +m[1].split(',')[3] : 0.25; }
  function drawOverlays() {
    if (!overlays.length) return;
    const g = document.createElementNS(NS, 'g'); g.setAttribute('class', 'an-roi');
    let html = '';
    for (const o of overlays) {
      const v = st.ov[o.name]; if (!v.show) continue;
      const col = solid(o.color), a = Math.max(0.08, alpha(o.color));
      for (const r of o.rects) {
        if (r.w < 3) { // narrow region (cut site) → line
          if (v.lines) html += `<line x1="${r.x + r.w / 2}" x2="${r.x + r.w / 2}" y1="0" y2="${BH}" stroke="${col}" stroke-width="1.5" stroke-dasharray="5 3"/>`;
        } else html += `<rect x="${r.x}" y="0" width="${r.w}" height="${BH}" fill="${col}" fill-opacity="${a}"/><rect x="${r.x}" y="0" width="${r.w}" height="8" fill="${col}" fill-opacity="${Math.min(1, a * 3)}"/>`;
        if (v.labels && r.label) html += `<text x="${r.x + r.w / 2}" y="${r.w < 3 ? 32 : 20}" text-anchor="middle" font-family="Arial, sans-serif" font-size="11" font-weight="600" fill="${col}">${esc(r.label)}</text>`;
      }
    }
    g.innerHTML = html; baseG.appendChild(g);
  }

  // ---------- render items ----------
  function drawItems() {
    itemsG.innerHTML = '';
    st.items.forEach((it, i) => {
      const g = document.createElementNS(NS, 'g'); g.dataset.i = i; g.style.cursor = 'move';
      itemsG.appendChild(g);              // must be in the document before text is measured
      if (it.type === 'text') renderText(g, it);
      else if (it.type === 'arrow') renderArrow(g, it);
      else renderBox(g, it);
      if (sel === i) addHandles(g, it, i);
    });
  }
  function renderText(g, it) {
    const lines = String(it.text).split('\n');
    const t = document.createElementNS(NS, 'text');
    t.setAttribute('font-family', it.font); t.setAttribute('font-size', it.size); t.setAttribute('fill', it.color);
    if (it.bold) t.setAttribute('font-weight', '700'); if (it.italic) t.setAttribute('font-style', 'italic'); if (it.underline) t.setAttribute('text-decoration', 'underline');
    const anchor = { left: 'start', center: 'middle', right: 'end' }[it.align];
    t.setAttribute('text-anchor', anchor);
    lines.forEach((ln, k) => { const ts = document.createElementNS(NS, 'tspan'); ts.setAttribute('x', 0); ts.setAttribute('dy', k ? it.size * it.lineHeight : it.size * 0.8); ts.textContent = ln || ' '; t.appendChild(ts); });
    g.appendChild(t);
    const bb = t.getBBox(); const p = it.pad;
    if (it.bgOpacity > 0 || it.borderWidth > 0) {
      const r = document.createElementNS(NS, 'rect');
      r.setAttribute('x', bb.x - p); r.setAttribute('y', bb.y - p); r.setAttribute('width', bb.width + 2 * p); r.setAttribute('height', bb.height + 2 * p); r.setAttribute('rx', 2);
      r.setAttribute('fill', it.bg); r.setAttribute('fill-opacity', it.bgOpacity);
      if (it.borderWidth > 0) { r.setAttribute('stroke', it.border); r.setAttribute('stroke-width', it.borderWidth); }
      g.insertBefore(r, t);
    }
    g.setAttribute('transform', `translate(${it.x} ${it.y}) rotate(${it.rotate || 0})`);
    it._bb = { x: bb.x - p, y: bb.y - p, w: bb.width + 2 * p, h: bb.height + 2 * p };
  }
  function renderArrow(g, it) {
    const l = document.createElementNS(NS, 'line');
    Object.entries({ x1: it.x1, y1: it.y1, x2: it.x2, y2: it.y2, stroke: it.color, 'stroke-width': it.width, 'stroke-linecap': 'round' }).forEach(([k, v]) => l.setAttribute(k, v));
    if (it.dashed) l.setAttribute('stroke-dasharray', `${it.width * 3} ${it.width * 2}`);
    const id = 'ah' + Math.random().toString(36).slice(2, 8);
    const mk = document.createElementNS(NS, 'marker');
    Object.entries({ id, viewBox: '0 0 10 10', refX: 8, refY: 5, markerUnits: 'userSpaceOnUse', orient: 'auto-start-reverse', markerWidth: it.head, markerHeight: it.head }).forEach(([k, v]) => mk.setAttribute(k, v));
    mk.innerHTML = `<path d="M0,0L10,5L0,10z" fill="${it.color}"/>`;
    g.appendChild(mk);
    l.setAttribute('marker-end', `url(#${id})`); if (it.headStart) l.setAttribute('marker-start', `url(#${id})`);
    const hit = l.cloneNode(); hit.setAttribute('stroke', 'transparent'); hit.setAttribute('stroke-width', Math.max(12, it.width + 8)); hit.removeAttribute('marker-end'); hit.removeAttribute('marker-start'); hit.removeAttribute('stroke-dasharray'); hit.classList.add('an-ui');
    g.appendChild(l); g.appendChild(hit);
  }
  function renderBox(g, it) {
    const r = document.createElementNS(NS, 'rect');
    Object.entries({ x: Math.min(it.x1, it.x2), y: Math.min(it.y1, it.y2), width: Math.abs(it.x2 - it.x1), height: Math.abs(it.y2 - it.y1), rx: it.radius, fill: it.fill, 'fill-opacity': it.fillOpacity, stroke: it.color, 'stroke-width': it.width }).forEach(([k, v]) => r.setAttribute(k, v));
    if (it.dashed) r.setAttribute('stroke-dasharray', `${it.width * 3} ${it.width * 2}`);
    if (it.fillOpacity === 0) r.setAttribute('pointer-events', 'stroke');
    g.appendChild(r);
  }
  function addHandles(g, it, i) {
    const mk = (x, y, which, abs = true) => { const c = document.createElementNS(NS, 'circle'); c.setAttribute('cx', x); c.setAttribute('cy', y); c.setAttribute('r', 5); c.setAttribute('fill', '#fff'); c.setAttribute('stroke', '#2a78d6'); c.setAttribute('stroke-width', 2); c.classList.add('an-ui', 'an-handle'); c.dataset.h = which; c.style.cursor = 'crosshair'; (abs ? itemsG : g).appendChild(c); };
    if (it.type === 'text') {
      const b = it._bb; const r = document.createElementNS(NS, 'rect');
      Object.entries({ x: b.x - 2, y: b.y - 2, width: b.w + 4, height: b.h + 4, fill: 'none', stroke: '#2a78d6', 'stroke-dasharray': '4 3', 'stroke-width': 1 }).forEach(([k, v]) => r.setAttribute(k, v));
      r.classList.add('an-ui'); g.appendChild(r);
    } else { mk(it.x1, it.y1, '1'); mk(it.x2, it.y2, '2'); }
  }

  // ---------- interaction ----------
  const toSvg = e => { const pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY; return pt.matrixTransform(svg.getScreenCTM().inverse()); };
  let drag = null;
  stage.addEventListener('pointerdown', e => {
    const h = e.target.closest('.an-handle');
    const g = e.target.closest('.an-items > g');
    if (h) { drag = { mode: 'h' + h.dataset.h, i: sel, p0: toSvg(e) }; }
    else if (g) { sel = +g.dataset.i; drag = { mode: 'move', i: sel, p0: toSvg(e), it0: { ...st.items[sel] } }; drawItems(); buildProps(); }
    else { sel = null; drawItems(); buildProps(); return; }
    stage.setPointerCapture(e.pointerId); e.preventDefault();
  });
  stage.addEventListener('pointermove', e => {
    if (!drag) return;
    const p = toSvg(e), it = st.items[drag.i];
    if (drag.mode === 'move') {
      const dx = p.x - drag.p0.x, dy = p.y - drag.p0.y, o = drag.it0;
      if (it.type === 'text') { it.x = o.x + dx; it.y = o.y + dy; } else { it.x1 = o.x1 + dx; it.y1 = o.y1 + dy; it.x2 = o.x2 + dx; it.y2 = o.y2 + dy; }
    } else if (drag.mode === 'h1') { it.x1 = p.x; it.y1 = p.y; } else if (drag.mode === 'h2') { it.x2 = p.x; it.y2 = p.y; }
    drawItems();
  });
  stage.addEventListener('pointerup', () => { if (drag) { drag = null; buildProps(); } });
  stage.addEventListener('dblclick', e => { const g = e.target.closest('.an-items > g'); if (g && st.items[+g.dataset.i].type === 'text') { sel = +g.dataset.i; buildProps(); const ta = $('#anCtl textarea'); ta?.focus(); ta?.select(); } });
  const onKey = e => {
    if (!document.body.contains(wrap)) return document.removeEventListener('keydown', onKey);
    if (sel == null || /INPUT|TEXTAREA|SELECT/.test(document.activeElement?.tagName)) return;
    const it = st.items[sel]; const d = e.shiftKey ? 10 : 1;
    const mv = (dx, dy) => { if (it.type === 'text') { it.x += dx; it.y += dy; } else { it.x1 += dx; it.x2 += dx; it.y1 += dy; it.y2 += dy; } drawItems(); };
    if (e.key === 'Delete' || e.key === 'Backspace') { st.items.splice(sel, 1); sel = null; drawItems(); buildProps(); e.preventDefault(); }
    else if (e.key === 'ArrowLeft') { mv(-d, 0); e.preventDefault(); } else if (e.key === 'ArrowRight') { mv(d, 0); e.preventDefault(); }
    else if (e.key === 'ArrowUp') { mv(0, -d); e.preventDefault(); } else if (e.key === 'ArrowDown') { mv(0, d); e.preventDefault(); }
  };
  document.addEventListener('keydown', onKey);

  // ---------- add items ----------
  const centre = () => { const m = st.margin, c = st.crop; return { x: m.l + (BW - c.l - c.r) / 2, y: m.t + (BH - c.t - c.b) / 2 }; };
  $('#anText').onclick = () => { const c = centre(); st.items.push({ ...TEXT_DEFAULT, ...(lastTextStyle || {}), text: 'Label', x: c.x - 30, y: c.y }); sel = st.items.length - 1; drawItems(); buildProps(); const ta = $('#anCtl textarea'); ta?.focus(); ta?.select(); };
  $('#anArrow').onclick = () => { const c = centre(); st.items.push({ ...ARROW_DEFAULT, x1: c.x - 80, y1: c.y - 50, x2: c.x, y2: c.y }); sel = st.items.length - 1; drawItems(); buildProps(); };
  $('#anBox').onclick = () => { const c = centre(); st.items.push({ ...BOX_DEFAULT, x1: c.x - 60, y1: c.y - 30, x2: c.x + 60, y2: c.y + 30 }); sel = st.items.length - 1; drawItems(); buildProps(); };

  // ---------- properties panel ----------
  function buildProps() {
    const ctl = $('#anCtl'); const it = sel != null ? st.items[sel] : null;
    const m = st.margin;
    let html = `<fieldset><legend>Canvas</legend>
      <div class="fe-grid">${['t', 'r', 'b', 'l'].map(k => `<label class="fe"><span>Margin ${({ t: 'top', r: 'right', b: 'bottom', l: 'left' })[k]}</span><input type="number" data-m="${k}" value="${m[k]}" min="0" max="2000" step="10"></label>`).join('')}</div>
      <label class="fe"><span>Background</span><input type="color" id="anBg" value="${st.bg}"></label>
      <label class="fe"><span>Transparent (PNG/SVG)</span><input type="checkbox" id="anTr" ${st.transparent ? 'checked' : ''}></label>
      <div class="hint">Margins add empty space around the figure for labels outside it.</div>
      <div class="fe-grid" style="margin-top:6px">${['t', 'r', 'b', 'l'].map(k => `<label class="fe"><span>Crop ${({ t: 'top', r: 'right', b: 'bottom', l: 'left' })[k]}</span><input type="number" data-c="${k}" value="${st.crop[k]}" min="0" max="5000" step="5"></label>`).join('')}</div>
      <div class="rulebar" style="margin:4px 0"><button id="anAutoCrop" title="Cut the figure to its drawn content (removes empty space at the edges)">Auto-crop to content</button><button id="anNoCrop">No crop</button></div></fieldset>
      ${overlays.length ? `<fieldset><legend>Highlighted regions (as on screen)</legend>${overlays.map((o, k) => `<div style="margin:4px 0"><b style="font-size:12px"><span class="swatch" style="background:${solid(o.color)}"></span> ${esc(o.name)}</b> <span class="muted">(${o.rects.length})</span><br>
        <label class="status"><input type="checkbox" data-ov="${k}" data-f="show" ${st.ov[o.name].show ? 'checked' : ''}> show</label>
        <label class="status"><input type="checkbox" data-ov="${k}" data-f="labels" ${st.ov[o.name].labels ? 'checked' : ''}> labels</label>
        <label class="status"><input type="checkbox" data-ov="${k}" data-f="lines" ${st.ov[o.name].lines ? 'checked' : ''}> cut lines</label></div>`).join('')}</fieldset>` : ''}
      <fieldset><legend>Items (${st.items.length})</legend>${st.items.length ? st.items.map((x, i) => `<div class="an-li ${i === sel ? 'on' : ''}" data-li="${i}">${x.type === 'text' ? '𝐓 ' + esc(String(x.text).split('\n')[0].slice(0, 28)) : x.type === 'arrow' ? '➔ Arrow' : '▭ Box'}</div>`).join('') : '<div class="hint">Add text boxes, arrows or boxes with the buttons above.</div>'}</fieldset>`;
    if (it) {
      const opt = (k, list) => `<select data-p="${k}">${list.map(([v, l]) => `<option value="${esc(v)}" ${String(it[k]) === String(v) ? 'selected' : ''}>${esc(l)}</option>`).join('')}</select>`;
      const num = (k, a, b, s = 1) => `<input type="number" data-p="${k}" value="${it[k]}" min="${a}" max="${b}" step="${s}">`;
      const col = k => `<input type="color" data-p="${k}" value="${it[k]}">`;
      const chk = (k, l) => `<label class="status" style="margin-right:8px"><input type="checkbox" data-p="${k}" ${it[k] ? 'checked' : ''}> ${l}</label>`;
      const F = (l, h) => `<label class="fe"><span>${l}</span>${h}</label>`;
      if (it.type === 'text') html += `<fieldset><legend>Text box</legend>
        <textarea data-p="text" rows="3" style="width:100%;font:13px var(--mono)">${esc(it.text)}</textarea>
        ${F('Font', `<select data-p="font">${[...new Set([it.font, ...FONTS])].map(f => `<option value='${esc(f)}' ${f === it.font ? 'selected' : ''}>${esc(f.split(',')[0].replace(/"/g, ''))}</option>`).join('')}</select>`)}
        ${F('…or installed font', `<input type="text" id="anFontFree" placeholder="e.g. Myriad Pro">`)}
        <div class="fe-grid">${F('Size', num('size', 4, 300))}${F('Colour', col('color'))}</div>
        <div style="margin:4px 0">${chk('bold', '<b>Bold</b>')}${chk('italic', '<i>Italic</i>')}${chk('underline', '<u>Underline</u>')}</div>
        <div class="fe-grid">${F('Align', opt('align', [['left', 'Left'], ['center', 'Centre'], ['right', 'Right']]))}${F('Line height', num('lineHeight', 0.8, 3, 0.1))}</div>
        <div class="fe-grid">${F('Background', col('bg'))}${F('Opacity', num('bgOpacity', 0, 1, 0.1))}${F('Border', col('border'))}${F('Border w', num('borderWidth', 0, 20, 0.5))}${F('Padding', num('pad', 0, 50))}${F('Rotate °', num('rotate', -180, 180, 5))}</div>
        <button id="anStyleAll" title="Apply this font/size/colour style to all text boxes">Apply style to all text boxes</button></fieldset>`;
      else html += `<fieldset><legend>${it.type === 'arrow' ? 'Arrow' : 'Box'}</legend>
        <div class="fe-grid">${F('Colour', col('color'))}${F('Line width', num('width', 0.5, 30, 0.5))}
        ${it.type === 'arrow' ? F('Head size', num('head', 4, 60)) : F('Fill', col('fill')) + F('Fill opacity', num('fillOpacity', 0, 1, 0.1)) + F('Corner', num('radius', 0, 50))}</div>
        <div style="margin:4px 0">${chk('dashed', 'Dashed')}${it.type === 'arrow' ? chk('headStart', 'Head at both ends') : ''}</div>
        <div class="hint">Drag the round handles to set the ${it.type === 'arrow' ? 'start and end' : 'corners'}.</div></fieldset>`;
      html += `<div class="rulebar"><button id="anDup">Duplicate</button><button id="anFront">Bring to front</button><button id="anDel" style="color:var(--del)">Delete</button></div>`;
    }
    ctl.innerHTML = html;
    ctl.querySelectorAll('[data-m]').forEach(el => el.oninput = () => { const k = el.dataset.m, nv = Math.max(0, +el.value || 0), dv = nv - st.margin[k]; st.margin[k] = nv; if (k === 'l' || k === 't') st.items.forEach(x => { if (x.type === 'text') { if (k === 'l') x.x += dv; else x.y += dv; } else { if (k === 'l') { x.x1 += dv; x.x2 += dv; } else { x.y1 += dv; x.y2 += dv; } } }); buildStage(); });
    ctl.querySelector('#anBg').oninput = e => { st.bg = e.target.value; buildStage(); };
    ctl.querySelectorAll('[data-c]').forEach(el => el.onchange = () => { const k = el.dataset.c, nv = Math.max(0, +el.value || 0), dv = nv - st.crop[k]; st.crop[k] = nv; if (k === 'l' || k === 't') shiftItems(k, -dv); buildStage(); });
    ctl.querySelector('#anAutoCrop').onclick = () => { autoCrop(); buildStage(); buildProps(); };
    ctl.querySelector('#anNoCrop').onclick = () => { shiftItems('l', st.crop.l); shiftItems('t', st.crop.t); st.crop = { t: 0, r: 0, b: 0, l: 0 }; buildStage(); buildProps(); };
    ctl.querySelectorAll('[data-ov]').forEach(el => el.onchange = () => { st.ov[overlays[+el.dataset.ov].name][el.dataset.f] = el.checked; buildStage(); });
    ctl.querySelector('#anTr').onchange = e => { st.transparent = e.target.checked; buildStage(); };
    ctl.querySelectorAll('[data-li]').forEach(el => el.onclick = () => { sel = +el.dataset.li; drawItems(); buildProps(); });
    if (!it) return;
    ctl.querySelectorAll('[data-p]').forEach(el => el.addEventListener(el.tagName === 'TEXTAREA' || el.type === 'text' || el.type === 'number' || el.type === 'color' ? 'input' : 'change', () => {
      const k = el.dataset.p; it[k] = el.type === 'checkbox' ? el.checked : el.type === 'number' ? +el.value : el.value;
      if (it.type === 'text') lastTextStyle = pickStyle(it);
      drawItems();
      if (k === 'text') { const li = ctl.querySelector(`[data-li="${sel}"]`); if (li) li.textContent = '𝐓 ' + String(it.text).split('\n')[0].slice(0, 28); }
    }));
    const ff = ctl.querySelector('#anFontFree'); if (ff) ff.onchange = () => { if (ff.value.trim()) { it.font = `"${ff.value.trim()}", Arial, sans-serif`; lastTextStyle = pickStyle(it); drawItems(); buildProps(); } };
    const sa = ctl.querySelector('#anStyleAll'); if (sa) sa.onclick = () => { const s = pickStyle(it); st.items.forEach(x => { if (x.type === 'text') Object.assign(x, s); }); drawItems(); };
    ctl.querySelector('#anDup').onclick = () => { const c = JSON.parse(JSON.stringify(it)); if (c.type === 'text') { c.x += 15; c.y += 15; } else { c.x1 += 15; c.x2 += 15; c.y1 += 15; c.y2 += 15; } st.items.push(c); sel = st.items.length - 1; drawItems(); buildProps(); };
    ctl.querySelector('#anFront').onclick = () => { st.items.push(st.items.splice(sel, 1)[0]); sel = st.items.length - 1; drawItems(); buildProps(); };
    ctl.querySelector('#anDel').onclick = () => { st.items.splice(sel, 1); sel = null; drawItems(); buildProps(); };
  }
  function shiftItems(k, dv) { st.items.forEach(x => { if (x.type === 'text') { if (k === 'l') x.x += dv; else x.y += dv; } else { if (k === 'l') { x.x1 += dv; x.x2 += dv; } else { x.y1 += dv; x.y2 += dv; } } }); }
  // find the extent of non-blank drawing in the base figure (ignores the white page background)
  function autoCrop() {
    const inner = baseG.querySelector('svg'); if (!inner) return;
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    const white = c => { c = (c || '').toLowerCase().replace(/\s+/g, ''); return c === 'white' || c === '#fff' || c === '#ffffff' || c === 'rgb(255,255,255)' || c === 'rgba(255,255,255,1)'; };
    const none = c => { c = (c || '').toLowerCase(); return c === 'none' || c === 'transparent' || /rgba\([^)]*,\s*0\)$/.test(c); };
    inner.querySelectorAll('path,rect,line,polyline,polygon,text,circle,image,ellipse').forEach(el => {
      if (el.closest('clipPath,defs,mask,marker,pattern,symbol')) return;
      const f = el.getAttribute('fill') ?? el.style.fill, sk = el.getAttribute('stroke') ?? el.style.stroke;
      const fillVisible = el.tagName === 'text' || el.tagName === 'image' || !(white(f) || none(f));
      const strokeVisible = sk && !none(sk) && !white(sk);
      if (!fillVisible && !strokeVisible) return;
      if (el.tagName === 'line' && !strokeVisible) return;
      let b; try { b = el.getBBox(); } catch { return; }
      if (!b || (b.width === 0 && b.height === 0)) return;
      const ctm = el.getCTM(), ictm = inner.getCTM().inverse();
      const pts = [[b.x, b.y], [b.x + b.width, b.y + b.height]].map(([x, y]) => { const p = svg.createSVGPoint(); p.x = x; p.y = y; return p.matrixTransform(ctm).matrixTransform(ictm); });
      for (const p of pts) { x0 = Math.min(x0, p.x); y0 = Math.min(y0, p.y); x1 = Math.max(x1, p.x); y1 = Math.max(y1, p.y); }
    });
    if (!isFinite(x0)) return;
    const pad = 6;
    const nc = { l: Math.max(0, Math.floor(x0 - pad)), t: Math.max(0, Math.floor(y0 - pad)), r: Math.max(0, Math.floor(BW - x1 - pad)), b: Math.max(0, Math.floor(BH - y1 - pad)) };
    shiftItems('l', st.crop.l - nc.l); shiftItems('t', st.crop.t - nc.t);
    st.crop = nc;
  }
  const pickStyle = it => ({ font: it.font, size: it.size, color: it.color, bold: it.bold, italic: it.italic, underline: it.underline, bg: it.bg, bgOpacity: it.bgOpacity, border: it.border, borderWidth: it.borderWidth, pad: it.pad });

  // ---------- export ----------
  function exportSvg() {
    const prev = sel; sel = null; drawItems();
    const c = svg.cloneNode(true);
    c.querySelectorAll('.an-ui').forEach(n => n.remove());
    c.querySelectorAll('[style]').forEach(n => { if (/cursor/.test(n.getAttribute('style'))) n.removeAttribute('style'); });
    c.querySelectorAll('[data-i]').forEach(n => n.removeAttribute('data-i'));
    c.removeAttribute('style'); c.removeAttribute('id');
    const vb = c.getAttribute('viewBox').split(' ').map(Number);
    c.setAttribute('width', vb[2]); c.setAttribute('height', vb[3]);
    const out = new XMLSerializer().serializeToString(c);
    sel = prev; drawItems();
    return { svg: out, w: vb[2], h: vb[3] };
  }
  const fname = ext => `${opts.filename || 'figure'}.${ext}`;
  const dl = (name, blob) => { const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500); };
  const busy = (btn, fn) => async () => { const t = btn.textContent; btn.disabled = true; btn.textContent = 'Rendering…'; try { await fn(); } catch (e) { alert('Export failed: ' + (e.message || e)); } finally { btn.disabled = false; btn.textContent = t; } };
  $('#anSvg').onclick = () => dl(fname('svg'), new Blob(['<?xml version="1.0" encoding="UTF-8"?>\n' + exportSvg().svg], { type: 'image/svg+xml' }));
  $('#anPng').onclick = busy($('#anPng'), async () => { const e = exportSvg(); dl(fname('png'), await svgToPngBlob(e.svg, e.w, e.h, dpi(), st.transparent)); });
  $('#anJpg').onclick = busy($('#anJpg'), async () => {
    const e = exportSvg(); const png = await svgToPngBlob(e.svg, e.w, e.h, dpi(), false);
    const img = await createImageBitmap(png); const cv = document.createElement('canvas'); cv.width = img.width; cv.height = img.height;
    const ctx = cv.getContext('2d'); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, cv.width, cv.height); ctx.drawImage(img, 0, 0);
    dl(fname('jpg'), await new Promise(r => cv.toBlob(r, 'image/jpeg', 0.95)));
  });
  $('#anPdf').onclick = () => {
    const e = exportSvg(); const w = window.open('', '_blank'); if (!w) return alert('Allow pop-ups to export PDF.');
    const mmW = e.w / 96 * 25.4, mmH = e.h / 96 * 25.4;
    w.document.write(`<!DOCTYPE html><html><head><title>${esc(fname('pdf'))}</title><style>@page{size:${mmW}mm ${mmH}mm;margin:0}html,body{margin:0}svg{display:block;width:${mmW}mm;height:${mmH}mm}</style></head><body>${e.svg}<script>setTimeout(()=>print(),400)<\/script></body></html>`);
    w.document.close();
  };
  $('#anDpi').onchange = () => { $('#anDpiCustom').style.display = $('#anDpi').value === 'custom' ? '' : 'none'; st.dpi = dpi(); updateInfo(); };
  $('#anDpiCustom').oninput = () => { st.dpi = dpi(); updateInfo(); };
  $('#anZoom').onchange = buildStage;
  $('#anClose').onclick = () => wrap.remove();
  buildStage(); buildProps();
}
