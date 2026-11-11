// GenoView sessions: save everything except the data itself to a small JSON file; re-open it later by
// pointing at the data folder/files once (browsers never reveal file paths, so files are matched by name + size).
// Forward/backward compatible: every field is optional on load, unknown fields are ignored, `version` is migrated.

export const SESSION_FORMAT = 'genoview-session';
export const SESSION_VERSION = 1;

// UI controls whose values are part of a session (missing ids are skipped, so old/new builds both work)
const CONTROL_IDS = ['pLib', 'pRegion', 'pCustom', 'pMapq', 'pMinSV', 'pSupport', 'pAF', 'pK', 'pDup',
  'fSize', 'fROI', 'fShowHidden', 'gMode', 'gPooled', 'gCalls',
  'cTarget', 'cNuc', 'cReads', 'cKinds', 'cMode', 'cWin', 'cROI', 'cFlank', 'cMapq', 'cMinIndel', 'cHp', 'cMerge', 'cIgnore'];
const TOGGLE_IDS = ['tgCenter', 'tgCursor', 'tgLabels', 'tgNames', 'tgArrows', 'tgIdeo', 'tgIdeoBox'];

export function initSession(ctx) {
  const { state, $, $$, toast, download, esc, fmt, fmeta } = ctx;

  // ---------------- save ----------------
  async function fileText(f) { try { return await f.text(); } catch { return ''; } }
  async function collect() {
    const br = state.browser;
    const idx = new Map(state.tracks.map((t, i) => [t, i]));
    const fileToIdx = new Map(state.tracks.filter(t => t.file && !t.generated).map(t => [t.file, idx.get(t)]));
    const tracks = [];
    for (const r of state.tracks) {
      if (r.isOverview) continue;                       // rebuilt automatically
      const t = r.igv || {};
      const base = { id: r.id, name: r.name, format: r.format, color: r.color, height: Math.round(t.height || r.height || 0) || undefined, mode: r.mode,
        hideNames: !!(t.gvHideNames ?? r.hideNames), hideArrows: !!(t.gvHideArrows ?? r.hideArrows), compact: r.compact || undefined };
      if (r.format === 'bam') Object.assign(base, { showCoverage: t.showCoverage !== false, showAlignments: t.showAlignments !== false, coverageHeight: t.coverageTrackHeight, sortMode: r.sortMode || '', displayMode: t.alignmentTrack?.displayMode });
      if (r.sample) base.sample = { use: r.sample.use, group: r.sample.group, groupEdited: !!r.sample.groupEdited, isControl: !!r.sample.isControl };
      if (!r.generated && r.file) { base.file = fmeta(r.file); if (r.index) base.index = fmeta(r.index); }
      else {
        // generated tracks: small, so their content travels inside the session
        const cfg = r.cfg || {};
        const g = { type: cfg.type || null, format: cfg.format || r.format, cfg: {} };
        for (const k of ['displayMode', 'itemRgb', 'autoscale', 'height', 'color']) if (cfg[k] !== undefined) g.cfg[k] = cfg[k];
        if (cfg.type === 'merged') {
          g.members = [];
          for (const m of cfg.tracks || []) {
            if (fileToIdx.has(m.url)) g.members.push({ ref: fileToIdx.get(m.url), name: m.name, color: m.color, format: m.format });
            else if (m.url instanceof Blob) g.members.push({ text: await fileText(m.url), fileName: m.url.name, name: m.name, color: m.color, format: m.format });
          }
        } else if (cfg.url instanceof Blob) { g.text = await fileText(cfg.url); g.fileName = cfg.url.name; }
        else continue;
        base.generated = g;
      }
      tracks.push(base);
    }
    const controls = {}; for (const id of CONTROL_IDS) { const el = document.getElementById(id); if (el) controls[id] = el.type === 'checkbox' ? el.checked : el.value; }
    const toggles = {}; for (const id of TOGGLE_IDS) { const el = document.getElementById(id); if (el) toggles[id] = el.classList.contains('on'); }
    const figStyles = {}; try { for (let i = 0; i < localStorage.length; i++) { const k = localStorage.key(i); if (k.startsWith('genoview.fig.')) figStyles[k] = localStorage.getItem(k); } } catch {}
    return {
      format: SESSION_FORMAT, version: SESSION_VERSION, savedAt: new Date().toISOString(), app: 'GenoView',
      genome: state.genomeSrc || null,
      locus: br ? (() => { const c = br.currentLoci?.(); return Array.isArray(c) ? c.join(' ') : (c || null); })() : null,
      tracks, controls, toggles,
      hideRules: { rules: state.hideRules || [], analyses: !!state.hideRulesInAnalyses },
      sv: { ran: !!state.result, region: state.result ? `${state.result.chr}:${state.result.start + 1}-${state.result.end}` : null, hidden: [...(state.hidden || [])], typeFilter: [...(state.typeFilter || [])], maskIntervals: state.maskIntervals || [] },
      crispr: state.crisprSession?.get() || null,
      figStyles, annotations: ctx.getAnnotationStore?.() || {},
      dockHeight: document.getElementById('work')?.style.getPropertyValue('--dockH') || '',
      fullPage: !!document.getElementById('work')?.classList.contains('full'),
    };
  }

  async function saveSession() {
    if (!state.browser) return toast('Nothing to save yet — load a genome first.', true);
    let sess; try { sess = await collect(); } catch (e) { console.error(e); return toast('Could not save the session: ' + e.message, true, 9000); }
    const def = (state.lastSessionName || `genoview_session_${new Date().toISOString().slice(0, 10)}`).replace(/\.genoview\.json$/i, '');
    const name = prompt('Session name (the file is saved to your Downloads folder):', def);
    if (name === null) return;
    const fname = (name.trim() || def).replace(/[\\/:*?"<>|]+/g, '_') + '.genoview.json';
    const text = JSON.stringify(sess, null, 1);
    download(fname, text, 'application/json');
    state.lastSessionName = fname;
    remember(fname, text);
    toast(`Session saved as ${fname} — keep it with your data folder.`, false, 6000);
  }
  function remember(name, text) {
    try {
      const list = JSON.parse(localStorage.getItem('genoview.sessions') || '[]').filter(x => x.name !== name);
      list.unshift({ name, savedAt: new Date().toISOString(), text });
      while (list.length && JSON.stringify(list).length > 3e6) list.pop();
      localStorage.setItem('genoview.sessions', JSON.stringify(list.slice(0, 12)));
    } catch {}
  }

  // ---------------- open ----------------
  function migrate(s) {
    if (!s || typeof s !== 'object') throw new Error('Not a GenoView session file.');
    if (s.format && s.format !== SESSION_FORMAT) throw new Error('This JSON is not a GenoView session.');
    const v = +s.version || 1;
    if (v > SESSION_VERSION) toast(`This session was saved by a newer GenoView (format v${v}); loading what this version understands.`, false, 8000);
    // future migrations go here: if (v < 2) { … }
    return { tracks: [], controls: {}, toggles: {}, hideRules: {}, sv: {}, figStyles: {}, annotations: {}, ...s };
  }
  function requiredFiles(s) {
    const need = [];
    if (s.genome?.type === 'files' || s.genome?.type === 'bam-header') for (const f of s.genome.files || []) need.push({ ...f, role: 'genome' });
    for (const t of s.tracks || []) { if (t.file) need.push({ ...t.file, role: t.name }); if (t.index) need.push({ ...t.index, role: (t.name || '') + ' index' }); }
    const seen = new Set(); return need.filter(n => !seen.has(n.name) && seen.add(n.name));
  }

  function openSessionDialog(preset) {
    document.getElementById('sessModal')?.remove();
    const m = document.createElement('div'); m.id = 'sessModal'; m.className = 'figmodal';
    let recent = []; try { recent = JSON.parse(localStorage.getItem('genoview.sessions') || '[]'); } catch {}
    m.innerHTML = `<div class="figbox" style="width:min(820px,94vw);height:auto;max-height:90vh">
      <div class="fighead"><b>Sessions</b><span style="flex:1"></span><button id="ssClose">Close</button></div>
      <div style="padding:12px 16px;overflow:auto" id="ssBody">
        <div class="rulebar"><button id="ssSave" class="primary">Save current session…</button><button id="ssOpen">Open session file…</button><input type="file" id="ssFile" accept=".json" hidden></div>
        <div class="hint">A session stores tracks, order, colours, heights, display options, position, hide rules, CRISPR / SV settings and masks, figure styles and annotations — not the data files. To reopen, choose the session and then your data folder (or the files); they are matched by name and size. Sessions from older or newer GenoView versions open too.</div>
        ${recent.length ? `<h4 style="margin:14px 0 6px">Recent sessions (this browser)</h4><table class="rules"><tbody>${recent.map((r, i) => `<tr><td><b>${esc(r.name)}</b></td><td class="muted">${new Date(r.savedAt).toLocaleString()}</td><td style="text-align:right"><button data-rec="${i}">Open</button> <button class="ghost x" data-recdel="${i}" title="Forget">×</button></td></tr>`).join('')}</tbody></table>` : ''}
        <div id="ssStep2"></div>
      </div></div>`;
    document.body.appendChild(m);
    const q = sel => m.querySelector(sel);
    q('#ssClose').onclick = () => m.remove();
    m.addEventListener('mousedown', e => { if (e.target === m) m.remove(); });
    q('#ssSave').onclick = () => { m.remove(); saveSession(); };
    q('#ssOpen').onclick = () => q('#ssFile').click();
    q('#ssFile').onchange = async e => { const f = e.target.files[0]; if (!f) return; try { step2(migrate(JSON.parse(await f.text())), f.name); } catch (err) { toast(err.message, true); } };
    m.querySelectorAll('[data-rec]').forEach(b => b.onclick = () => { try { step2(migrate(JSON.parse(recent[+b.dataset.rec].text)), recent[+b.dataset.rec].name); } catch (err) { toast(err.message, true); } });
    m.querySelectorAll('[data-recdel]').forEach(b => b.onclick = () => { recent.splice(+b.dataset.recdel, 1); try { localStorage.setItem('genoview.sessions', JSON.stringify(recent)); } catch {} b.closest('tr').remove(); });
    if (preset) step2(preset.sess, preset.name);

    function step2(sess, name) {
      const need = requiredFiles(sess);
      const pool = new Map();   // name -> File
      const box = q('#ssStep2');
      const draw = () => {
        const found = need.filter(n => pool.has(n.name));
        box.innerHTML = `<h4 style="margin:16px 0 6px">${esc(name)} <span class="muted" style="font-weight:400">· saved ${sess.savedAt ? new Date(sess.savedAt).toLocaleString() : '?'} · ${sess.tracks.length} tracks${sess.genome?.type === 'ucsc' ? ' · genome ' + esc(sess.genome.db) + ' (online)' : ''}</span></h4>
          ${need.length ? `<div class="rulebar"><button id="ssFolder" class="primary">Choose data folder…</button><button id="ssFiles">Choose files…</button><input type="file" id="ssDir" webkitdirectory hidden><input type="file" id="ssMulti" multiple hidden><span class="status">or drop the files here · ${found.length}/${need.length} found</span></div>
          <table class="rules"><thead><tr><th>File</th><th>Used for</th><th>Status</th></tr></thead><tbody>${need.map(n => { const f = pool.get(n.name); const ok = f && (!n.size || f.size === n.size);
            return `<tr><td class="loc">${esc(n.name)}</td><td class="muted">${esc(n.role)}</td><td>${f ? (ok ? '<span style="color:var(--ins)">✓ found</span>' : '<span style="color:var(--cnv)">✓ found (size differs)</span>') : '<span style="color:var(--del)">missing</span>'}</td></tr>`; }).join('')}</tbody></table>` : '<div class="hint">This session needs no local files.</div>'}
          <div class="rulebar" style="margin-top:10px"><span style="flex:1"></span><label class="status"><input type="checkbox" id="ssRerun" checked> re-run analyses (Variants &amp; SVs, CRISPR)</label>
          <button id="ssGo" class="primary" ${need.length && !found.length && sess.genome?.type !== 'ucsc' ? 'disabled' : ''}>${found.length < need.length ? 'Load with what was found' : 'Load session'}</button></div>`;
        const add = files => { for (const f of files) { const want = need.find(n => n.name === f.name); if (!want) continue; const cur = pool.get(f.name); if (!cur || (want.size && f.size === want.size && cur.size !== want.size)) pool.set(f.name, f); } draw(); };
        if (q('#ssFolder')) { q('#ssFolder').onclick = () => q('#ssDir').click(); q('#ssFiles').onclick = () => q('#ssMulti').click(); q('#ssDir').onchange = e => add(e.target.files); q('#ssMulti').onchange = e => add(e.target.files); }
        box.ondragover = e => e.preventDefault();
        box.ondrop = e => { e.preventDefault(); e.stopPropagation(); add(e.dataTransfer.files); };
        q('#ssGo').onclick = async () => { const rerun = q('#ssRerun').checked; m.remove(); state.lastSessionName = name; await applySession(sess, pool, rerun); };
      };
      draw();
    }
  }

  // ---------------- apply ----------------
  async function applySession(s, pool, rerun = true) {
    const busy = msg => toast(msg, false, 60000);
    try {
      busy('Loading session…');
      // genome
      const g = s.genome;
      if (g?.type === 'ucsc') {
        const base = `https://hgdownload.soe.ucsc.edu/goldenPath/${g.db}/bigZips/`;
        await ctx.setGenome({ id: g.db, name: g.db, twoBitURL: base + g.db + '.2bit', chromSizesURL: base + g.db + '.chrom.sizes' });
        state.genomeSrc = g;
      } else if (g?.type === 'files') {
        const files = (g.files || []).map(f => pool.get(f.name)).filter(Boolean);
        if (!files.length) throw new Error('The genome file is missing — choose the folder that contains ' + (g.files?.[0]?.name || 'it') + '.');
        await ctx.loadGenomeFromFiles(files);
      } else if (g?.type === 'bam-header') {
        const [bf, bi] = (g.files || []).map(f => pool.get(f.name));
        if (!bf) throw new Error('Missing ' + g.files?.[0]?.name);
        await ctx.genomeFromBam(bf, bi);
      } else if (!state.browser) throw new Error('This session has no genome information.');
      // clear existing tracks
      for (const r of [...state.tracks]) { try { if (r.igv) state.browser.removeTrack(r.igv); } catch {} }
      state.tracks = []; ctx.renderTrackList();
      // tracks, in saved order
      const idMap = {}; const recByIdx = []; const skipped = [];
      for (const [i, t] of (s.tracks || []).entries()) {
        let rec = null;
        try {
          if (t.file) {
            const f = pool.get(t.file.name); if (!f) { skipped.push(t.name || t.file.name); continue; }
            const idx = t.index ? pool.get(t.index.name) : ctx.findIndex(f, [...pool.values()]);
            rec = await ctx.addFileTrack(f, idx || null, { name: t.name, format: t.format, color: t.color, deferList: true, order: i });
          } else if (t.generated) {
            const G = t.generated;
            if (G.type === 'merged') {
              const tracks = (G.members || []).map(mb => mb.ref != null ? (recByIdx[mb.ref] ? { url: recByIdx[mb.ref].file, indexURL: recByIdx[mb.ref].index || undefined, format: mb.format || recByIdx[mb.ref].format, name: mb.name, color: mb.color } : null)
                : { url: new File([mb.text || ''], mb.fileName || 'track.bedgraph'), format: mb.format, name: mb.name, color: mb.color }).filter(Boolean);
              if (tracks.length) {
                const cfg = { type: 'merged', name: t.name, tracks, height: t.height || G.cfg.height || 110, autoscale: G.cfg.autoscale !== false, order: (i + 1) * 10 };
                const it = await state.browser.loadTrack(cfg);
                rec = { id: state.nextId++, name: t.name, format: 'merged', igv: it, color: t.color || '#52514e', height: cfg.height, mode: 'full', generated: true, cfg };
                state.tracks.push(rec);
              }
            } else if (G.text != null) {
              rec = await ctx.addFileTrack(new File([G.text], G.fileName || 'track.bed'), null, { name: t.name, format: G.format || t.format, color: t.color, generated: true, deferList: true, cfg: G.cfg || {} });
            }
          }
        } catch (e) { console.warn(e); skipped.push(t.name); }
        recByIdx[i] = rec;
        if (rec) idMap[t.id] = rec.id;
      }
      // per-track settings
      for (const [i, t] of (s.tracks || []).entries()) {
        const r = recByIdx[i]; if (!r) continue;
        try {
          if (t.sample && r.sample) { Object.assign(r.sample, { use: t.sample.use, isControl: t.sample.isControl, groupEdited: t.sample.groupEdited }); r.sample.group = t.sample.group; r.sample.ready?.then(() => { r.sample.group = t.sample.group; r.sample.use = t.sample.use && !r.sample.err; }); }
          if (t.color && r.igv && t.color !== r.color) ctx.setColor(r, t.color);
          if (r.format === 'bam') {
            if (t.showCoverage === false) ctx.setBamPart(r, 'showCoverage', false);
            if (t.showAlignments === false) ctx.setBamPart(r, 'showAlignments', false);
            if (t.coverageHeight) ctx.setCoverageHeight(r, t.coverageHeight);
            if (t.displayMode && r.igv?.alignmentTrack?.setDisplayMode) { r.igv.alignmentTrack.setDisplayMode(t.displayMode); r.mode = { EXPANDED: 'pack', SQUISHED: 'squish', FULL: 'full', COLLAPSED: 'squish' }[t.displayMode] || r.mode; }
            if (t.sortMode) r.sortMode = t.sortMode;
          }
          if (t.height) ctx.setHeight(r, t.height);
          if (t.compact && t.compact !== 'normal') ctx.setCompact?.(r, t.compact);
          if (t.hideNames) ctx.setFeatureOpt(r, 'gvHideNames', true);
          if (t.hideArrows) ctx.setFeatureOpt(r, 'gvHideArrows', true);
          if (t.mode && t.mode !== r.mode && r.format !== 'bam') await ctx.setMode(r, t.mode);
          if (t.mode === 'hide') await ctx.setMode(r, 'hide');
        } catch (e) { console.warn('track setting', e); }
      }
      ctx.applyOrder();
      // display toggles
      for (const [id, on] of Object.entries(s.toggles || {})) { const el = document.getElementById(id); if (el && el.classList.contains('on') !== !!on) el.click(); }
      // controls
      for (const [id, v] of Object.entries(s.controls || {})) {
        const el = document.getElementById(id); if (!el) continue;
        if (el.type === 'checkbox') el.checked = !!v; else el.value = v;
        if (id === 'pRegion' || id === 'cMode') el.dispatchEvent(new Event('change'));
      }
      // hide rules
      /* hide-indel rules no longer exist: sessions that contain them load without them */
      // SV state
      state.hidden = new Set(s.sv?.hidden || []); state.typeFilter = new Set(s.sv?.typeFilter || []); state.maskIntervals = s.sv?.maskIntervals || [];
      if (state.maskIntervals.length && document.getElementById('fMaskInfo')) document.getElementById('fMaskInfo').textContent = `${state.maskIntervals.length.toLocaleString()} known variants loaded (from session)`;
      // figure styles + annotations
      try { for (const [k, v] of Object.entries(s.figStyles || {})) localStorage.setItem(k, v); } catch {}
      ctx.setAnnotationStore?.(s.annotations || {});
      if (s.dockHeight) document.getElementById('work')?.style.setProperty('--dockH', s.dockHeight);
      if (s.fullPage && !document.getElementById('work')?.classList.contains('full')) document.getElementById('btnFull')?.click();
      // locus
      if (s.locus) await ctx.goto(s.locus);
      ctx.renderTrackList();
      state.crisprSession?.apply(s.crispr, idMap);
      // wait for BAM readers, then re-run analyses
      await Promise.all(state.tracks.filter(t => t.sample?.ready).map(t => t.sample.ready));
      for (const r of state.tracks) if (r.sortMode) ctx.sortReads(r);
      if (rerun && s.sv?.ran && s.sv.region) {
        const reg = document.getElementById('pRegion'); const cus = document.getElementById('pCustom');
        const prevMode = reg.value, prevCustom = cus.value;
        reg.value = 'custom'; reg.dispatchEvent(new Event('change')); cus.value = s.sv.region;
        await ctx.runAnalysis();
        reg.value = prevMode; reg.dispatchEvent(new Event('change')); cus.value = prevCustom;
      }
      if (rerun && s.crispr?.ran) await state.crisprSession?.run();
      if (s.locus) await ctx.goto(s.locus);
      toast(`Session loaded${skipped.length ? ` — ${skipped.length} track(s) skipped (files not found): ${skipped.slice(0, 4).join(', ')}${skipped.length > 4 ? '…' : ''}` : ''}.`, !!skipped.length, 8000);
    } catch (e) { console.error(e); toast('Could not load session: ' + (e.message || e), true, 10000); }
  }

  // ---------------- wiring ----------------
  document.getElementById('btnSession').onclick = () => openSessionDialog();
  // dropping a .genoview.json anywhere opens it
  document.addEventListener('drop', async e => {
    const f = [...(e.dataTransfer?.files || [])].find(x => /\.genoview\.json$|session.*\.json$/i.test(x.name));
    if (!f) return;
    e.stopImmediatePropagation(); e.preventDefault();
    try { openSessionDialog({ sess: migrate(JSON.parse(await f.text())), name: f.name }); } catch (err) { toast(err.message, true); }
  }, true);
  state.session = { collect, applySession, migrate, saveSession };
}
