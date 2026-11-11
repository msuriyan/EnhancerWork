import { yieldNow } from './yield.js';
import { sampleSpectrum, revcomp, collectReads, classify, findGuides, NUCLEASES, alleleAlignment, editsToTSV, allelesToTSV, eventKey, parseSizeFilter } from './crispr.js';
import { indelHistogram, footprintChart, FRAME_COLORS, sampleBarSVG } from './crisprCharts.js';
import { chromAlias } from './analysis.js';
import { renderSpectrum, poolSpecs, DEFAULT_STYLE } from './spectrum.js';
import { openFigureEditor } from './figureEditor.js';
import { renderBars, openBarEditor, barStyle } from './barplot.js';

export function initCrispr(ctx) {
  const { openAnnotator, handleFiles, state, $, $$, toast, download, serializeSVG, svgToPNG, goto, currentLocus, esc, fmt } = ctx;
  const BAR_METRICS = [['readsWithDelPct', '% reads with a deletion', 'Reads with deletion (%)'], ['editedPct', '% edited reads', 'Edited reads (%)'], ['frameshiftPct', '% frameshift reads', 'Frameshift reads (%)'], ['inframePct', '% in-frame reads', 'In-frame reads (%)'], ['largePct', '% reads with deletion ≥ 50 bp', 'Reads with deletion ≥ 50 bp (%)'], ['wtPct', '% unmodified reads', 'Unmodified reads (%)']];
  const cs = {
    sites: [], raw: {}, rawKey: '', results: {}, sample: null, target: null, ref: null, mask: new Map(),
    bar: { metric: 'readsWithDelPct', order: 'input', color: '#4c72b0', byGroup: false, values: true, minDel: 1, maxDel: 0 }, view: { range: 30, bin: 1, basis: 'all', showWT: true, cumulative: false, foot: 40, top: 25 },
  };
  state.crispr = cs;
  cs.use = {};   // track id -> analyse?

  // ---------- BAM picker ----------
  function bamRecs() { return state.tracks.filter(r => r.format === 'bam' && r.sample); }
  function renderBams() {
    const el = $('#cBams'); if (!el) return;
    const recs = bamRecs();
    const status = r => r.sample.err ? `<span style="color:var(--del)">${esc(r.sample.err)}</span>` : !r.sample.bam ? '<span class="muted">reading…</span>' : r.sample.inMemory ? '<span class="muted">no index · in memory</span>' : '<span class="muted">indexed</span>';
    el.innerHTML = `<span class="status"><b>BAMs to analyse:</b></span>` + (recs.length ? recs.map(r => `<label class="sitechip"><input type="checkbox" data-bam="${r.id}" ${cs.use[r.id] !== false && !r.sample.err ? 'checked' : ''} ${r.sample.err ? 'disabled' : ''}> ${esc(r.name)} ${status(r)}</label>`).join('') : '<span class="muted">none loaded</span>') +
      `<button id="cAddBam" title="Choose BAM files (with or without .bai/.csi) to load and analyse">Add BAM…</button>` +
      (recs.length > 1 ? `<button class="ghost" id="cAllBam">all</button><button class="ghost" id="cNoBam">none</button>` : '');
    $$('[data-bam]', el).forEach(i => i.onchange = () => { cs.use[i.dataset.bam] = i.checked; if (Object.keys(cs.raw).length) runOrReclassify(); });
    $('#cAddBam').onclick = () => $('#cBamInput').click();
    if ($('#cAllBam')) { $('#cAllBam').onclick = () => { recs.forEach(r => cs.use[r.id] = true); renderBams(); }; $('#cNoBam').onclick = () => { recs.forEach(r => cs.use[r.id] = false); renderBams(); }; }
  }
  $('#cBamInput').onchange = async e => { const f = [...e.target.files]; e.target.value = ''; await handleFiles(f); renderBams(); };
  state.onSamplesChanged = renderBams;
  // session support
  state.crisprSession = {
    get: () => ({ ran: Object.keys(cs.raw).length > 0, mask: [...cs.mask.entries()], use: { ...cs.use }, view: { ...cs.view }, sample: cs.sample, target: cs.target, studyOpen: !!cs.studyOpen,
      sitesOff: cs.sites.filter(x => !x.on).map(x => x.cut), bar: { ...cs.bar } }),
    apply: (o, idMap) => {
      if (!o) return;
      cs.mask = new Map(o.mask || []);
      cs.use = {}; for (const [k, v] of Object.entries(o.use || {})) if (idMap[k] != null) cs.use[idMap[k]] = v;
      Object.assign(cs.view, o.view || {}); Object.assign(cs.bar, o.bar || {}); cs.studyOpen = !!o.studyOpen;
      cs._restore = { sample: idMap[o.sample] != null ? String(idMap[o.sample]) : null, target: o.target, sitesOff: o.sitesOff || [] };
      renderBams();
    },
    run: async () => {
      await run(true);
      const r = cs._restore; if (!r) return;
      if (r.sitesOff.length) { cs.sites.forEach(x => { if (r.sitesOff.includes(x.cut)) x.on = false; }); renderSites(); markSites(); await runOrReclassify(); }
      if (r.sample && cs.results[r.sample]) cs.sample = r.sample;
      if (r.target && cs.targets?.find(t => t.id === r.target)) cs.target = r.target;
      cs._restore = null; if (Object.keys(cs.results).length) render();
    },
  };
  state.onRulesChanged = () => { if (Object.keys(cs.raw).length) reclassify(); };
  state.onTrackRenamed = (rec) => { const raw = cs.raw[rec.id]; renderBams(); if (raw) { raw.name = rec.name; reclassify(); } };
  renderBams();

  $('#cNuc').innerHTML = Object.entries(NUCLEASES).map(([k, v]) => `<option value="${k}">${v.label}</option>`).join('');
  $('#cReads').onchange = e => { cs.readsSet = true; applyReadPreset(e.target.value); if (Object.keys(cs.raw).length) run(true); };
  function applyReadPreset(v) {
    const long = v === 'long';
    $('#cWin').value = long ? 10 : 5; $('#cMinIndel').value = long ? 3 : 1; $('#cHp').checked = long; $('#cMerge').value = long ? 2 : 0;
    $('#cFlank').value = long ? 10000 : 3000;
  }
  const kindWords = () => $('#cKinds').value === 'D' ? { one: 'deletion', many: 'deletions', Cap: 'Deletion' } : { one: 'indel', many: 'indels', Cap: 'Indel' };
  $('#cKinds').onchange = e => {
    // deletions only: the size distribution is shown as % of deletions (adds up to 100%)
    try { const k = 'genoview.fig.crispr'; const st = JSON.parse(localStorage.getItem(k) || '{}') || {}; st.measure = e.target.value === 'D' ? 'pctEvents' : (st.measure === 'pctEvents' ? 'pct' : st.measure); localStorage.setItem(k, JSON.stringify(st)); } catch {}
    if (Object.keys(cs.raw).length) reclassify();
  };
  $('#cMode').onchange = e => { $('#cROI').style.display = e.target.value === 'roi' ? '' : 'none'; $('#cWinLbl').style.display = e.target.value === 'roi' ? 'none' : ''; if (Object.keys(cs.raw).length) runOrReclassify(); };
  // filters that only need re-classification (no BAM reading)
  for (const id of ['#cWin', '#cMinIndel', '#cHp', '#cMerge', '#cIgnore']) $(id).addEventListener('change', () => { if (Object.keys(cs.raw).length) runOrReclassify(); });
  $('#cROI').addEventListener('change', () => { if (Object.keys(cs.raw).length) runOrReclassify(); });

  async function getSeq(chr, s, e) {
    try { const q = await state.browser.genome.getSequence(chr, Math.max(0, s), e); return q && /[ACGTacgt]/.test(q) ? q : null; } catch { return null; }
  }

  // ---------- sites ----------
  $('#cFind').onclick = () => findSites(true);
  async function findSites(showList) {
    if (!state.browser) { toast('Load a genome first.', true); return []; }
    const t = $('#cTarget').value.trim();
    const tokens = t.split(/[\s;]+|,(?!\d{3}(?!\d))/).filter(Boolean);
    const sites = [];
    const l = currentLocus();
    const guides = [];
    for (const tok of tokens) {
      const m = /^([^:\s]+):([\d,]+)$/.exec(tok);
      if (m) { const chr = state.browser.genome.getChromosomeName?.(m[1]) || m[1]; sites.push({ chr, cut: +m[2].replace(/,/g, ''), strand: '.', start: null, end: null, pam: '', manual: true }); }
      else { const g = tok.toUpperCase().replace(/[^ACGTN]/g, ''); if (g.length >= 15) guides.push(g); }
    }
    if (guides.length) {
      const pad = Math.max(0, Math.min(1_000_000, 2_000_000 - (l.end - l.start)) / 2);
      const s = Math.max(0, Math.floor(l.start - pad)), e = Math.ceil(l.end + pad);
      const seq = await getSeq(l.chr, s, e);
      if (!seq) { toast('No reference sequence loaded — enter cut positions as chr:pos instead.', true); }
      else {
        const miss = [];
        for (const g of guides) {
          const r = resolveGuide(seq, s, g, $('#cNuc').value);
          if (!r.hits.length) miss.push(g.slice(0, 10) + '…');
          for (const h of r.hits) sites.push({ ...h, chr: l.chr, guide: r.guide, note: r.note });
          if (r.note) toast(`${g.slice(0, 12)}…: ${r.note}`, false, 7000);
        }
        if (miss.length) toast(`Not found with a ${NUCLEASES[$('#cNuc').value].pam} PAM within ±1 Mb of the view: ${miss.join(', ')}`, true, 8000);
      }
    }
    if (!sites.length) { if (!tokens.length) toast('Enter one or more guides (space/comma separated) or cut positions like chr1:100,031.', true); renderSites(); return []; }
    // one chromosome only, sorted
    const chr = sites[0].chr;
    const same = sites.filter(x => x.chr === chr).sort((a, b) => a.cut - b.cut);
    if (same.length < sites.length) toast('Sites on other chromosomes are ignored; analyse them separately.', true);
    same.forEach((x, i) => { x.n = i + 1; x.on = true; });
    cs.sites = same; cs.siteKey = t + $('#cNuc').value;
    renderSites(); markSites();
    return same;
  }
  // accept the guide as protospacer (either strand), as its reverse complement, with the PAM attached, or without any PAM match
  function resolveGuide(seq, off, g, nuc) {
    const N = NUCLEASES[nuc]; const P = N.pam.length;
    const IU = { A: 'A', C: 'C', G: 'G', T: 'T', N: '[ACGT]', R: '[AG]', Y: '[CT]', S: '[GC]', W: '[AT]', K: '[GT]', M: '[AC]', B: '[CGT]', D: '[AGT]', H: '[ACT]', V: '[ACG]' };
    const pamRe = new RegExp('^' + N.pam.split('').map(c => IU[c] || c).join('') + '$');
    const hasPam = x => x.length >= 17 + P && pamRe.test(N.side === 3 ? x.slice(-P) : x.slice(0, P));
    const strip = x => N.side === 3 ? x.slice(0, -P) : x.slice(P);
    const tries = [];
    const gr = revcomp(g);
    if (hasPam(g)) tries.push([strip(g), 'PAM removed']);
    if (hasPam(gr)) tries.push([strip(gr), 'reverse complement with PAM — used protospacer ' + strip(gr)]);
    tries.push([g, ''], [gr, 'entered as reverse complement — used the protospacer ' + gr]);
    for (const [gg, note] of tries) { const hits = findGuides(seq, off, gg, nuc); if (hits.length) return { hits, guide: gg, note }; }
    // last resort: sequence match on either strand without a PAM; cut placed as for the chosen nuclease
    for (const [gg, note0] of tries) {
      const S = seq.toUpperCase(); const hits = [];
      let i = -1; while ((i = S.indexOf(gg, i + 1)) >= 0) hits.push({ strand: '+', start: i + off, end: i + gg.length + off, pam: '?', cut: N.side === 3 ? i + gg.length + N.cut + off : i + N.cut + off });
      const rc = revcomp(gg); i = -1; while ((i = S.indexOf(rc, i + 1)) >= 0) hits.push({ strand: '-', start: i + off, end: i + gg.length + off, pam: '?', cut: N.side === 3 ? i - N.cut + off : i + gg.length - N.cut + off });
      if (hits.length) return { hits, guide: gg, note: `no ${N.pam} PAM next to it${note0 ? ' (' + note0 + ')' : ''} — cut placed as for ${N.label}; check the position` };
    }
    return { hits: [], guide: g, note: '' };
  }

  function renderSites() {
    const el = $('#cSites');
    if (!cs.sites.length) { el.innerHTML = ''; return; }
    el.innerHTML = cs.sites.map((x, i) => `<label class="sitechip"><input type="checkbox" data-i="${i}" ${x.on ? 'checked' : ''}> <b>cut ${x.n}</b> <span class="loc">${esc(x.chr)}:${fmt(x.cut)}|${fmt(x.cut + 1)}</span> ${x.manual ? '<span class="muted">manual</span>' : `<span class="muted"><b>${x.strand === '+' ? 'forward (+)' : 'reverse (−)'}</b> ${esc(x.guide.slice(0, 6))}…${esc(x.pam)}</span>${x.note ? ` <span title="${esc(x.note)}" style="color:var(--cnv)">⚠</span>` : ''}`}</label>`).join('') +
      (cs.sites.length > 1 ? `<span class="muted">distance cut 1→${cs.sites.length}: ${fmt(cs.sites[cs.sites.length - 1].cut - cs.sites[0].cut)} bp</span>` : '');
    $$('input', el).forEach(i => i.onchange = () => { cs.sites[+i.dataset.i].on = i.checked; markSites(); if (Object.keys(cs.raw).length) runOrReclassify(); });
  }
  async function markSites() {
    try {
      state.browser.clearROIs();
      const f = [];
      for (const x of cs.sites.filter(s => s.on)) { f.push({ chr: x.chr, start: x.cut - 1, end: x.cut + 1, name: 'cut ' + x.n }); if (x.start != null) f.push({ chr: x.chr, start: x.start, end: x.end, name: `guide ${x.n} (${x.strand})` }); }
      if (f.length) await state.browser.loadROI([{ name: 'CRISPR targets', color: 'rgba(42,120,214,0.15)', features: f }]);
    } catch (e) { console.warn(e); }
  }

  // ---------- targets (what is counted) ----------
  function parseROIs(chr) {
    return $('#cROI').value.split(/[\s;]+|,(?!\d{3}(?!\d))/).filter(Boolean).map(t => {
      const m = /^(?:([^:\s]+):)?([\d,]+)-([\d,]+)$/.exec(t); if (!m) return null;
      return { chr: m[1] || chr, s: +m[2].replace(/,/g, '') - 1, e: +m[3].replace(/,/g, '') };
    }).filter(Boolean);
  }
  function buildTargets() {
    const on = cs.sites.filter(s => s.on);
    const W = +$('#cWin').value;
    const cutsAll = on.map(s => s.cut);
    const T = [];
    if ($('#cMode').value === 'roi') {
      const rois = parseROIs(on[0]?.chr || currentLocus().chr).filter(r => r.e > r.s);
      if (!rois.length) { toast('Enter regions of interest like chr1:100,020-100,045 (several allowed).', true); return []; }
      rois.forEach((r, i) => {
        const cuts = cutsAll.filter(c => c >= r.s && c <= r.e);
        T.push({ id: 'roi' + i, label: `ROI ${i + 1} ${fmt(r.s + 1)}-${fmt(r.e)}`, cuts, areas: [[r.s, r.e]], origin: cuts[0] ?? r.s, chr: r.chr });
      });
      if (rois.length > 1) T.push({ id: 'roiAll', label: 'all ROIs combined', cuts: cutsAll.filter(c => rois.some(r => c >= r.s && c <= r.e)), areas: rois.map(r => [r.s, r.e]), origin: cutsAll[0] ?? rois[0].s, combined: true });
    } else {
      on.forEach(s => T.push({ id: 'cut' + s.n, label: `cut ${s.n}`, cuts: [s.cut], areas: [[s.cut - W, s.cut + W]], origin: s.cut }));
      if (on.length > 1) T.push({ id: 'all', label: `all ${on.length} cuts combined`, cuts: cutsAll, areas: cutsAll.map(c => [c - W, c + W]), origin: cutsAll[0], combined: true });
    }
    return T;
  }

  function classifyOpts() {
    return { kinds: $('#cKinds').value, minIndel: +$('#cMinIndel').value, hpFilter: $('#cHp').checked, mergeTol: +$('#cMerge').value, largeDel: 50,
      ignoreSize: parseSizeFilter($('#cIgnore').value), mask: new Set(cs.mask.keys()), ref: cs.ref,
      ruleMask: state.hideRulesInAnalyses ? (chr => e => state.ruleHidesIndel(chr, e.pos, e.len, e.kind))(cs.sites[0]?.chr || currentLocus()?.chr) : null };
  }

  // ---------- run ----------
  $('#cRun').onclick = () => run(true);
  async function runOrReclassify() {
    const need = regionNeeded(); if (!need) return;
    if (cs.rawKey === need.key) reclassify(); else run(false);
  }
  function regionNeeded() {
    const on = cs.sites.filter(s => s.on);
    const flank = +$('#cFlank').value;
    let lo = Math.min(...on.map(s => s.cut)), hi = Math.max(...on.map(s => s.cut));
    if ($('#cMode').value === 'roi') for (const r of parseROIs(on[0]?.chr)) { lo = Math.min(lo, r.s); hi = Math.max(hi, r.e); }
    if (!isFinite(lo)) return null;
    const chr = on[0]?.chr || currentLocus().chr;
    return { chr, s: Math.max(0, lo - flank), e: hi + flank, key: `${chr}:${Math.max(0, lo - flank)}-${hi + flank}:${$('#cMapq').value}:${bamRecs().filter(r => cs.use[r.id] !== false && r.sample.bam).map(r => r.id).join(',')}` };
  }

  async function run(force) {
    if (!state.browser) return toast('Load a genome and BAMs first.', true);
    if (force || !cs.sites.length || cs.siteKey !== $('#cTarget').value.trim() + $('#cNuc').value) { if (!(await findSites(true)).length && $('#cMode').value !== 'roi') return; }
    if (!cs.sites.some(s => s.on) && $('#cMode').value !== 'roi') return toast('Tick at least one cut site.', true);
    const chosen = bamRecs().filter(r => cs.use[r.id] !== false && !r.sample.err);
    if (!bamRecs().length) return toast('No BAM loaded yet — use “Add BAM…” here or drop BAM files (index optional) on the page.', true, 7000);
    if (!chosen.length) return toast('Tick at least one BAM under “BAMs to analyse”.', true);
    if (chosen.some(r => !r.sample.bam)) { $('#cStatus').textContent = 'Waiting for BAMs to finish reading…'; await Promise.all(chosen.map(r => r.sample.ready)); }
    const recs = chosen.filter(r => r.sample.bam);
    if (!recs.length) return toast('The selected BAMs could not be read: ' + chosen.map(r => r.name + ' — ' + (r.sample.err || 'unknown error')).join('; '), true, 9000);
    const need = regionNeeded(); if (!need) return toast('No cut site or region of interest given.', true);
    if (need.e - need.s > 2_000_000) return toast('Cut sites / regions span > 2 Mb; analyse them in smaller groups.', true);
    $('#cRun').disabled = true;
    try {
      const seq = await getSeq(need.chr, need.s - 30, need.e + 30);
      cs.ref = seq ? { seq, start: Math.max(0, need.s - 30) } : null;
      cs.raw = {};
      for (const [i, r] of recs.entries()) {
        $('#cStatus').textContent = `Reading ${r.name} (${i + 1}/${recs.length})…`;
        await yieldNow(true);
        const bamChr = chromAlias(need.chr, r.sample.chroms);
        cs.raw[r.id] = bamChr ? await collectReads({ name: r.name, bam: r.sample.bam, bamChr }, { regStart: need.s, regEnd: need.e, minMapq: +$('#cMapq').value })
          : { name: r.name, error: `no ${need.chr} in BAM` };
      }
      cs.rawKey = need.key;
      // first run: detect short vs long reads from aligned read lengths (median > 1 kb = long)
      if (!cs.readsSet) {
        cs.readsSet = true;
        const lens = Object.values(cs.raw).flatMap(r => (r.reads || []).slice(0, 2000).map(x => x.s1 - x.s0)).sort((a, b) => a - b);
        const guess = lens.length ? (lens[lens.length >> 1] > 1000 ? 'long' : 'short') : null;
        if (guess && guess !== $('#cReads').value) {
          $('#cReads').value = guess; applyReadPreset(guess);
          toast(`${guess === 'long' ? 'Long' : 'Short'} reads detected — ${guess === 'long' ? 'long-read settings applied (± 10 bp window, homopolymer filter, ± 2 bp merge, 10 kb flank)' : 'short-read settings applied'}.`);
          $('#cRun').disabled = false; return run(true);
        }
      }
      reclassify();
    } catch (e) { console.error(e); toast('CRISPR analysis failed: ' + e.message, true, 9000); }
    finally { $('#cRun').disabled = false; }
  }

  function reclassify() {
    const targets = buildTargets(); if (!targets.length) return;
    const o = classifyOpts();
    cs.targets = targets; cs.results = {};
    const t0 = performance.now();
    const single = targets.filter(t => !t.combined);
    const allAreas = single.flatMap(t => t.areas), allCuts = [...new Set(single.flatMap(t => t.cuts))];
    cs.spec = {};
    for (const [id, raw] of Object.entries(cs.raw)) {
      cs.results[id] = raw.error ? { error: raw.error, name: raw.name } : Object.fromEntries(targets.map(t => [t.id, classify(raw, t, o)]));
      if (!raw.error) {
        cs.spec[id] = { __all: sampleSpectrum(raw, allAreas, allCuts, o) };
        for (const t of targets) cs.spec[id][t.id] = sampleSpectrum(raw, t.areas, t.cuts, o);
      }
    }
    // all samples combined (reads pooled) — the study-wide view
    const okIds = Object.keys(cs.raw).filter(id => !cs.raw[id].error);
    if (okIds.length > 1) {
      const raws = okIds.map(id => cs.raw[id]);
      const pr = { name: `All ${okIds.length} samples combined`, regStart: raws[0].regStart, regEnd: raws[0].regEnd, total: raws.reduce((a, r) => a + r.total, 0), lowMapq: raws.reduce((a, r) => a + r.lowMapq, 0), reads: raws.flatMap(r => r.reads) };
      cs.poolName = pr.name;
      cs.results.__pool = Object.fromEntries(targets.map(t => [t.id, classify(pr, t, o)]));
    } else cs.poolName = null;
    const ids = Object.keys(cs.results);
    if (!cs.sample || !cs.results[cs.sample]) cs.sample = cs.results.__pool ? '__pool' : ids[0];
    if (!cs.target || !targets.find(t => t.id === cs.target)) cs.target = targets[0].id;
    const nOn = cs.sites.filter(s => s.on).length;
    $('#cStatus').textContent = `${targets.length} target${targets.length > 1 ? 's' : ''} · ${nOn} cut${nOn !== 1 ? 's' : ''} · ${ids.filter(i => i !== '__pool').length} sample${ids.filter(i => i !== '__pool').length > 1 ? 's' : ''}${cs.results.__pool ? ' + combined' : ''} · recalculated in ${Math.round(performance.now() - t0)} ms`;
    render();
  }

  // ---------- size spectrum (whole sample / study) ----------
  const scopeLabel = sc => sc === '__all' ? (cs.targets.filter(t => !t.combined).length > 1 ? `all ${cs.targets.filter(t => !t.combined).length} ${$('#cMode').value === 'roi' ? 'regions' : 'cuts'} pooled` : cs.targets[0]?.label || 'all') : cs.targets.find(t => t.id === sc)?.label || sc;
  function getSpec(sel) {
    const ids = Object.keys(cs.spec);
    const sc = sel.scope || '__all';
    if (sel.sample === '__pool' || sel.sample === '__mean') {
      const specs = ids.map(id => cs.spec[id][sc]).filter(Boolean);
      const sp = poolSpecs(specs, sel.sample === '__pool' ? 'reads' : 'mean', `${specs.length} samples`);
      sp.subtitle = sel.sample === '__pool'
        ? `${specs.length} samples, reads pooled · n = ${sp.total.toLocaleString('en-US')} reads (${sp.edited.toLocaleString('en-US')} with ${kindWords().many}) · ${scopeLabel(sc)}`
        : `Mean of ${specs.length} samples (each weighted equally) · ${scopeLabel(sc)}`;
      return sp;
    }
    if (!cs.spec[sel.sample]) return getSpec({ ...sel, sample: '__pool' });
    const sp = { ...cs.spec[sel.sample][sc] };
    sp.subtitle = `${sp.name} · n = ${sp.total.toLocaleString('en-US')} reads (${sp.edited.toLocaleString('en-US')} with ${kindWords().many}) · ${scopeLabel(sc)}`;
    return sp;
  }
  function editorSelects(sample, scope) {
    const ids = Object.keys(cs.spec);
    const sOpts = ids.map(id => ({ value: id, label: cs.raw[id].name }));
    if (ids.length > 1) sOpts.unshift({ value: '__pool', label: `All ${ids.length} samples combined — reads pooled` }, { value: '__mean', label: `All ${ids.length} samples — mean ± SD of samples` });
    const scOpts = [{ value: '__all', label: 'Whole sample (all cuts / regions)' }, ...cs.targets.map(t => ({ value: t.id, label: t.label }))];
    return [{ id: 'sample', label: 'Sample', options: sOpts, value: sample }, { id: 'scope', label: 'Indels from', options: scOpts, value: scope }];
  }
  function openEditor(sample, scope = '__all') {
    openFigureEditor({ key: 'crispr', title: `${kindWords().Cap} size spectrum — figure editor`, selects: editorSelects(sample, scope), getSpec,
      defaults: { xLabel: 'Indel length (bp)' },
      filename: sel => `${kindWords().one}_spectrum_${sel.sample === '__mean' ? 'study_mean' : sel.sample === '__pool' ? 'study_pooled' : (cs.raw[sel.sample]?.name || 'sample').replace(/\.bam$/i, '')}_${sel.scope === '__all' ? 'all' : sel.scope}` });
  }
  function loadSavedStyle() { try { return { ...DEFAULT_STYLE, ...(JSON.parse(localStorage.getItem('genoview.fig.crispr') || 'null') || {}) }; } catch { return { ...DEFAULT_STYLE }; } }
  function barValue(r, d) {
    if (d.metric === 'readsWithDelPct' && ((d.minDel || 1) > 1 || d.maxDel > 0)) {
      const lo = d.minDel || 1, hi = d.maxDel > 0 ? d.maxDel : Infinity, n = r.stats.spanning || 0;
      if (!n) return 0;
      return 100 * (r.alleleList || []).filter(a => (a.events || []).some(e => e.kind === 'D' && e.len >= lo && e.len <= hi)).reduce((s, a) => s + a.count, 0) / n;
    }
    return r.stats[d.metric] || 0;
  }
  function barRows(d = cs.bar) {
    const ids = Object.keys(cs.results).filter(id => id !== '__pool' && !cs.results[id].error && cs.results[id][cs.target]);
    let rows = ids.map(id => { const r = cs.results[id][cs.target]; const rec = state.tracks.find(t => String(t.id) === String(id));
      return { id: cs.raw[id].name, name: cs.raw[id].name.replace(/\.bam$/i, ''), value: barValue(r, d), n: r.stats.spanning, group: rec?.sample?.group || '' }; });
    const o = d.order || 'input';
    if (o === 'value') rows.sort((a, b) => b.value - a.value); else if (o === 'valueAsc') rows.sort((a, b) => a.value - b.value); else if (o === 'name') rows.sort((a, b) => a.name.localeCompare(b.name)); else if (o === 'group') rows.sort((a, b) => a.group.localeCompare(b.group) || a.name.localeCompare(b.name));
    return rows;
  }
  const barLabel = d => d.metric === 'readsWithDelPct' && ((d.minDel || 1) > 1 || d.maxDel > 0) ? `Reads with deletion ${d.maxDel > 0 ? `${d.minDel || 1}–${d.maxDel}` : '≥ ' + d.minDel} bp (%)` : BAR_METRICS.find(m => m[0] === d.metric)[2];
  // inline preview uses the editor's saved style, with the quick controls on top
  function barSVG() {
    const S = barStyle(); S.color = cs.bar.color; S.showValues = cs.bar.values; S.colorMode = cs.bar.byGroup ? 'group' : 'single';
    return renderBars(barRows(), S, { ylabel: barLabel(cs.bar) });
  }
  function barCardHTML() {
    const n = Object.keys(cs.results).filter(id => id !== '__pool' && !cs.results[id].error).length;
    if (n < 2) return '';
    const b = cs.bar;
    return `<div class="card2" style="margin-bottom:12px"><div class="graphhead" style="margin-bottom:4px"><h5 style="margin:0">Per-sample bar plot · ${esc(cs.targets.find(t => t.id === cs.target)?.label || '')}</h5>
      <label class="status">value <select id="cbMetric">${BAR_METRICS.map(([v, l]) => `<option value="${v}" ${v === b.metric ? 'selected' : ''}>${l}</option>`).join('')}</select></label>
      <label class="status" title="Count only reads whose deletion is within this size range (0 = no upper limit)" ${b.metric === 'readsWithDelPct' ? '' : 'hidden'}>deletion ≥ <input type="number" id="cbMin" value="${b.minDel || 1}" min="1" style="width:60px"> ≤ <input type="number" id="cbMax" value="${b.maxDel || 0}" min="0" style="width:60px"> bp</label>
      <label class="status">order <select id="cbOrder">${[['input', 'as loaded'], ['value', 'high → low'], ['valueAsc', 'low → high'], ['name', 'by name'], ['group', 'by group']].map(([v, l]) => `<option value="${v}" ${v === b.order ? 'selected' : ''}>${l}</option>`).join('')}</select></label>
      <label class="status">colour <input type="color" id="cbColor" value="${b.color}"></label>
      <label class="status" title="Colour bars by the BAM's group (set in the track list)"><input type="checkbox" id="cbGroup" ${b.byGroup ? 'checked' : ''}> by group</label>
      <label class="status"><input type="checkbox" id="cbVals" ${b.values ? 'checked' : ''}> values</label>
      <span style="flex:1"></span><button id="cbTsv">TSV</button><button id="cbEdit" class="primary" title="Full editor: hide/rename x-axis names, colours, fonts, sizes, axis, text boxes; export PNG at any dpi / SVG / PDF">Edit &amp; export…</button></div>
      <div id="cbChart" class="chart" style="overflow-x:auto"></div></div>`;
  }
  function drawBars() { const el = $('#cbChart'); if (!el) return; el.innerHTML = barSVG(); }
  function bindBars() {
    if (!$('#cbChart')) return;
    const upd = () => drawBars();
    $('#cbMetric').onchange = e => { cs.bar.metric = e.target.value; $('#cbMin').parentElement.hidden = cs.bar.metric !== 'readsWithDelPct'; upd(); };
    $('#cbMin').onchange = e => { cs.bar.minDel = Math.max(1, +e.target.value || 1); upd(); };
    $('#cbMax').onchange = e => { cs.bar.maxDel = Math.max(0, +e.target.value || 0); upd(); };
    $('#cbOrder').onchange = e => { cs.bar.order = e.target.value; upd(); };
    $('#cbColor').oninput = e => { cs.bar.color = e.target.value; upd(); };
    $('#cbGroup').onchange = e => { cs.bar.byGroup = e.target.checked; upd(); };
    $('#cbVals').onchange = e => { cs.bar.values = e.target.checked; upd(); };
    $('#cbTsv').onclick = () => download(`per_sample_${cs.bar.metric}.tsv`, 'sample\tgroup\tspanning_reads\t' + barLabel(cs.bar).replace(/\s+/g, '_') + '\n' + barRows().map(r => `${r.name}\t${r.group}\t${r.n}\t${r.value.toFixed(4)}`).join('\n') + '\n', 'text/tab-separated-values');
    $('#cbEdit').onclick = () => {
      // seed the editor style from the quick controls
      try { const S = barStyle(); S.color = cs.bar.color; S.showValues = cs.bar.values; S.colorMode = cs.bar.byGroup ? 'group' : 'single'; localStorage.setItem('genoview.fig.bars', JSON.stringify(S)); } catch {}
      openBarEditor({ metrics: BAR_METRICS, data: { metric: cs.bar.metric, minDel: cs.bar.minDel || 1, maxDel: cs.bar.maxDel || 0, order: cs.bar.order },
        getRows: d => barRows(d),
        onData: d => { Object.assign(cs.bar, d); const m = $('#cbMetric'); if (m) { m.value = d.metric; $('#cbMin').value = d.minDel || 1; $('#cbMax').value = d.maxDel || 0; $('#cbOrder').value = d.order; $('#cbMin').parentElement.hidden = d.metric !== 'readsWithDelPct'; } drawBars(); },
        filename: `per_sample_${cs.bar.metric}${cs.bar.minDel > 1 ? '_ge' + cs.bar.minDel + 'bp' : ''}` });
    };
    // when the editor is closed, refresh the inline preview with its style
    if (!cs._barCloseHook) { cs._barCloseHook = true; document.addEventListener('click', e => { if (e.target?.id === 'beClose' || e.target?.id === 'beReset') setTimeout(drawBars, 0); }); }
    drawBars();
  }

  function studyTable(all) {
    const ids = Object.keys(cs.results).filter(id => id !== '__pool' && !cs.results[id].error && cs.results[id][cs.target]);
    if (ids.length < 2) return '';
    const keys = [['editedPct', 'Edited'], ['frameshiftPct', 'Frameshift'], ['inframePct', 'In-frame'], ['wtPct', 'Unmodified'], ['largePct', 'Large ≥50'], ...(all[cs.target].cuts.length > 1 ? [['dropoutPct', 'Dropout']] : []), ['alleles5', 'Alleles ≥5%']];
    const vals = ids.map(id => cs.results[id][cs.target].stats);
    const mean = k => vals.reduce((a, v) => a + v[k], 0) / vals.length;
    const sd = k => { const m = mean(k); return Math.sqrt(vals.reduce((a, v) => a + (v[k] - m) ** 2, 0) / Math.max(1, vals.length - 1)); };
    const f = (k, v) => k === 'alleles5' ? (+v).toFixed(k === 'alleles5' && v % 1 ? 1 : 0) : v.toFixed(1) + '%';
    return `<details class="card2" style="margin-bottom:10px" ${cs.studyOpen ? 'open' : ''} id="cStudy"><summary><b>Study summary · ${ids.length} samples</b> <span class="muted">(${esc(all[cs.target].target.label)}) · mean ± SD across samples</span>
      ${keys.slice(0, 4).map(([k, l]) => `<span class="pill" style="margin-left:6px">${l} ${mean(k).toFixed(1)}${k === 'alleles5' ? '' : '%'} ± ${sd(k).toFixed(1)}</span>`).join('')}
      <button id="cStudyFig" class="primary" style="margin-left:10px">Study size spectrum…</button></summary>
      <div style="overflow-x:auto;margin-top:6px"><table class="calls"><thead><tr><th>Sample</th><th class="num">Reads</th>${keys.map(([, l]) => `<th class="num">${l}</th>`).join('')}</tr></thead><tbody>
      ${ids.map((id, i) => `<tr data-samp="${id}" class="${id === cs.sample ? 'sel' : ''}"><td>${esc(cs.raw[id].name)}</td><td class="num">${vals[i].spanning.toLocaleString('en-US')}</td>${keys.map(([k]) => `<td class="num">${f(k, vals[i][k])}</td>`).join('')}</tr>`).join('')}
      <tr style="font-weight:600"><td>Mean ± SD</td><td class="num">${Math.round(mean('spanning')).toLocaleString('en-US')}</td>${keys.map(([k]) => `<td class="num">${mean(k).toFixed(1)} ± ${sd(k).toFixed(1)}</td>`).join('')}</tr></tbody></table></div></details>`;
  }

  // ---------- report ----------
  function render() {
    const el = $('#cReport');
    const ids = Object.keys(cs.results);
    const all = cs.results[cs.sample];
    const v = cs.view;
    const nameOf = id => id === '__pool' ? `★ ${cs.poolName} (reads pooled)` : cs.raw[id].name;
    const sampOpts = [...ids.filter(i => i === '__pool'), ...ids.filter(i => i !== '__pool')].map(id => `<option value="${id}" ${id === cs.sample ? 'selected' : ''}>${esc(nameOf(id))}</option>`).join('');
    if (all.error) { el.innerHTML = `<div class="sampbar"><label class="status">Sample <select id="cSamp">${sampOpts}</select></label></div><div class="empty">${esc(all.error)}</div>`; $('#cSamp').onchange = e => { cs.sample = e.target.value; render(); }; return; }
    const res = all[cs.target];
    const st = res.stats; const p = x => x.toFixed(1);
    if (!res.spanning) {
      el.innerHTML = `<div class="sampbar"><label class="status">Sample <select id="cSamp">${sampOpts}</select></label></div><div class="empty">No reads in <b>${esc(res.name)}</b> span ${esc(res.target.label)} (${fmt(res.total)} reads at the locus, ${fmt(res.notSpanning)} end inside the counting region). Use a smaller window / region, or check the coordinates.</div>`;
      $('#cSamp').onchange = e => { cs.sample = e.target.value; render(); };
      if (cs.targets.length > 1) { /* still allow choosing another target */ el.insertAdjacentHTML('beforeend', `<div class="chips">${cs.targets.map(t => `<button class="chip" data-t2="${t.id}">${esc(t.label)}</button>`).join('')}</div>`); $$('[data-t2]').forEach(b => b.onclick = () => { cs.target = b.dataset.t2; render(); }); }
      return;
    }
    const multi = cs.targets.length > 1;
    const tRows = cs.targets.map(t => { const r = all[t.id]; const s2 = r.stats; const top = r.alleleList.find(a => a.cls !== 'wt');
      return `<tr data-t="${t.id}" class="${t.id === cs.target ? 'sel' : ''}"><td><b>${esc(t.label)}</b>${t.combined ? ' <span class="muted" title="Reads must span every counting window">(reads spanning all)</span>' : ''}</td><td class="num">${fmt(s2.spanning)}</td><td class="num"><b>${p(s2.editedPct)}%</b></td><td class="num">${p(s2.frameshiftPct)}%</td><td class="num">${p(s2.inframePct)}%</td><td class="num">${p(s2.largePct)}%</td>${multi ? `<td class="num">${s2.dropoutPct ? p(s2.dropoutPct) + '%' : '–'}</td>` : ''}<td class="num">${s2.alleles5}</td><td>${top ? esc(top.label) + ` <span class="muted">${top.pct.toFixed(1)}%</span>` : '<span class="muted">–</span>'}</td></tr>`; }).join('');
    const maskChips = [...cs.mask.entries()].map(([k, lab]) => `<span class="chip" style="background:var(--line-soft)">${esc(lab)} <button class="ghost x" data-unmask="${esc(k)}" title="Count again">×</button></span>`).join('');
    const filt = res.masked + res.filteredSize + res.filteredHp + res.filteredNoise;
    el.innerHTML = `
      <div class="sampbar"><label class="status">Sample <select id="cSamp">${sampOpts}</select></label>
        <span class="status">${fmt(res.total)} reads at the locus · ${fmt(st.spanning)} span ${res.target.combined ? 'all windows' : 'the counting window'}${res.notSpanning ? ` · ${fmt(res.notSpanning)} end inside it (ignored)` : ''}${res.lowMapq ? ` · ${fmt(res.lowMapq)} low MAPQ` : ''}${filt ? ` · ignored indel events: ${[res.masked && res.masked + ' masked', res.filteredSize && res.filteredSize + ' by size', (res.filteredHp + res.filteredNoise) && (res.filteredHp + res.filteredNoise) + ' as seq. error'].filter(Boolean).join(', ')}` : ''}</span>
        <span style="flex:1"></span><button id="cGoCut">Go to target</button><button id="cExHist">Histogram TSV</button><button id="cExAll">Alleles TSV</button></div>
      ${studyTable(all)}
      ${barCardHTML()}
      ${cs.targets.length > 1 ? `<div style="overflow-x:auto;margin-bottom:8px"><table class="calls" id="cTargets"><thead><tr><th>Target (click to view)</th><th class="num">Reads</th><th class="num">Edited</th><th class="num">Frameshift</th><th class="num">In-frame</th><th class="num">Large ≥50</th>${multi ? '<th class="num">Dropout</th>' : ''}<th class="num">Alleles ≥5%</th><th>Top edited allele</th></tr></thead><tbody>${tRows}</tbody></table>
        ${cs.targets.some(t => t.combined) ? '<div class="hint">Combined targets count only reads that cross every window; with reads short relative to the cut distance, dropout alleles are over-represented — prefer per-cut numbers or long reads for frequencies.</div>' : ''}</div>` : ''}
      ${cs.mask.size ? `<div class="chips"><span class="status">Hidden as polymorphism / not of interest:</span>${maskChips}<button id="cUnmaskAll">show all</button></div>` : ''}
      <div class="tiles">
        <div class="tile"><div class="v">${p(st.editedPct)}<small>%</small></div><div class="l">edited</div></div>
        <div class="tile"><div class="v" style="color:#c2501f">${p(st.frameshiftPct)}<small>%</small></div><div class="l">frameshift</div></div>
        <div class="tile"><div class="v" style="color:#1c5cab">${p(st.inframePct)}<small>%</small></div><div class="l">in-frame</div></div>
        <div class="tile"><div class="v">${p(st.wtPct)}<small>%</small></div><div class="l">unmodified</div></div>
        <div class="tile"><div class="v">${p(st.largePct)}<small>%</small></div><div class="l">large deletion ≥ 50 bp</div></div>
        ${res.cuts.length > 1 ? `<div class="tile"><div class="v">${p(st.dropoutPct)}<small>%</small></div><div class="l">dropout between cuts</div></div>` : ''}
        ${st.complexPct > 0 ? `<div class="tile"><div class="v">${p(st.complexPct)}<small>%</small></div><div class="l">complex / unresolved</div></div>` : ''}
        <div class="tile"><div class="v">${st.alleles5} <small>/ ${st.alleles1}</small></div><div class="l">edited alleles ≥ 5% / ≥ 1%</div></div>
        <div class="tile"><div class="v">${st.meanNet > 0 ? '+' : ''}${st.meanNet.toFixed(1)}<small> bp</small></div><div class="l">mean net indel (edited)</div></div>
      </div>
      <div class="card2" style="margin-bottom:12px"><div class="graphhead" style="margin-bottom:4px"><h5 style="margin:0">${kindWords().Cap} size spectrum${$('#cKinds').value === 'D' ? ' (% of deletions, sums to 100%)' : ''} · whole sample (all ${$('#cMode').value === 'roi' ? 'regions' : 'cuts'} pooled), small → large</h5>
        <span class="status">every indel counted once per read · deletions and insertions separately</span><span style="flex:1"></span>
        <button id="cSpecEdit" class="primary">Edit figure &amp; export…</button></div>
        <div id="cSpec" class="chart" style="overflow-x:auto"></div></div>
      <div class="crgrid">
        <div class="card2"><div class="graphhead" style="margin-bottom:2px"><h5 style="margin:0">Net indel per read · ${esc(res.target.label)}</h5>
          <label class="status">y <select id="cvBasis"><option value="all">% of spanning reads</option><option value="edited">% of edited reads</option></select></label>
          <label class="status">range ± <select id="cvRange">${[10, 20, 30, 50, 100, 200, 500].map(x => `<option ${x === v.range ? 'selected' : ''}>${x}</option>`).join('')}</select></label>
          <label class="status">bin <select id="cvBin">${[1, 3, 5, 10].map(x => `<option ${x === v.bin ? 'selected' : ''}>${x}</option>`).join('')}</select> bp</label>
          <label class="status"><input type="checkbox" id="cvWT" ${v.showWT ? 'checked' : ''}> unmodified bar</label>
          <label class="status"><input type="checkbox" id="cvCum" ${v.cumulative ? 'checked' : ''}> cumulative</label>
          <span style="flex:1"></span><button id="cHSVG">SVG</button><button id="cHPNG">PNG</button><button id="cHAnn">Annotate…</button></div>
          <div id="cHist" class="chart"></div></div>
        <div class="card2"><div class="graphhead" style="margin-bottom:2px"><h5 style="margin:0">Deletion footprint</h5>
          <label class="status">± <select id="cvFoot">${[20, 40, 100, 300, 1000, 3000, 10000].filter(x => x <= +$('#cFlank').value).map(x => `<option ${x === v.foot ? 'selected' : ''}>${x}</option>`).join('')}</select> bp</label>
          <span style="flex:1"></span><button id="cFSVG">SVG</button><button id="cFPNG">PNG</button><button id="cFAnn">Annotate…</button></div>
          <div id="cFoot" class="chart"></div></div>
      </div>
      <div class="card2" style="margin-top:12px"><div class="graphhead" style="margin-bottom:4px"><h5 style="margin:0">Alleles · ${esc(res.target.label)}</h5><span class="status">reference (grey) over read allele · <span class="seqcell"><span class="d">-</span></span> deleted · <span class="seqcell"><span class="i">acg</span></span> inserted · <span class="seqcell"><span class="c">|</span></span> cut · click <b>⊘</b> on an indel to hide it as a polymorphism / not of interest</span>
          <span style="flex:1"></span><label class="status">show <select id="cvTop">${[10, 25, 50, 200].map(x => `<option ${x === v.top ? 'selected' : ''}>${x}</option>`).join('')}</select></label></div>
        <div style="overflow-x:auto"><table class="calls"><thead><tr><th class="num">#</th><th class="num">% reads</th><th class="num">reads</th><th>Allele (indel events)</th><th>Frame</th><th>Sequence around cut${res.cuts.length > 1 ? 's' : ''} (±25 bp)</th></tr></thead><tbody id="cAlleles"></tbody></table></div></div>`;

    $('#cvBasis').value = v.basis;
    $('#cSpecEdit').onclick = () => openEditor(cs.sample, '__all');
    bindBars();
    if ($('#cStudyFig')) $('#cStudyFig').onclick = e => { e.preventDefault(); openEditor('__mean', '__all'); };
    if ($('#cStudy')) { $('#cStudy').ontoggle = () => { cs.studyOpen = $('#cStudy').open; }; $$('#cStudy tr[data-samp]').forEach(tr => tr.onclick = () => { cs.sample = tr.dataset.samp; render(); }); }
    $('#cSamp').onchange = e => { cs.sample = e.target.value; render(); };
    $$('#cTargets tr[data-t]').forEach(tr => tr.onclick = () => { cs.target = tr.dataset.t; render(); });
    $$('[data-unmask]').forEach(b => b.onclick = () => { cs.mask.delete(b.dataset.unmask); reclassify(); });
    if ($('#cUnmaskAll')) $('#cUnmaskAll').onclick = () => { cs.mask.clear(); reclassify(); };
    $('#cGoCut').onclick = () => { const lo = Math.min(...res.areas.map(a => a[0])), hi = Math.max(...res.areas.map(a => a[1])); goto(`${cs.sites[0]?.chr || currentLocus().chr}:${fmt(Math.max(1, lo - 60))}-${fmt(hi + 60)}`); };
    const base = `${res.name}_${res.target.label.replace(/[^\w.-]+/g, '_')}`;
    $('#cExHist').onclick = () => download(`${base}_indel_distribution.tsv`, editsToTSV(res), 'text/tab-separated-values');
    $('#cExAll').onclick = () => download(`${base}_alleles.tsv`, allelesToTSV(res), 'text/tab-separated-values');
    const bind = (id, key, conv = x => x) => $(id).onchange = e => { v[key] = conv(e.target.type === 'checkbox' ? e.target.checked : e.target.value); draw(); };
    bind('#cvBasis', 'basis'); bind('#cvRange', 'range', Number); bind('#cvBin', 'bin', Number); bind('#cvWT', 'showWT'); bind('#cvCum', 'cumulative'); bind('#cvFoot', 'foot', Number); bind('#cvTop', 'top', Number);
    $('#cHSVG').onclick = () => download(`${base}_indel_histogram.svg`, serializeSVG($('#cHist svg')), 'image/svg+xml');
    $('#cHPNG').onclick = () => svgToPNG($('#cHist svg'), `${base}_indel_histogram.png`);
    $('#cFSVG').onclick = () => download(`${base}_deletion_footprint.svg`, serializeSVG($('#cFoot svg')), 'image/svg+xml');
    $('#cFPNG').onclick = () => svgToPNG($('#cFoot svg'), `${base}_deletion_footprint.png`);
    const ann = (sel, key, title, fn) => () => { const el = $(sel); if (!el) return; openAnnotator({ svg: serializeSVG(el).replace(/^<\?xml[^>]*>\s*/, ''), key, title, filename: fn }); };
    $('#cHAnn').onclick = ann('#cHist svg', 'crispr-hist', 'Annotate & export · net indel histogram', `${base}_indel_histogram`);
    $('#cFAnn').onclick = ann('#cFoot svg', 'crispr-foot', 'Annotate & export · deletion footprint', `${base}_deletion_footprint`);
    draw();
  }

  const evLabel = (e, chr) => e.kind === 'D' ? `−${e.len} bp del at ${chr}:${fmt(e.pos + 1)}` : `+${e.len} bp ins${e.seq && e.len <= 10 ? ' ' + e.seq : ''} at ${chr}:${fmt(e.pos)}^${fmt(e.pos + 1)}`;
  function draw() {
    const all = cs.results[cs.sample]; if (!all || all.error) return;
    const res = all[cs.target]; const v = cs.view;
    const chr = cs.sites[0]?.chr || currentLocus().chr;
    indelHistogram($('#cHist'), res, v);
    if ($('#cSpec') && (cs.spec[cs.sample] || cs.sample === '__pool')) {
      const st = loadSavedStyle(); const w = Math.max(480, Math.min(1100, $('#cSpec').clientWidth || 800));
      const sp = getSpec({ sample: cs.sample, scope: '__all' });
      $('#cSpec').innerHTML = renderSpectrum(sp, { ...st, width: w, height: Math.round(w * 0.42) }).svg;
    }
    footprintChart($('#cFoot'), res, cs.sites.filter(s => s.on), +$('#cvFoot').value || v.foot);
    const rows = res.alleleList.slice(0, v.top);
    const max = Math.max(...rows.map(a => a.pct), 1);
    $('#cAlleles').innerHTML = rows.map((a, i) => {
      const col = a.cls === 'wt' ? FRAME_COLORS.wt : a.frame === 'in-frame' ? FRAME_COLORS['in-frame'] : a.frame === 'frameshift' ? FRAME_COLORS.frameshift : '#8a8984';
      let seq = '';
      const cutsShow = res.cuts.length ? res.cuts : [res.origin];
      if (a.net != null && cs.ref) {
        const al = alleleAlignment(a, cs.ref, cutsShow, 25);
        const top = [...al.top].map((c, k) => c === '|' ? '<span class="c">|</span>' : al.marks[k] === 'i' ? '<span class="r">-</span>' : `<span class="r">${c}</span>`).join('');
        const bot = [...al.bot].map((c, k) => c === '|' ? '<span class="c">|</span>' : al.marks[k] === 'd' ? '<span class="d">-</span>' : al.marks[k] === 'i' ? `<span class="i">${c}</span>` : c).join('');
        seq = `<div class="seqcell">${top}\n${bot}</div>`;
      } else seq = `<span class="muted">${a.net == null ? 'not aligned through the cut' : 'no reference sequence loaded'}</span>`;
      const evs = a.events.length ? a.events.map((e, k) => `<span class="evchip" title="${esc(evLabel(e, chr))}">${esc(e.kind === 'D' ? `−${e.len >= 1000 ? (e.len / 1000).toFixed(1) + ' kb' : e.len + ' bp'}${e.split ? ' (split)' : ''}` : `+${e.len}${e.seq && e.len <= 6 ? ' ' + e.seq : ''}`)} <span class="muted">@${e.rel > 0 ? '+' : ''}${e.rel}</span><button class="ghost mask" data-a="${i}" data-e="${k}" title="Hide this indel everywhere (polymorphism / not of interest)">⊘</button></span>`).join(' ') : `<span>${esc(a.label)}</span>`;
      const drop = a.dropout ? ' <span class="pill" style="background:#fbe0df">dropout</span>' : '';
      return `<tr><td class="num">${i + 1}</td><td class="num"><span class="minibar" style="width:${Math.max(2, 46 * a.pct / max)}px;background:${col}"></span><b>${a.pct.toFixed(2)}</b></td><td class="num">${fmt(a.count)}${a.merged ? ` <span class="muted" title="${a.merged} reads with ±${$('#cMerge').value} bp breakpoint jitter merged in">+m</span>` : ''}</td>
        <td>${evs}${drop}</td><td><span class="pill"><i style="background:${col}"></i>${a.frame === 'wt' ? 'unmodified' : a.frame}</span></td><td>${seq}</td></tr>`;
    }).join('') + (res.alleleList.length > v.top ? `<tr><td colspan="6" class="muted">${res.alleleList.length - v.top} rarer alleles not shown (in the TSV export)</td></tr>` : '');
    $$('#cAlleles button.mask').forEach(b => b.onclick = () => {
      const e = rows[+b.dataset.a].events[+b.dataset.e];
      cs.mask.set(eventKey(e), evLabel(e, chr));
      reclassify();
      toast(`Hidden: ${evLabel(e, chr)}. Reads carrying only this change now count as unmodified.`);
    });
  }
  window.addEventListener('resize', () => { if ($('#pane-crispr').classList.contains('on') && cs.sample) draw(); });
  ctx.onTab = name => { if (name === 'crispr') { renderBams(); if (cs.sample) draw(); } };
}
