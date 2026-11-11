import igv from 'igv/dist/igv.esm.js';
import { BamFile } from '@gmod/bam';
import { BlobFile } from 'generic-filehandle2';
import { unzip } from '@gmod/bgzf-filehandle';
import { chromAlias, analyzeRegion, openSampleHeader, callsToVCF, callsToBED, callsToBEDPE, callsToTSV, depthToBedGraph, junctionsToBED } from './analysis.js';
import { buildFai } from './faidx.js';
import { initCrispr } from './crisprUI.js';
import { parseSizeFilter } from './crispr.js';
import { computeOverview, overviewToBedGraph } from './overview.js';
import { InMemoryBam } from './memBam.js';
import { openFigureEditor } from './figureEditor.js';
import { openAnnotator, getAnnotationStore, setAnnotationStore } from './annotate.js';
import { initSession } from './session.js';
import { poolSpecs } from './spectrum.js';
import { depthChart, evidenceChart, TYPE_COLORS, SERIES } from './charts.js';

const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => [...el.querySelectorAll(s)];
const fmt = n => Math.round(n).toLocaleString('en-US');
const fmtLen = n => { n = Math.abs(n); return n >= 1e6 ? (n / 1e6).toFixed(2) + ' Mb' : n >= 1e3 ? (n / 1e3).toFixed(n >= 1e4 ? 0 : 1) + ' kb' : n + ' bp'; };
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

const state = {
  browser: null, genome: null, tracks: [], nextId: 1, result: null, selected: null,
  sort: { key: 'start', dir: 1 }, typeFilter: new Set(), abort: null, colorIdx: 0, hidden: new Set(), maskIntervals: [],
};

// ---------------- toast ----------------
let toastTimer;
function toast(msg, err = false, ms = 4500) {
  const t = $('#toast'); t.textContent = msg; t.className = 'toast' + (err ? ' err' : ''); t.style.display = 'block';
  clearTimeout(toastTimer); toastTimer = setTimeout(() => t.style.display = 'none', ms);
}

// ---------------- file classification ----------------
const IDX_EXT = ['.bai', '.csi', '.crai', '.tbi', '.idx', '.fai', '.gzi'];
const lower = f => f.name.toLowerCase();
const isIndex = f => IDX_EXT.some(e => lower(f).endsWith(e));
const isFasta = f => /\.(fa|fasta|fna|fas|faa|seq)(\.gz|\.bgz)?$/.test(lower(f));
const is2bit = f => lower(f).endsWith('.2bit');
const isSizes = f => /(chrom\.sizes|\.sizes|\.genome)(\.txt)?$/.test(lower(f));
const isCyto = f => /cytoband/i.test(f.name);

function formatOf(f) {
  const n = lower(f).replace(/\.(gz|bgz)$/, '');
  const ext = n.split('.').pop();
  const map = { bam: 'bam', cram: 'cram', sam: 'sam', bw: 'bigwig', bigwig: 'bigwig', bb: 'bigbed', bigbed: 'bigbed', bed: 'bed', narrowpeak: 'narrowpeak', broadpeak: 'broadpeak',
    bedgraph: 'bedgraph', bg: 'bedgraph', wig: 'wig', gff: 'gff3', gff3: 'gff3', gtf: 'gtf', vcf: 'vcf', bedpe: 'bedpe', seg: 'seg', interact: 'interact', maf: 'maf', gwas: 'gwas', psl: 'psl', bigmaf: 'bigbed', tdf: 'tdf', txt: 'bed' };
  return map[ext] || ext;
}

function findIndex(f, files) {
  const n = f.name;
  const cands = [n + '.bai', n.replace(/\.bam$/i, '.bai'), n + '.csi', n + '.crai', n.replace(/\.cram$/i, '.crai'), n + '.tbi', n + '.idx'];
  for (const c of cands) { if (c === n) continue; const x = files.find(g => g.name === c); if (x) return x; }
  return null;
}

// ---------------- loading ----------------
async function handleFiles(fileList) {
  const files = [...fileList];
  if (!files.length) return;
  const genomeFiles = files.filter(f => isFasta(f) || is2bit(f) || (isSizes(f) && !files.some(isFasta)));
  try {
    if (genomeFiles.length) await loadGenomeFromFiles(files);
    const trackFiles = files.filter(f => !isIndex(f) && !isFasta(f) && !is2bit(f) && !isSizes(f) && !isCyto(f));
    // index files without their data file in this drop: attach to already-loaded tracks
    const orphanIdx = files.filter(f => isIndex(f) && !/\.(fai|gzi)$/i.test(f.name) && !trackFiles.some(t => findIndex(t, [f])));
    for (const f of orphanIdx) if (!(await attachLateIndex(f)) && !trackFiles.length) toast(`No loaded file matches the index ${f.name}.`, true);
    if (!trackFiles.length) return;
    if (!state.browser) {
      const bam = trackFiles.find(f => formatOf(f) === 'bam');
      if (!bam) { toast('Load a genome first (FASTA, .2bit, chrom.sizes) — or include a BAM to derive one from its header.', true, 7000); return; }
      await genomeFromBam(bam, findIndex(bam, files));
    }
    const failed = [];
    const t0 = performance.now();
    toast(`Loading ${trackFiles.length} file${trackFiles.length > 1 ? 's' : ''}…`, false, 30000);
    // load in parallel, keep the drop order in the track list
    const results = await Promise.allSettled(trackFiles.map((f, i) => addFileTrack(f, findIndex(f, files), { order: i, deferList: true })));
    results.forEach((r, i) => { if (r.status === 'rejected') { console.error(r.reason); failed.push(trackFiles[i].name + ': ' + (r.reason?.message || r.reason)); } });
    state.tracks.sort((a, b) => a.id - b.id);
    renderTrackList();
    if (!failed.length) toast(`Loaded ${trackFiles.length} file${trackFiles.length > 1 ? 's' : ''} in ${((performance.now() - t0) / 1000).toFixed(1)} s`, false, 2500);
    if (failed.length) toast('Could not load ' + failed.join(' · '), true, 10000);
  } catch (e) { console.error(e); toast(e.message || String(e), true, 8000); }
}

const fmeta = f => f ? { name: f.name, size: f.size, lastModified: f.lastModified } : null;
async function loadGenomeFromFiles(files) {
  const fa = files.find(isFasta), tb = files.find(is2bit), sizes = files.find(isSizes);
  const cyto = files.find(isCyto);
  let ref;
  if (fa) {
    let fasta = fa, fai = files.find(f => f.name === fa.name + '.fai');
    const gzi = files.find(f => f.name === fa.name + '.gzi');
    const gz = /\.(gz|bgz)$/i.test(fa.name);
    if (gz && !(fai && gzi)) {
      toast('Decompressing ' + fa.name + ' in memory… (bgzip + .fai + .gzi avoids this for large genomes)', false, 10000);
      let blob;
      try { blob = new Blob([await unzip(new Uint8Array(await fa.arrayBuffer()))]); }
      catch (e) { const ds = fa.stream().pipeThrough(new DecompressionStream('gzip')); blob = await new Response(ds).blob(); }
      fasta = new File([blob], fa.name.replace(/\.(gz|bgz)$/i, ''));
      fai = null;
    }
    if (!fai) {
      toast('Indexing ' + fasta.name + ' …', false, 60000);
      const text = await buildFai(fasta, p => toast(`Indexing ${fasta.name} … ${(p * 100).toFixed(0)}%`, false, 60000));
      fai = new File([text], fasta.name + '.fai');
      toast('Built FASTA index for ' + fasta.name);
    }
    ref = { id: fa.name, name: fa.name.replace(/\.(fa|fasta|fna|fas)(\.gz|\.bgz)?$/i, ''), fastaURL: fasta, indexURL: fai };
    if (gz && gzi && fasta === fa) ref.compressedIndexURL = gzi;
    ref._chroms = (await fai.text()).trim().split('\n').map(l => { const [n, len] = l.split('\t'); return { name: n, len: +len }; });
  } else if (tb) {
    ref = { id: tb.name, name: tb.name.replace(/\.2bit$/i, ''), twoBitURL: tb, format: '2bit' };
  } else if (sizes) {
    ref = { id: sizes.name, name: sizes.name.replace(/\.(chrom\.)?sizes.*$/i, ''), format: 'chromsizes', fastaURL: sizes };
    ref._chroms = (await sizes.text()).trim().split('\n').map(l => { const [n, len] = l.split(/\s+/); return { name: n, len: +len }; });
  }
  if (cyto) ref.cytobandURL = cyto;
  const src = { type: 'files', files: [fa, fa && files.find(f => f.name === fa.name + '.fai'), fa && files.find(f => f.name === fa.name + '.gzi'), tb, !fa && !tb && sizes, cyto].filter(Boolean).map(fmeta) };
  await setGenome(ref);
  state.genomeSrc = src;
}

async function genomeFromBam(bamFile, bai) {
  const bam = new BamFile({ bamFilehandle: new BlobFile(bamFile), ...(bai ? (bai.name.endsWith('.csi') ? { csiFilehandle: new BlobFile(bai) } : { baiFilehandle: new BlobFile(bai) }) : {}) });
  await bam.getHeader();
  const chroms = (bam.indexToChr || []).map(c => ({ name: c.refName, len: c.length }));
  if (!chroms.length) throw new Error('Could not read reference names from ' + bamFile.name);
  const sizes = new File([chroms.map(c => `${c.name}\t${c.len}`).join('\n') + '\n'], 'from_bam.chrom.sizes');
  toast('No genome loaded: using chromosome sizes from ' + bamFile.name + ' (no reference sequence).', false, 7000);
  await setGenome({ id: 'bam-header', name: bamFile.name + ' header', format: 'chromsizes', fastaURL: sizes, _chroms: chroms });
  state.genomeSrc = { type: 'bam-header', files: [fmeta(bamFile), fmeta(bai)].filter(Boolean) };
}

async function setGenome(ref) {
  ensureKeepAlive();
  const chroms = ref._chroms; delete ref._chroms;
  const first = chroms?.[0];
  const locus = first ? `${first.name}:1-${Math.min(first.len, 50000)}` : undefined;
  if (!state.browser) {
    $('#welcome').style.display = 'none';
    state.browser = await igv.createBrowser($('#igvDiv'), {
      reference: ref, locus, loadDefaultGenomes: false, search: false, showNavigation: false, showCursorGuide: false,
      showCenterGuide: false, showTrackLabels: true, showSampleNames: true, showChromosomeWidget: true,
      tracks: [],
    });
    state.browser.on('locuschange', () => { updateLocusBox(); resortAll(); });
    $('#igvDiv').addEventListener('mouseup', () => setTimeout(syncOrderFromBrowser, 50));
    state.browser.on('trackorderchanged', () => setTimeout(syncOrderFromBrowser));
    state.browser.on('trackremoved', tracks => { for (const t of (Array.isArray(tracks) ? tracks : [tracks])) { const r = state.tracks.find(x => x.igv === t); if (r && !r.hiding) removeTrackRecord(r, false); } });
  } else {
    await state.browser.loadGenome(ref);
    if (locus) await state.browser.search(locus);
    state.tracks = [];
    renderTrackList();
  }
  state.genome = { name: ref.name, chroms: chroms || state.browser.genome.chromosomeNames.map(n => ({ name: n, len: state.browser.genome.getChromosome(n)?.bpLength })) };
  $('#genomeName').textContent = ref.name;
  updateLocusBox();
}

function nextColor() { const c = SERIES[state.colorIdx % SERIES.length]; state.colorIdx++; return c; }

async function addFileTrack(file, index, extra = {}) {
  const myId = state.nextId++;
  const format = extra.format || formatOf(file);
  const name = extra.name || file.name;
  const needsIndex = ['bam', 'cram'].includes(format) || (/\.(gz|bgz)$/i.test(file.name) && ['vcf', 'bed', 'gff3', 'gtf', 'bedgraph'].includes(format));
  const cfg = { url: file, name, format, ...extra.cfg };
  if (index) cfg.indexURL = index; else if (needsIndex) { cfg.indexed = false; }
  if (format === 'cram' && !index) throw new Error(file.name + ': CRAM needs a .crai index');
  if (format === 'bam' && !index && file.size > MAX_UNINDEXED) toast(`${file.name} has no index; it is too large to display or analyse without one (samtools index).`, true, 9000);
  const color = extra.color || (['bigwig', 'bedgraph', 'wig'].includes(format) ? nextColor() : undefined);
  if (color) cfg.color = color;
  if (format === 'bam') Object.assign(cfg, window.__extra || {}, { showSoftClips: false, viewAsPairs: false, height: 220, displayMode: (state.tracks.filter(t => t.format === 'bam').length >= 1 || extra.order > 0) ? 'SQUISHED' : 'EXPANDED' });
  cfg.order = myId * 10;
  const t = await state.browser.loadTrack(cfg);
  const isSig = ['bigwig', 'bedgraph', 'wig'].includes(format);
  const rec = { id: myId, name, format, file, index, igv: t, color: color || defaultColor(format), height: t?.height || 50, mode: isSig ? 'full' : format === 'bam' ? (cfg.displayMode === 'SQUISHED' ? 'squish' : 'pack') : 'pack', generated: !!extra.generated, cfg };
  if (format === 'bam') {
    rec.sample = { use: true, group: file.name.replace(/\.(bam|sorted|markdup|dedup)+$/ig, '').replace(/\.(bam)$/i, ''), isControl: false };
    attachSampleReader(rec, index);
  }
  state.tracks.push(rec);
  if (!extra.deferList) renderTrackList();
  return rec;
}

// open a reader for analysis: indexed BamFile, or (no index) the whole file in memory
const MAX_UNINDEXED = 400e6;
function attachSampleReader(rec, index) {
  const smp = rec.sample; const file = rec.file;
  delete smp.err; smp.bam = null; smp.inMemory = false;
  let bam;
  if (index) bam = new BamFile({ bamFilehandle: new BlobFile(file), ...(index.name.endsWith('.csi') ? { csiFilehandle: new BlobFile(index) } : { baiFilehandle: new BlobFile(index) }) });
  else if (file.size <= MAX_UNINDEXED) { bam = new InMemoryBam(file); smp.inMemory = true; rec.ovStatus = 'no index — reading into memory…'; }
  else { smp.err = `no .bai/.csi index and too large (${(file.size / 1e6).toFixed(0)} MB) to read without one — drop its index, or run samtools index`; smp.use = false; smp.ready = Promise.resolve(); return; }
  smp.ready = (async () => {
    try {
      const h = await openSampleHeader(bam);
      smp.bam = bam; Object.assign(smp, h);
      if (h.sampleName && !smp.groupEdited) smp.group = h.sampleName;
      rec.ovStatus = 'scanning…'; renderTrackList(); state.onSamplesChanged?.();
      computeOverview(bam, { fileSize: index ? file.size : 0 }).then(ov => onOverview(rec, ov)).catch(e => { console.warn(e); rec.ovStatus = 'no overview'; renderTrackList(); });
    } catch (e) { console.warn(e); smp.use = false; smp.err = 'could not read: ' + (e.message || e); renderTrackList(); state.onSamplesChanged?.(); }
  })();
}
// an index dropped after its BAM: attach it and reload the track with the index
async function attachLateIndex(idxFile) {
  const rec = state.tracks.find(r => ['bam', 'cram'].includes(r.format) && !r.index && findIndex(r.file, [idxFile]));
  if (!rec) return false;
  rec.index = idxFile;
  if (rec.cfg) { rec.cfg.indexURL = idxFile; delete rec.cfg.indexed; }
  if (rec.igv) { rec.hiding = true; try { state.browser.removeTrack(rec.igv); } finally { rec.hiding = false; } rec.igv = await state.browser.loadTrack({ ...rec.cfg, color: rec.color }); }
  if (rec.sample) attachSampleReader(rec, idxFile);
  toast(`Index ${idxFile.name} attached to ${rec.name}.`);
  renderTrackList(); state.onSamplesChanged?.();
  return true;
}

// ---------------- read-location overview ----------------
let ovTimer = null, lastJumpAt = 0;
function onOverview(rec, ov) {
  rec.overview = ov; rec.ovStatus = '';
  if (!state.tracks.includes(rec)) return;
  // let small BAMs draw reads when zoomed out (scaled so ~20k reads are fetched at most)
  if (ov.scan && rec.igv) {
    let vw = Infinity;
    for (const d of Object.values(ov.perChrom)) vw = Math.min(vw, d.len * Math.min(1, 20000 / Math.max(1, d.reads)));
    vw = Math.max(30000, Math.round(isFinite(vw) ? vw : 30000));
    rec.igv.visibilityWindow = vw; if (rec.cfg) rec.cfg.visibilityWindow = vw;
    rec.igv.trackView?.updateViews?.();
  }
  renderTrackList();
  clearTimeout(ovTimer); ovTimer = setTimeout(() => { rebuildOverviewTrack(); maybeJumpToReads(rec); }, 150);
}
function overviewRecs() { return state.tracks.filter(t => t.overview && Object.keys(t.overview.perChrom).length); }
async function rebuildOverviewTrack() {
  const recs = overviewRecs();
  const old = state.tracks.find(t => t.isOverview);
  if (old) removeTrackRecord(old, true);
  if (!recs.length || !state.browser) return;
  const tracks = recs.map((r, i) => ({ url: new File([overviewToBedGraph(r.overview, r.name)], r.name + '.overview.bedgraph'), format: 'bedgraph', name: r.name, color: SERIES[i % SERIES.length] }));
  const cfg = { type: 'merged', name: 'Read density · all BAMs (any zoom)', tracks, height: 60, autoscale: true, order: 1 };
  try {
    const t = await state.browser.loadTrack(cfg);
    state.tracks.unshift({ id: 0, name: cfg.name, format: 'merged', igv: t, color: '#52514e', height: 60, mode: 'full', generated: true, isOverview: true, cfg });
    renderTrackList();
  } catch (e) { console.warn(e); }
}
function viewHasReads() {
  const l = currentLocus(); if (!l) return true;
  return overviewRecs().some(r => r.overview.hotspots.some(h => (h.chr === l.chr || chromAlias(l.chr, new Set([h.chr])) === h.chr) && h.start < l.end && h.end > l.start));
}
function hotspotLocus(h) {
  const chr = state.browser.genome.getChromosomeName?.(h.chr) || h.chr;
  const pad = Math.max(50, Math.round((h.end - h.start) * 0.1));
  return `${chr}:${Math.max(1, h.start - pad + 1)}-${h.end + pad}`;
}
function maybeJumpToReads(rec) {
  if (Date.now() - lastJumpAt < 3000) return;
  if (viewHasReads()) return;
  const h = rec.overview.hotspots[0]; if (!h) return;
  lastJumpAt = Date.now();
  goto(hotspotLocus(h));
  toast(`No reads in view — jumped to ${rec.name}'s reads at ${h.chr}:${fmt(h.start + 1)}-${fmt(h.end)}. All read regions are listed under the track.`, false, 7000);
}

function defaultColor(fmt) { return { bam: '#8a8984', cram: '#8a8984', vcf: '#2a78d6', gff3: '#184f95', gtf: '#184f95', bed: '#184f95', bigbed: '#184f95' }[fmt] || '#184f95'; }

// ---------------- track list UI (UCSC-style controls) ----------------
const MODES = { hide: null, dense: 'COLLAPSED', squish: 'SQUISHED', pack: 'EXPANDED', full: 'FULL' };
function renderTrackList() {
  const el = $('#trackList');
  if (!state.tracks.length) { el.innerHTML = '<div class="empty">No tracks loaded.</div>'; $('#overlaySec').style.display = 'none'; return; }
  el.innerHTML = '';
  for (const r of state.tracks) {
    const d = document.createElement('div'); d.className = 'trow';
    const isSignal = ['bigwig', 'bedgraph', 'wig', 'merged'].includes(r.format);
    const modeOpts = (r.format === 'bam' ? ['hide', 'squish', 'pack', 'full'] : isSignal ? ['hide', 'full'] : ['hide', 'dense', 'squish', 'pack'])
      .map(m => `<option ${m === r.mode ? 'selected' : ''}>${m}</option>`).join('');
    d.draggable = true; d.dataset.tid = r.id;
    d.innerHTML = `<div class="r1"><span class="grip" title="Drag to reorder">⠿</span><span class="kind">${esc(r.format)}</span><span class="nm" title="${esc(r.name)}${r.file && r.file.name !== r.name ? '\nfile: ' + esc(r.file.name) : ''}\nclick to rename">${esc(r.name)}</span>
      <button class="ghost ren" title="Rename (display name used in the browser, analyses, figures and exports)">✎</button><button class="ghost mv" data-mv="-1" title="Move up">▲</button><button class="ghost mv" data-mv="1" title="Move down">▼</button><button class="ghost x" title="Remove track">×</button></div>
      <div class="r2"><select class="mode" title="Visibility (UCSC style)">${modeOpts}</select>
      <input type="color" class="col" value="${r.color}" title="Colour">
      <label title="Panel height (px)">h <input type="number" class="hgt" value="${Math.round(r.igv?.height || r.height || 50)}" min="15" max="3000" step="10"></label><button class="ghost hstep" data-f="0.66" title="Smaller">−</button><button class="ghost hstep" data-f="1.5" title="Larger">+</button>
      ${r.format === 'bam' ? `<label title="Height of the coverage panel inside this BAM track">cov h <input type="number" class="covh" value="${Math.round(r.igv?.coverageTrackHeight ?? 50)}" min="10" max="400" step="5" style="width:48px"></label>` : ''}
      <button class="ghost cen" title="Zoom to this track's data" style="display:none">⤢</button></div>`;
    if (['bed', 'gff3', 'gtf', 'bigbed', 'narrowpeak', 'broadpeak', 'gwas', 'psl', 'interact', 'bedpe', 'vcf'].includes(r.format) || (r.igv && r.igv.type === 'annotation')) {
      const t = r.igv || {};
      d.insertAdjacentHTML('beforeend', `<div class="sv"><label title="Gene / feature names under the features"><input type="checkbox" class="pnames" ${t.gvHideNames ?? r.hideNames ? '' : 'checked'}> names</label><label title="Strand arrows on genes / features"><input type="checkbox" class="parrows" ${t.gvHideArrows ?? r.hideArrows ? '' : 'checked'}> strand arrows</label>
        <label title="Row height / spacing; the panel shrinks to fit its content">rows <select class="pcompact">${['normal', 'compact', 'tight'].map(v => `<option ${(r.compact || 'normal') === v ? 'selected' : ''}>${v}</option>`).join('')}</select></label></div>`);
    }
    if (r.format === 'bam') {
      const ov = r.overview;
      let html = '';
      if (r.ovStatus) html = `<span class="muted">reads: ${esc(r.ovStatus)}</span>`;
      else if (ov) {
        const hs = ov.hotspots;
        html = `<details class="hot"><summary>${ov.scan ? fmt(ov.total) + ' reads' : '~' + fmt(ov.total) + ' reads (index estimate)'} in <b>${hs.length}</b> region${hs.length === 1 ? '' : 's'}</summary>` +
          hs.slice(0, 200).map((h, i) => `<a href="#" data-hot="${i}" class="loc">${esc(h.chr)}:${fmt(h.start + 1)}-${fmt(h.end)}</a> <span class="muted">${ov.scan ? 'peak ' + fmt(h.peak) + '×' : ''}</span>`).join('<br>') + (hs.length > 200 ? `<br><span class="muted">+${hs.length - 200} more</span>` : '') + '</details>';
      }
      const t = r.igv;
      const showCov = t ? t.showCoverage !== false : r.cfg?.showCoverage !== false, showAl = t ? t.showAlignments !== false : r.cfg?.showAlignments !== false;
      const sorts = [['', 'default order'], ['GV_DEL_MAX', 'largest deletion first'], ['GV_DEL_MIN', 'smallest deletion first'], ['GV_INDEL_MAX', 'largest indel first'], ['GV_RANDOM', 'random'], ['BASE', 'base at centre'], ['strand', 'strand'], ['START', 'start position'], ['MQ', 'mapping quality'], ['INSERT_SIZE', 'insert size']];
      d.insertAdjacentHTML('beforeend', `<div class="sv"><label title="Coverage panel above the reads"><input type="checkbox" class="pcov" ${showCov ? 'checked' : ''}> coverage</label><label title="Individual reads"><input type="checkbox" class="paln" ${showAl ? 'checked' : ''}> reads</label>
        <label title="Order the read rows; deletion/indel/random sorts look at the whole view and follow you as you move">sort <select class="psort">${sorts.map(([v, l]) => `<option value="${v}" ${(r.sortMode || '') === v ? 'selected' : ''}>${l}</option>`).join('')}</select></label>${r.sortMode === 'GV_RANDOM' ? '<button class="ghost preshuf" title="Shuffle again">↻</button>' : ''}${html ? '<br>' + html : ''}</div>`);
    }
    if (r.sample) {
      const s = r.sample;
      d.insertAdjacentHTML('beforeend', `<div class="sv">${s.err ? `<span class="muted">${esc(s.err)}</span>` : `
        ${s.inMemory ? '<span class="muted" title="No .bai/.csi: the whole BAM is read into memory. Drop the index any time to switch.">no index · in memory</span>' : ''}
        <label title="Include in variant / SV analysis"><input type="checkbox" class="use" ${s.use ? 'checked' : ''}> analyse</label>
        <label title="BAMs with the same group are merged (reads pooled)">group <input type="text" class="grp" value="${esc(s.group)}"></label>
        <label title="Reference sample for log2 ratio and depth CNV"><input type="radio" name="ctrl" class="ctl" ${s.isControl ? 'checked' : ''}> control</label>`}</div>`);
    }
    $('.x', d).onclick = () => { removeTrackRecord(r, true); if (r.overview) rebuildOverviewTrack(); };
    $$('.mv', d).forEach(b => b.onclick = () => moveTrack(r, +b.dataset.mv));
    $$('.hstep', d).forEach(b => b.onclick = () => { const h = Math.max(15, Math.round((r.igv?.height || r.height || 50) * +b.dataset.f)); setHeight(r, h); $('.hgt', d).value = h; });
    if ($('.covh', d)) $('.covh', d).onchange = e => setCoverageHeight(r, +e.target.value);
    d.addEventListener('dragstart', e => { if (e.target.closest('input,select,button')) { e.preventDefault(); return; } dragTid = r.id; d.classList.add('dragging'); e.dataTransfer.effectAllowed = 'move'; e.dataTransfer.setData('text/plain', 'track'); });
    d.addEventListener('dragend', () => { d.classList.remove('dragging'); $$('.trow').forEach(x => x.classList.remove('dropbefore', 'dropafter')); });
    d.addEventListener('dragover', e => { if (dragTid == null) return; e.preventDefault(); e.stopPropagation(); const after = e.offsetY > d.clientHeight / 2; d.classList.toggle('dropafter', after); d.classList.toggle('dropbefore', !after); });
    d.addEventListener('dragleave', () => d.classList.remove('dropbefore', 'dropafter'));
    d.addEventListener('drop', e => { if (dragTid == null) return; e.preventDefault(); e.stopPropagation(); const after = e.offsetY > d.clientHeight / 2; dropTrack(dragTid, r.id, after); dragTid = null; });
    const startRename = () => {
      const nm = $('.nm', d); const inp = document.createElement('input'); inp.type = 'text'; inp.value = r.name; inp.className = 'nmedit';
      nm.replaceWith(inp); inp.focus(); inp.select();
      let done = false;
      const finish = ok => { if (done) return; done = true; const v = inp.value.trim(); if (ok && v && v !== r.name) renameTrack(r, v); else renderTrackList(); };
      inp.onkeydown = e => { if (e.key === 'Enter') finish(true); if (e.key === 'Escape') finish(false); };
      inp.onblur = () => finish(true);
    };
    $('.nm', d).onclick = startRename; $('.ren', d).onclick = startRename;
    if ($('.pcompact', d)) $('.pcompact', d).onchange = e => { setCompact(r, e.target.value); renderTrackList(); };
    if ($('.pnames', d)) { $('.pnames', d).onchange = e => setFeatureOpt(r, 'gvHideNames', !e.target.checked); $('.parrows', d).onchange = e => setFeatureOpt(r, 'gvHideArrows', !e.target.checked); }
    if ($('.psort', d)) { $('.psort', d).onchange = e => { r.sortMode = e.target.value; r.shuffleSeed = (r.shuffleSeed || 0) + 1; sortReads(r); renderTrackList(); }; }
    if ($('.preshuf', d)) $('.preshuf', d).onclick = () => { r.shuffleSeed = (r.shuffleSeed || 0) + 1; sortReads(r, true); };
    if ($('.pcov', d)) { $('.pcov', d).onchange = e => setBamPart(r, 'showCoverage', e.target.checked); $('.paln', d).onchange = e => setBamPart(r, 'showAlignments', e.target.checked); }
    $$('[data-hot]', d).forEach(a => a.onclick = e => { e.preventDefault(); goto(hotspotLocus(r.overview.hotspots[+a.dataset.hot])); });
    $('.mode', d).onchange = e => setMode(r, e.target.value);
    $('.col', d).oninput = e => setColor(r, e.target.value);
    $('.hgt', d).onchange = e => setHeight(r, +e.target.value);
    if (r.sample && !r.sample.err) {
      $('.use', d).onchange = e => { r.sample.use = e.target.checked; };
      $('.grp', d).onchange = e => { r.sample.group = e.target.value.trim() || r.name; };
      $('.ctl', d).onchange = () => { for (const t of state.tracks) if (t.sample) t.sample.isControl = (t === r); };
      $('.ctl', d).ondblclick = e => { e.target.checked = false; r.sample.isControl = false; };
    }
    el.appendChild(d);
  }
  // overlay candidates
  const sig = state.tracks.filter(r => ['bigwig', 'bedgraph', 'wig'].includes(r.format));
  $('#overlaySec').style.display = sig.length >= 2 ? '' : 'none';
  state.onSamplesChanged?.();
  $('#overlayList').innerHTML = sig.map(r => `<label class="status"><input type="checkbox" value="${r.id}"> <span class="swatch" style="background:${r.color}"></span> ${esc(r.name)}</label>`).join('');
}

function trackView(r) { return r.igv?.trackView || state.browser.trackViews.find(tv => tv.track === r.igv); }

async function setMode(r, mode) {
  const prev = r.mode; r.mode = mode;
  if (mode === 'hide') {
    if (r.igv) { r.hiding = true; try { state.browser.removeTrack(r.igv); } finally { r.hiding = false; } r.igv = null; }
    return;
  }
  if (!r.igv) {               // re-show: reload with current settings, same order
    const cfg = { ...r.cfg, color: r.color, height: r.height };
    if (r.cfg?.type !== 'merged') cfg.displayMode = MODES[mode] === 'FULL' && !r.cfg?.format?.match(/bam|cram/) ? 'EXPANDED' : (MODES[mode] || 'EXPANDED');
    r.igv = await state.browser.loadTrack(cfg);
    if (r.hideNames) r.igv.gvHideNames = true; if (r.hideArrows) r.igv.gvHideArrows = true;
    if (prev === 'hide' && r.format !== 'bam') return;
  }
  const tv = trackView(r); if (!tv) return;
  const m = MODES[mode];
  const t = r.igv;
  if (t.alignmentTrack) {
    t.alignmentTrack.setDisplayMode ? t.alignmentTrack.setDisplayMode(m === 'COLLAPSED' ? 'SQUISHED' : m) : (t.alignmentTrack.displayMode = m);
  } else if ('displayMode' in t && m !== 'FULL') { t.displayMode = m; }
  tv.checkContentHeight?.(); tv.repaintViews?.(); tv.updateViews?.();
}

function renameTrack(r, name) {
  const old = r.name;
  r.name = name;
  if (r.cfg) r.cfg.name = name;
  if (r.igv) { try { r.igv.name = name; } catch (e) { console.warn(e); } r.igv.trackView?.repaintViews?.(); }
  if (r.sample && !r.sample.groupEdited && (r.sample.group === old || r.sample.group === old.replace(/\.bam$/i, '') || r.sample.group === r.sample.sampleName)) r.sample.group = name.replace(/\.bam$/i, '');
  renderTrackList();
  state.onTrackRenamed?.(r, old);
  if (r.overview) rebuildOverviewTrack();
  toast(`Renamed “${old}” → “${name}”.`, false, 2500);
}

// show / hide a whole igv track view (ideogram, ruler…) on screen and in exports
function setTrackViewVisible(tv, on) {
  if (!tv) return;
  tv.track.gvHidden = !on;
  const els = [tv.axis, tv.outerScroll, tv.dragHandle, tv.gearContainer, tv.sampleNameViewport?.viewport, tv.sampleInfoViewport?.viewport, ...(tv.viewports || []).map(v => v.viewportElement)];
  for (const el of els) if (el) el.style.display = on ? '' : 'none';
}
function setIdeogram(on) { if (!state.browser) return; setTrackViewVisible(state.browser.trackViews.find(t => t.track.type === 'ideogram' || t.track.id === 'ideogram'), on); }

// compact annotation tracks: smaller rows / margins and height fitted to the content
const COMPACT = { normal: { e: 30, s: 15, m: 10 }, compact: { e: 18, s: 8, m: 4 }, tight: { e: 12, s: 5, m: 2 } };
function setCompact(r, level) {
  const t = r.igv; if (!t || !('expandedRowHeight' in t)) return;
  const c = COMPACT[level] || COMPACT.normal;
  r.compact = level;
  t.expandedRowHeight = c.e; t.squishedRowHeight = c.s; t.margin = c.m;
  t.minHeight = 0; t.autoHeight = true;
  if (r.cfg) Object.assign(r.cfg, { expandedRowHeight: c.e, squishedRowHeight: c.s, margin: c.m });
  try {
    const feats = t.trackView.viewports.flatMap(v => v.cachedFeatures?.features || (Array.isArray(v.cachedFeatures) ? v.cachedFeatures : []));
    const h = Math.max(c.e + c.m, t.computePixelHeight ? t.computePixelHeight(feats) : c.e + c.m);
    t.trackView.setTrackHeight(h, true); r.height = h;
  } catch (e) { t.trackView.checkContentHeight?.(); }
  t.trackView.repaintViews?.();
}
function fitBam(r, cap) {
  const t = r.igv; if (!t?.alignmentTrack) return;
  if (cap < 200) t.alignmentTrack.setDisplayMode?.('SQUISHED');
  let h = 0;
  for (const vp of t.trackView.viewports) { try { h = Math.max(h, t.computePixelHeight(vp.cachedFeatures)); } catch {} }
  if (!h) h = (t.showCoverage ? t.coverageTrackHeight : 0) + 10;
  setHeight(r, Math.max(30, Math.min(cap, Math.ceil(h) + 4)));
  r.mode = cap < 200 ? 'squish' : r.mode;
}
function compactAll(level) {
  for (const r of state.tracks) {
    if (r.format === 'bam') fitBam(r, level === 'tight' ? 100 : level === 'compact' ? 160 : 220);
    if (r.igv && 'expandedRowHeight' in r.igv && r.format !== 'bam') setCompact(r, level);
    if (level !== 'normal' && ['bigwig', 'bedgraph', 'wig', 'merged'].includes(r.format) && (r.igv?.height || 0) > 45) setHeight(r, 40);
  }
  renderTrackList();
}
$('#compactSel').onchange = e => compactAll(e.target.value);

function setFeatureOpt(r, prop, hide) {
  r[prop === 'gvHideNames' ? 'hideNames' : 'hideArrows'] = hide;
  if (r.igv) { r.igv[prop] = hide; r.igv.trackView?.repaintViews(); }
}
// hide names / arrows on every annotation track at once
function setAllFeatureOpt(prop, hide) {
  for (const r of state.tracks) if (r.igv && r.format !== 'bam' && !['bigwig', 'bedgraph', 'wig', 'merged'].includes(r.format)) setFeatureOpt(r, prop, hide);
  renderTrackList();
}

// sort read rows in every visible panel of a BAM track
function sortReads(r, reshuffle = false) {
  const t = r.igv; const at = t?.alignmentTrack; if (!at) return;
  const mode = r.sortMode;
  for (const vp of t.trackView.viewports) {
    const ac = vp.cachedFeatures; if (!ac || !ac.sortRows) continue;
    const rf = vp.referenceFrame;
    const start = Math.floor(rf.start), end = Math.ceil(rf.start + (vp.viewportElement?.clientWidth || 1000) * rf.bpPerPixel);
    const pos = Math.floor((start + end) / 2);
    if (!mode) { at.sortObject = undefined; continue; }
    if (mode === 'GV_RANDOM' && reshuffle) for (const g of ac.packedGroups?.values?.() || []) for (const row of g.rows || []) delete row.__rnd;
    const opt = { chr: rf.chr, position: pos, option: mode, direction: mode === 'BASE' || mode === 'strand' ? true : true, sortAsPairs: at.viewAsPairs, start, end };
    at.sortObject = opt;
    try { ac.sortRows(opt); } catch (e) { console.warn(e); }
    vp.repaint();
  }
}
let resortTimer = null;
function resortAll() { clearTimeout(resortTimer); resortTimer = setTimeout(() => { for (const r of state.tracks) if (r.sortMode) sortReads(r); }, 350); }

// ---------------- panel order & size ----------------
let dragTid = null;
function applyOrder() {
  state.tracks.forEach((t, i) => { const o = (i + 1) * 10; if (t.igv) t.igv.order = o; if (t.cfg) t.cfg.order = o; });
  try { state.browser.reorderTracks(); } catch (e) { console.warn(e); }
  renderTrackList();
}
function moveTrack(r, dir) {
  const i = state.tracks.indexOf(r), j = i + dir;
  if (i < 0 || j < 0 || j >= state.tracks.length) return;
  [state.tracks[i], state.tracks[j]] = [state.tracks[j], state.tracks[i]];
  applyOrder();
}
function dropTrack(fromId, toId, after) {
  if (fromId === toId) return;
  const from = state.tracks.find(t => t.id == fromId); if (!from) return;
  state.tracks.splice(state.tracks.indexOf(from), 1);
  let k = state.tracks.findIndex(t => t.id == toId); if (after) k++;
  state.tracks.splice(k, 0, from);
  applyOrder();
}
// keep the list in step when tracks are dragged inside the browser itself
function syncOrderFromBrowser() {
  const order = state.browser.trackViews.map(tv => tv.track);
  const idx = t => { const k = order.indexOf(t.igv); return k < 0 ? 1e9 : k; };
  const before = state.tracks.map(t => t.id).join();
  state.tracks.sort((a, b) => idx(a) - idx(b));
  if (state.tracks.map(t => t.id).join() !== before) renderTrackList();
}
function setCoverageHeight(r, h) {
  const t = r.igv; if (!t || !t.alignmentTrack) return;
  const old = t.coverageTrackHeight; t.coverageTrackHeight = Math.max(10, h); if (t.coverageTrack) t.coverageTrack.height = t.coverageTrackHeight;
  if (r.cfg) r.cfg.coverageTrackHeight = t.coverageTrackHeight;
  t.alignmentTrack.setTop?.(t.coverageTrack, t.showCoverage);
  t.trackView.setTrackHeight(Math.max(30, (t.height || 0) + (t.coverageTrackHeight - old)), true);
  t.trackView.checkContentHeight?.(); t.trackView.repaintViews();
}

function setBamPart(r, prop, on) {
  if (r.cfg) r.cfg[prop] = on;
  const t = r.igv; if (!t) return;
  t[prop] = on;
  if (!t.autoHeight) { const h = (t.showCoverage ? t.coverageTrackHeight : 0) + (t.showAlignments ? t.alignmentTrack.height : 0); t.trackView.setTrackHeight(Math.max(20, h)); }
  t.trackView.checkContentHeight?.(); t.trackView.repaintViews();
}

function setColor(r, c) {
  r.color = c; const t = r.igv; if (!t) return;
  if (t.alignmentTrack) { t.coverageTrack && (t.coverageTrack.color = c); t.alignmentTrack.color = c; }
  else if (typeof t.setColor === 'function' && t.type !== 'merged') { try { t.color = c; } catch { t.setColor(c); } }
  else t.color = c;
  trackView(r)?.repaintViews();
  renderOverlayChips();
}
function renderOverlayChips() { $$('#overlayList label').forEach(l => { const r = state.tracks.find(t => t.id == $('input', l).value); if (r) $('.swatch', l).style.background = r.color; }); }

function setHeight(r, h) { r.height = h; if (r.cfg) r.cfg.height = h; if (!r.igv) return; r.igv.autoHeight = false; if (r.igv.alignmentTrack) { r.igv.height = h; } const tv = trackView(r); if (tv) { tv.setTrackHeight(h, true); tv.repaintViews?.(); } }

function removeTrackRecord(r, fromUI) {
  state.tracks = state.tracks.filter(t => t !== r);
  if (fromUI && r.igv) try { state.browser.removeTrack(r.igv); } catch (e) { console.warn(e); }
  renderTrackList();
}

$('#btnOverlay').onclick = async () => {
  const ids = $$('#overlayList input:checked').map(i => +i.value);
  const recs = state.tracks.filter(t => ids.includes(t.id));
  if (recs.length < 2) return toast('Tick at least two signal tracks.', true);
  const tracks = recs.map(r => ({ url: r.file, format: r.format, name: r.name, color: r.color, ...(r.index ? { indexURL: r.index } : {}) }));
  const cfg = { type: 'merged', name: recs.map(r => r.name).join(' + '), tracks, height: 120, autoscale: true, order: state.nextId * 10 };
  const t = await state.browser.loadTrack(cfg);
  state.tracks.push({ id: state.nextId++, name: 'Overlay: ' + recs.map(r => r.name).join(' + '), format: 'merged', igv: t, color: '#52514e', height: 120, mode: 'full', generated: true, cfg });
  renderTrackList();
};

// ---------------- navigation ----------------
function currentLocus() {
  const rf = state.browser?.referenceFrameList?.[0];
  if (!rf) return null;
  const start = Math.max(0, Math.floor(rf.start));
  let end = Math.floor(rf.end ?? (rf.start + rf.bpPerPixel * (state.browser.calculateViewportWidth?.(1) || 1000)));
  const len = chrLength(rf.chr); if (len) end = Math.min(end, len);
  return { chr: rf.chr, start, end };
}
let locusDirty = false;
$('#locus').addEventListener('input', () => { locusDirty = true; });
$('#locus').addEventListener('blur', () => { setTimeout(() => { if (document.activeElement !== $('#btnGo')) { locusDirty = false; updateLocusBox(); } }, 200); });
function updateLocusBox(force) {
  const l = currentLocus(); if (!l) return;
  const len = chrLength(l.chr);
  const atWhole = len && (l.end - l.start) >= len * 0.999;
  $$('[data-zoom]').forEach(b => { const z = b.dataset.zoom; if (z === 'chr' || (+z > 0 && +z < 1)) b.disabled = !!atWhole; });
  if (locusDirty && !force) return;          // never overwrite what the user is typing
  $('#locus').value = `${l.chr}:${fmt(l.start + 1)}-${fmt(l.end)}`;
  $('#locusSize').textContent = fmtLen(l.end - l.start);
}
// keep every view inside its chromosome: never wider than the chromosome, never past its ends
function chrLength(chr) { try { return state.browser.genome.getChromosome(chr)?.bpLength || null; } catch { return null; } }
function clampRange(chr, s, e) {
  const len = chrLength(chr); if (!len) return { chr, s: Math.max(0, s), e };
  let w = Math.max(1, Math.min(e - s, len));
  s = Math.max(0, Math.min(Math.round(s), len - w));
  return { chr, s, e: s + w };
}
function clampLocusString(str) {
  const parts = str.trim().split(/\s+/);
  return parts.map(p => {
    const m = /^([^:]+):([\d,]+)-([\d,]+)$/.exec(p);
    if (!m) return p;
    const chr = state.browser.genome.getChromosomeName?.(m[1]) || m[1];
    const r = clampRange(chr, +m[2].replace(/,/g, '') - 1, +m[3].replace(/,/g, ''));
    return `${chr}:${r.s + 1}-${r.e}`;
  }).join(' ');
}
async function goto(str) {
  if (!state.browser) return toast('Load a genome first.', true);
  locusDirty = false;
  try { str = clampLocusString(str); } catch {}
  try { const ok = await state.browser.search(str.trim()); if (ok === false) toast('Not found: ' + str, true); } catch (e) { toast('Not found: ' + str, true); }
}
$('#btnGo').onclick = () => goto($('#locus').value);
$('#locus').addEventListener('keydown', e => { if (e.key === 'Enter') goto(e.target.value); });
$$('[data-move]').forEach(b => b.onclick = () => {
  const l = currentLocus(); if (!l) return;
  const w = l.end - l.start, d = Math.round(w * +b.dataset.move);
  const r = clampRange(l.chr, l.start + d, l.start + d + w);
  goto(`${l.chr}:${r.s + 1}-${r.e}`);
});
$$('[data-zoom]').forEach(b => b.onclick = () => {
  const l = currentLocus(); if (!l) return;
  const z = b.dataset.zoom;
  if (z === 'chr') return goto(l.chr);
  const len = chrLength(l.chr) || Infinity;
  const c = (l.start + l.end) / 2;
  const w = Math.min(len, z === 'base' ? 60 : Math.max(20, (l.end - l.start) / +z));
  if (w >= len) return goto(l.chr);                 // zooming out stops at the whole chromosome
  const r = clampRange(l.chr, c - w / 2, c + w / 2);
  goto(`${l.chr}:${r.s + 1}-${r.e}`);
});
function toggle(btn, fn) { btn.onclick = () => { btn.classList.toggle('on'); fn(btn.classList.contains('on')); }; }
toggle($('#tgCenter'), on => { if (!state.browser) return; state.browser.setCenterLineVisibility?.(on) ?? state.browser.centerLineButton?.button?.click(); });
toggle($('#tgCursor'), on => { if (!state.browser) return; state.browser.setCursorGuideVisibility?.(on) ?? state.browser.cursorGuide?.setVisibility(on); });
toggle($('#tgNames'), on => setAllFeatureOpt('gvHideNames', !on));
toggle($('#tgIdeo'), on => setIdeogram(on));
toggle($('#tgIdeoBox'), on => { globalThis.__gvNoIdeoBox = !on; const tv = state.browser?.trackViews.find(t => t.track.type === 'ideogram'); tv?.repaintViews?.(); });
toggle($('#tgArrows'), on => setAllFeatureOpt('gvHideArrows', !on));
toggle($('#tgLabels'), on => { if (!state.browser) return; state.browser.setTrackLabelVisibility?.(on) ?? state.browser.trackLabelControl?.button?.click(); });

// ---------------- hide rules (read display, optionally analyses) ----------------
// rule = { on, where: 'chr1:1-100 chr1:500-900' | '', what: 'both'|'D'|'I'|'MM', size: '2, 5-10, >50' | '', note }
const normChr = c => String(c).replace(/^chr/i, '');
function loadRules() { try { return JSON.parse(localStorage.getItem('genoview.hideRules') || '[]'); } catch { return []; } }
function saveRules() { try { localStorage.setItem('genoview.hideRules', JSON.stringify({ rules: state.hideRules, analyses: state.hideRulesInAnalyses })); } catch {} }
// hide-indel rules were removed from the UI: start clean and drop any rules saved by older versions
state.hideRules = []; state.hideRulesInAnalyses = false; try { localStorage.removeItem('genoview.hideRules'); } catch {}
let compiledRules = [];
function compileRules() {
  compiledRules = state.hideRules.filter(r => r.on !== false).map(r => ({
    what: r.what || 'both',
    regs: parseRegions(r.where || '', '').map(g => ({ c: normChr(g.chr), s: g.s, e: g.e })),
    size: parseSizeFilter(r.size || ''),
  }));
}
// does any rule hide this indel? kind 'D' | 'I'
function ruleHidesIndel(chr, start, len, kind) {
  if (!compiledRules.length) return false;
  const c = normChr(chr), a = start, b = kind === 'I' ? start + 1 : start + len;
  for (const r of compiledRules) {
    if (r.what === 'MM' || (r.what !== 'both' && r.what !== kind)) continue;
    if (r.size && !r.size(len)) continue;
    if (r.regs.length && !r.regs.some(g => g.c === c && a < g.e && b > g.s)) continue;
    return true;
  }
  return false;
}
function ruleHidesMismatch(chr, pos) {
  const c = normChr(chr);
  return compiledRules.some(r => r.what === 'MM' && (!r.regs.length || r.regs.some(g => g.c === c && pos >= g.s && pos < g.e)));
}
state.ruleHidesIndel = (chr, start, len, kind) => state.hideRulesInAnalyses && ruleHidesIndel(chr, start, len, kind);
function applyHideIndel() {
  compileRules();
  const anyIndel = compiledRules.some(r => r.what !== 'MM'), anyMM = compiledRules.some(r => r.what === 'MM');
  globalThis.__gvHide = anyIndel ? (chr, start, len, kind) => ruleHidesIndel(chr, start, len, kind) : null;
  globalThis.__gvHideMM = anyMM ? (chr, pos) => ruleHidesMismatch(chr, pos) : null;
  for (const r of state.tracks) { const a = r.igv?.alignmentTrack; if (!a) continue; a.hideSmallIndels = false; r.igv.trackView?.repaintViews(); }
  const n = compiledRules.length;
  $('#btnHideRules')?.classList.toggle('on', n > 0);
  if ($('#hideRulesInfo')) $('#hideRulesInfo').textContent = n ? `${n} rule${n > 1 ? 's' : ''} active${state.hideRulesInAnalyses ? ' · also in analyses' : ''}` : '';
  if (state.hideRulesInAnalyses) { state.onRulesChanged?.(); if (state.result) renderCalls(); }
}

function openHideRules() {
  document.getElementById('rulesModal')?.remove();
  const m = document.createElement('div'); m.id = 'rulesModal'; m.className = 'figmodal';
  m.innerHTML = `<div class="figbox" style="width:min(1100px,96vw);height:auto;max-height:90vh">
    <div class="fighead"><b>Hide indels / mismatches in reads</b><span class="muted" style="margin-left:8px">each rule: where × what × which sizes — add as many as you need</span><span style="flex:1"></span><button id="rmClose" class="primary">Done</button></div>
    <div style="padding:10px 14px;overflow:auto">
      <table class="rules"><thead><tr><th style="width:34px">On</th><th>Where (regions, space-separated · empty = everywhere)</th><th style="width:170px">What</th><th style="width:150px">Sizes (bp)</th><th style="width:180px">Note</th><th style="width:34px"></th></tr></thead><tbody id="rmBody"></tbody></table>
      <div class="rulebar">
        <button id="rmAdd" class="primary">+ Add rule</button>
        <button id="rmView" title="New rule covering the region currently shown">+ Rule for current view</button>
        <button id="rmCrispr" title="One exact rule per indel you hid with ⊘ in the CRISPR tab">+ From CRISPR hidden indels</button>
        <button id="rmNoise" title="Hide all 1–2 bp indels everywhere (long-read noise)">+ Long-read noise (1–2 bp everywhere)</button>
        <span style="flex:1"></span>
        <button id="rmImport">Import TSV/BED…</button><input type="file" id="rmFile" accept=".tsv,.txt,.bed" hidden><button id="rmExport">Export</button><button id="rmClear">Clear all</button>
      </div>
      <label class="status"><input type="checkbox" id="rmAnalyses" ${state.hideRulesInAnalyses ? 'checked' : ''}> <b>Also apply to analyses</b> — rule-matched indels are excluded from CRISPR edit results and hidden from the Variants &amp; SVs list (re-run or change any filter to refresh)</label>
      <div class="hint" style="margin-top:8px">Sizes: <code>2</code> · <code>1-3</code> · <code>&lt;=5</code> · <code>&gt;50</code> · <code>2, 7, 10-12</code> · empty = any size. Regions: <code>chr1:100,030-100,040</code> (several separated by spaces) · empty = everywhere. A deletion is matched if it overlaps the region; an insertion if it sits inside it. Rules are kept in this browser.</div>
    </div></div>`;
  document.body.appendChild(m);
  const body = m.querySelector('#rmBody');
  let t = null; const live = () => { clearTimeout(t); t = setTimeout(applyHideIndel, 250); };
  const draw = () => {
    body.innerHTML = state.hideRules.length ? state.hideRules.map((r, i) => `<tr class="${r.on === false ? 'off' : ''}" data-i="${i}">
      <td><input type="checkbox" data-f="on" ${r.on === false ? '' : 'checked'}></td>
      <td><input type="text" data-f="where" value="${esc(r.where || '')}" placeholder="everywhere" style="font-family:var(--mono)"></td>
      <td><select data-f="what">${[['both', 'Deletions + insertions'], ['D', 'Deletions'], ['I', 'Insertions'], ['MM', 'Mismatches (SNPs)']].map(([v, l]) => `<option value="${v}" ${(r.what || 'both') === v ? 'selected' : ''}>${l}</option>`).join('')}</select></td>
      <td><input type="text" data-f="size" value="${esc(r.size || '')}" placeholder="any size" ${r.what === 'MM' ? 'disabled' : ''}></td>
      <td><input type="text" data-f="note" value="${esc(r.note || '')}" placeholder="e.g. parental 2-bp del"></td>
      <td><button class="ghost x" data-del="${i}" title="Delete rule">×</button></td></tr>`).join('') : '<tr><td colspan="6" class="empty">No rules yet — add one below.</td></tr>';
    body.querySelectorAll('[data-f]').forEach(el => el.addEventListener(el.type === 'text' ? 'input' : 'change', () => {
      const i = +el.closest('tr').dataset.i, f = el.dataset.f;
      state.hideRules[i][f] = el.type === 'checkbox' ? el.checked : el.value;
      if (f === 'on') el.closest('tr').classList.toggle('off', !el.checked);
      if (f === 'what') { el.closest('tr').querySelector('[data-f="size"]').disabled = el.value === 'MM'; }
      live();
    }));
    body.querySelectorAll('[data-del]').forEach(b => b.onclick = () => { state.hideRules.splice(+b.dataset.del, 1); draw(); live(); });
  };
  const add = r => { state.hideRules.push({ on: true, where: '', what: 'both', size: '', note: '', ...r }); draw(); live(); body.querySelector('tr:last-child [data-f="where"]')?.focus(); };
  m.querySelector('#rmAdd').onclick = () => add({});
  m.querySelector('#rmView').onclick = () => { const l = currentLocus(); if (l) add({ where: `${l.chr}:${fmt(l.start + 1)}-${fmt(l.end)}` }); };
  m.querySelector('#rmNoise').onclick = () => add({ where: '', what: 'both', size: '1-2', note: 'long-read noise' });
  m.querySelector('#rmCrispr').onclick = () => {
    const mask = state.crispr?.mask; if (!mask || !mask.size) return toast('No indels hidden with ⊘ in the CRISPR tab yet.', true);
    for (const k of mask.keys()) { const [kind, pos, len] = k.split(':'); const p = +pos, L = +len; const chr = state.crispr.sites[0]?.chr || currentLocus().chr;
      add({ where: kind === 'D' ? `${chr}:${p + 1}-${p + L}` : `${chr}:${p}-${p + 1}`, what: kind, size: String(L), note: mask.get(k) }); }
  };
  m.querySelector('#rmClear').onclick = () => { if (!state.hideRules.length || confirm('Delete all rules?')) { state.hideRules = []; draw(); live(); } };
  m.querySelector('#rmAnalyses').onchange = e => { state.hideRulesInAnalyses = e.target.checked; applyHideIndel(); state.onRulesChanged?.(); if (state.result) renderCalls(); };
  m.querySelector('#rmExport').onclick = () => download('genoview_hide_rules.tsv', '#chr\tstart\tend\twhat\tsizes\tnote\n' + state.hideRules.flatMap(r => {
    const regs = parseRegions(r.where || '', ''); const row = g => [g ? g.chr : '*', g ? g.s : '', g ? g.e : '', r.what || 'both', r.size || '', r.note || ''].join('\t');
    return regs.length ? regs.map(row) : [row(null)]; }).join('\n') + '\n', 'text/tab-separated-values');
  m.querySelector('#rmImport').onclick = () => m.querySelector('#rmFile').click();
  m.querySelector('#rmFile').onchange = async e => {
    let n = 0;
    for (const f of e.target.files) for (const line of (await f.text()).split('\n')) {
      if (!line.trim() || line[0] === '#' || /^(track|browser)/.test(line)) continue;
      const t = line.split('\t');
      const what = ({ del: 'D', deletion: 'D', d: 'D', ins: 'I', insertion: 'I', i: 'I', both: 'both', indel: 'both', mm: 'MM', mismatch: 'MM', snp: 'MM' })[(t[3] || 'both').toLowerCase()] || 'both';
      const where = t[0] === '*' || t[0] === '' ? '' : `${t[0]}:${+t[1] + 1}-${t[2]}`;
      state.hideRules.push({ on: true, where, what, size: t[4] || '', note: t[5] || f.name }); n++;
    }
    e.target.value = ''; draw(); live(); toast(`Imported ${n} rule${n === 1 ? '' : 's'}.`);
  };
  const close = () => { applyHideIndel(); m.remove(); };
  m.querySelector('#rmClose').onclick = close;
  m.addEventListener('mousedown', e => { if (e.target === m) close(); });
  draw();
  if (!state.hideRules.length) add({});
}
if ($('#btnHideRules')) $('#btnHideRules').onclick = openHideRules;
setTimeout(applyHideIndel);

// regions of interest (guides, cut sites, selected calls…) in export coordinates, matching what is on screen
async function collectROIOverlays() {
  const br = state.browser; const out = [];
  try {
    const cc = br.columnContainer.getBoundingClientRect();
    const tv = br.trackViews.find(t => t.viewports?.length) || br.trackViews[0];
    const frames = br.referenceFrameList.map((rf, i) => {
      const r = tv.viewports[i].viewportElement.getBoundingClientRect();
      return { rf, x0: r.left - cc.left, w: r.width, top: 0 };
    });
    for (const set of br.roiManager?.roiSets || []) {
      const rects = [];
      for (const f of frames) {
        const end = f.rf.start + f.w * f.rf.bpPerPixel;
        let feats = [];
        try { feats = await set.getFeatures(f.rf.chr, Math.floor(f.rf.start), Math.ceil(end)); } catch {}
        for (const ft of feats || []) {
          const a = Math.max(ft.start, f.rf.start), b = Math.min(ft.end, end);
          if (b <= a) continue;
          const x = f.x0 + (a - f.rf.start) / f.rf.bpPerPixel, w = (b - a) / f.rf.bpPerPixel;
          rects.push({ x, w, label: ft.name || '', clipL: f.x0, clipR: f.x0 + f.w });
        }
      }
      if (rects.length) out.push({ name: set.name || 'Regions', color: set.color || 'rgba(68,134,247,0.25)', rects });
    }
  } catch (e) { console.warn('ROI export', e); }
  return out;
}

$('#btnSVG').onclick = async () => {
  if (!state.browser) return toast('Nothing to export yet.', true);
  const btn = $('#btnSVG'); btn.disabled = true; btn.textContent = 'Capturing…';
  try {
    const svg = await state.browser.toSVG();
    const overlays = await collectROIOverlays();
    const l = currentLocus();
    openAnnotator({ svg, overlays, key: 'browser', title: `Annotate & export · ${l.chr}:${fmt(l.start + 1)}-${fmt(l.end)}`, filename: `genoview_${l.chr}_${l.start + 1}-${l.end}` });
  } catch (e) { console.error(e); toast('Could not capture the view: ' + e.message, true); }
  finally { btn.disabled = false; btn.textContent = 'Export view…'; }
};

// ---------------- files: drop / pick ----------------
const drop = $('#drop'), fi = $('#fileInput');
drop.onclick = () => { fi.removeAttribute('data-genome'); fi.click(); };
$('#btnGenome').onclick = () => { fi.click(); };
fi.onchange = () => { handleFiles(fi.files); fi.value = ''; };
['dragenter', 'dragover'].forEach(ev => document.addEventListener(ev, e => { e.preventDefault(); drop.classList.add('over'); }));
['dragleave', 'drop'].forEach(ev => document.addEventListener(ev, e => { e.preventDefault(); if (ev === 'drop' || e.target === document.documentElement) drop.classList.remove('over'); }));
document.addEventListener('drop', e => { drop.classList.remove('over'); if (e.dataTransfer?.files?.length) handleFiles(e.dataTransfer.files); });

// hosted UCSC assemblies
$('#btnHosted').onclick = () => { $('#hostedModal').style.display = 'flex'; };
$('#hostedCancel').onclick = () => { $('#hostedModal').style.display = 'none'; };
$('#hostedGo').onclick = async () => {
  const db = $('#hostedFree').value.trim() || $('#hostedSel').value;
  $('#hostedModal').style.display = 'none';
  const base = `https://hgdownload.soe.ucsc.edu/goldenPath/${db}/bigZips/`;
  try {
    toast('Loading ' + db + ' from UCSC…', false, 8000);
    await setGenome({ id: db, name: db, twoBitURL: base + db + '.2bit', chromSizesURL: base + db + '.chrom.sizes' });
    state.genomeSrc = { type: 'ucsc', db };
  } catch (e) { toast('Could not load ' + db + ' (' + (e.message || e) + '). Check the name and your connection.', true, 9000); }
};

// ---------------- dock: tabs + splitter ----------------
$$('.tab').forEach(t => t.onclick = () => {
  $$('.tab').forEach(x => x.classList.toggle('on', x === t));
  $$('.pane').forEach(p => p.classList.toggle('on', p.id === 'pane-' + t.dataset.tab));
  if (t.dataset.tab === 'depth') drawDepth();
  crCtx.onTab?.(t.dataset.tab);
});
function setFullPage(on) {
  // keep the browser at its real size off-screen, so positions / zoom stay valid while the panels use the window
  if (on) { const r = $('#igvWrap').getBoundingClientRect(); $('#work').style.setProperty('--igvW', r.width + 'px'); $('#work').style.setProperty('--igvH', Math.max(200, r.height) + 'px'); }
  $('#work').classList.toggle('full', on);
  $('#btnFull').textContent = on ? '⤡ Back to browser' : '⤢ Full page';
  $('#btnFull').classList.toggle('on', on);
  if (!on) setTimeout(() => { try { state.browser?.visibilityChange?.(); } catch {} window.dispatchEvent(new Event('resize')); }, 60);
  else setTimeout(() => window.dispatchEvent(new Event('resize')), 60);
}
$('#btnFull').onclick = () => setFullPage(!$('#work').classList.contains('full'));
function showTab(name) { $(`.tab[data-tab="${name}"]`).click(); }
// in-page help links (e.g. "How is this calculated?") open the Help tab at that section
document.addEventListener('click', e => { const a = e.target.closest?.('a[href^="#help"]'); if (!a) return; e.preventDefault(); showTab('help'); setTimeout(() => document.getElementById(a.getAttribute('href').slice(1))?.scrollIntoView({ block: 'start' }), 50); });
$('#btnHelp').onclick = () => showTab('help');
(() => {
  const sp = $('#splitter'), work = $('#work'); let drag = false;
  sp.onmousedown = e => { drag = true; e.preventDefault(); };
  window.addEventListener('mousemove', e => { if (!drag) return; const r = work.getBoundingClientRect(); const h = Math.min(r.height - 80, Math.max(90, r.bottom - e.clientY)); work.style.setProperty('--dockH', h + 'px'); });
  window.addEventListener('mouseup', () => { if (drag) { drag = false; if ($('#pane-depth').classList.contains('on')) drawDepth(); } });
})();

// ---------------- analysis ----------------
$('#pRegion').onchange = e => { $('#pCustom').style.display = e.target.value === 'custom' ? '' : 'none'; $('#btnRun').textContent = e.target.value === 'view' ? 'Analyse view' : 'Analyse region'; };
$('#pLib').onchange = e => { if (e.target.value === 'rna') toast('RNA-seq mode: N operations become splice junctions; intron-spanning pairs are not called as deletions.'); };

function parseRegion(str) {
  const m = /^\s*([^:\s]+)(?::([\d,]+)(?:-([\d,]+))?)?\s*$/.exec(str || '');
  if (!m) return null;
  const chr = m[1]; const len = state.browser.genome.getChromosome(chr)?.bpLength;
  if (!len && !m[2]) return null;
  const s = m[2] ? Math.max(0, +m[2].replace(/,/g, '') - 1) : 0;
  const e = m[3] ? +m[3].replace(/,/g, '') : (len || s + 1);
  return { chr: state.browser.genome.getChromosomeName?.(chr) || chr, start: s, end: Math.min(e, len || e) };
}

$('#btnCancel').onclick = () => state.abort?.abort();
$('#btnRun').onclick = runAnalysis;

async function runAnalysis() {
  if (!state.browser) return toast('Load a genome and BAM files first.', true);
  await Promise.all(state.tracks.filter(t => t.sample?.ready).map(t => t.sample.ready));
  const recs = state.tracks.filter(t => t.sample && t.sample.use && t.sample.bam);
  if (!recs.length) return toast('Tick “analyse” on at least one indexed BAM in the track list.', true);
  let reg;
  const mode = $('#pRegion').value;
  if (mode === 'view') reg = currentLocus();
  else if (mode === 'chr') { const l = currentLocus(); reg = { chr: l.chr, start: 0, end: state.browser.genome.getChromosome(l.chr).bpLength }; }
  else reg = parseRegion($('#pCustom').value);
  if (!reg || reg.end <= reg.start) return toast('Could not read the region.', true);
  if (reg.chr === 'all') return toast('Zoom into one chromosome first.', true);
  const span = reg.end - reg.start;
  if (span > 30e6 && !confirm(`Scanning ${fmtLen(span)} across ${recs.length} BAM(s) can take a while. Continue?`)) return;

  const samples = recs.map(r => ({ id: r.id, name: r.name, group: r.sample.group, isControl: r.sample.isControl, bam: r.sample.bam, chroms: r.sample.chroms, refNames: r.sample.refNames }));
  const options = {
    library: $('#pLib').value, minMapq: +$('#pMapq').value, minSV: +$('#pMinSV').value, minSupport: +$('#pSupport').value,
    minAF: +$('#pAF').value, insertSD: +$('#pK').value, skipDuplicates: $('#pDup').checked,
  };
  state.abort = new AbortController();
  $('#btnRun').disabled = true; $('#btnCancel').style.display = ''; $('#prog').style.display = 'block';
  const t0 = performance.now();
  try {
    const res = await analyzeRegion({ samples, chr: reg.chr, start: reg.start, end: reg.end, options, signal: state.abort.signal,
      onProgress: p => { $('#prog > div').style.width = (p.fraction * 100).toFixed(1) + '%'; $('#runStatus').textContent = `${p.sample} · ${fmt(p.at)}`; } });
    res.colors = Object.fromEntries(res.groups.map((g, i) => [g, SERIES[i % SERIES.length]]));
    res.elapsed = (performance.now() - t0) / 1000;
    state.result = res; state.selected = null; detailCache = null;
    $('#detail').style.display = 'none'; $('.callsplit').classList.remove('has');
    const missing = res.perSample.filter(p => !p.bamChr).map(p => p.name);
    const ins = res.perSample.filter(p => p.insert).map(p => `${p.name}: ${Math.round(p.insert.median)}±${Math.round(p.insert.mad)}`);
    $('#runStatus').innerHTML = `${res.calls.length} calls + ${res.cnv.length} depth segments in ${res.chr}:${fmt(res.start + 1)}-${fmt(res.end)} · ${res.elapsed.toFixed(1)} s · bin ${fmt(res.binSize)} bp` +
      (ins.length ? ` · <span title="${esc(ins.join('\n'))}">insert size ✓</span>` : ' · <span class="muted" title="No properly paired reads found: discordant-pair evidence off">single-end / no pairs</span>') +
      (missing.length ? ` · <span style="color:var(--del)">no ${esc(res.chr)} in ${esc(missing.join(', '))}</span>` : '');
    renderCalls(); renderJunctions();
    if ($('#pane-depth').classList.contains('on')) drawDepth();
    if (!res.calls.length && !res.cnv.length) toast('No variants passed the filters in this region.');
  } catch (e) {
    if (e.message !== 'cancelled') { console.error(e); toast('Analysis failed: ' + e.message, true, 9000); } else $('#runStatus').textContent = 'Cancelled.';
  } finally {
    $('#btnRun').disabled = false; $('#btnCancel').style.display = 'none'; $('#prog').style.display = 'none'; $('#prog > div').style.width = 0;
  }
}

// rows = breakpoint/CIGAR calls + depth-only segments
function allRows() {
  const r = state.result; if (!r) return [];
  const rows = r.calls.map(c => ({ ...c, source: c.small ? 'indel' : 'sv', klass: c.small ? 'small ' + c.type : c.type }));
  r.cnv.forEach((s, i) => rows.push({ id: 'cnv' + (i + 1), type: s.type === 'DEL' ? 'LOSS' : 'GAIN', klass: 'CNV', source: 'depth', chr: s.chr, start: s.start, end: s.end, svlen: s.end - s.start, score: Math.round(Math.abs(s.log2) * s.bins), group: s.group, log2: s.log2, cn: s.cn, supportedBy: s.supportedBy, perGroup: {} }));
  return rows;
}
const rowKey = r => `${r.klass}:${r.chr}:${r.start}:${r.end}:${r.group || ''}`;
function parseRegions(txt, defChr) {
  return String(txt || '').split(/[\s;]+|,(?!\d{3}(?!\d))/).filter(Boolean).map(t => { const m = /^(?:([^:\s]+):)?([\d,]+)-([\d,]+)$/.exec(t); return m ? { chr: m[1] || defChr, s: +m[2].replace(/,/g, '') - 1, e: +m[3].replace(/,/g, '') } : null; }).filter(Boolean);
}
function maskedBy(r) {
  if (!state.maskIntervals.length) return false;
  const len = Math.max(1, r.end - r.start);
  return state.maskIntervals.some(m => {
    if (m.chr !== r.chr && m.chr !== 'chr' + r.chr && 'chr' + m.chr !== r.chr) return false;
    if (r.small || len <= 50) return Math.abs(m.s - r.start) <= 3 && (!m.type || m.type === r.type || r.type === 'LOSS' || r.type === 'GAIN');
    const ov = Math.min(m.e, r.end) - Math.max(m.s, r.start);
    return ov > 0 && ov >= 0.5 * len && ov >= 0.5 * (m.e - m.s);
  });
}
// why a row is hidden ('' = visible)
function hiddenReason(r) {
  if (state.hidden.has(rowKey(r))) return 'hidden by you';
  const sz = parseSizeFilter($('#fSize').value);
  const len = r.type === 'TRA' ? null : Math.abs(r.svlen ?? (r.end - r.start));
  if (sz && len != null && sz(len)) return 'size filter';
  const rois = parseRegions($('#fROI').value, r.chr);
  if (rois.length && !rois.some(g => g.chr === r.chr && r.start < g.e && Math.max(r.end, r.start + 1) > g.s)) return 'outside regions of interest';
  if (maskedBy(r)) return 'known variant';
  if (state.hideRulesInAnalyses && (r.type === 'DEL' || r.type === 'INS' || r.type === 'LOSS') && state.ruleHidesIndel(r.chr, r.start, Math.max(1, Math.abs(r.svlen ?? (r.end - r.start))), r.type === 'INS' ? 'I' : 'D')) return 'hide rule';
  return '';
}
function visibleRows() { return allRows().filter(r => !hiddenReason(r)); }
function filteredResult() {
  const res = state.result; const keep = new Set(visibleRows().filter(r => r.source !== 'depth').map(r => r.id));
  return { ...res, calls: res.calls.filter(c => keep.has(c.id)), cnv: res.cnv.filter((c, i) => !hiddenReason(allRows().find(r => r.id === 'cnv' + (i + 1)))) };
}
['#fSize', '#fROI'].forEach(id => $(id).addEventListener('change', () => { if (state.result) { renderCalls(); if ($('#pane-depth').classList.contains('on')) drawDepth(); } }));
$('#fShowHidden').onchange = () => state.result && renderTableBody();
$('#fMaskBtn').onclick = () => $('#fMaskFile').click();
$('#fMaskFile').onchange = async e => {
  const ivs = [];
  for (const f of e.target.files) {
    const txt = /\.gz$/i.test(f.name) ? await new Response(f.stream().pipeThrough(new DecompressionStream('gzip'))).text() : await f.text();
    for (const line of txt.split('\n')) {
      if (!line || line[0] === '#' || line.startsWith('track') || line.startsWith('browser')) continue;
      const t = line.split('\t');
      if (/\.vcf(\.gz)?$/i.test(f.name)) {
        const pos = +t[1] - 1, ref = t[3] || '', alts = (t[4] || '').split(',');
        const end = /END=(\d+)/.exec(t[7] || ''); const svt = /SVTYPE=(\w+)/.exec(t[7] || '');
        for (const alt of alts) {
          if (end) ivs.push({ chr: t[0], s: pos, e: +end[1], type: svt?.[1] });
          else if (ref.length > alt.length) ivs.push({ chr: t[0], s: pos + 1, e: pos + ref.length, type: 'DEL' });
          else if (alt.length > ref.length) ivs.push({ chr: t[0], s: pos + 1, e: pos + 2, type: 'INS' });
          else ivs.push({ chr: t[0], s: pos, e: pos + ref.length });
        }
      } else if (t.length >= 3 && !isNaN(+t[1])) ivs.push({ chr: t[0], s: +t[1], e: +t[2] });
    }
  }
  state.maskIntervals = ivs; e.target.value = '';
  $('#fMaskInfo').innerHTML = ivs.length ? `${ivs.length.toLocaleString()} known variants loaded <button class="ghost x" id="fMaskClear" title="Remove mask">×</button>` : '';
  if ($('#fMaskClear')) $('#fMaskClear').onclick = () => { state.maskIntervals = []; $('#fMaskInfo').textContent = ''; if (state.result) renderCalls(); };
  if (state.result) renderCalls();
};

const KLASSES = ['small DEL', 'small INS', 'DEL', 'DUP', 'INV', 'INS', 'TRA', 'CNV'];
const klassColor = k => k === 'CNV' ? TYPE_COLORS.CNV : TYPE_COLORS[k.replace('small ', '')] || '#8a8984';

function renderCalls() {
  const res = state.result; const rows = visibleRows();
  const nHidden = allRows().length - rows.length;
  $('#nCalls').textContent = rows.length;
  $('#fHiddenInfo').textContent = nHidden ? `${nHidden} hidden` : '';
  const counts = {}; rows.forEach(r => counts[r.klass] = (counts[r.klass] || 0) + 1);
  $('#chips').innerHTML = KLASSES.filter(k => counts[k]).map(k => `<button class="chip ${state.typeFilter.has(k) ? 'off' : ''}" data-k="${k}"><span class="swatch" style="background:${klassColor(k)}"></span>${k} <b>${counts[k]}</b></button>`).join('') +
    `<span style="flex:1"></span><input type="text" id="callSearch" placeholder="filter…" style="width:120px" value="${esc(state.search || '')}">
     <button id="exSpec" class="primary" title="Size distribution of the visible calls, small → large; editable publication figure">Size spectrum…</button><button id="exVCF">VCF</button><button id="exBED">BED</button><button id="exBEDPE">BEDPE</button><button id="exTSV">TSV</button><button id="exShow" class="primary" title="Add calls as a BED track and highlight them as regions of interest">Show in browser</button>`;
  $$('#chips .chip').forEach(c => c.onclick = () => { const k = c.dataset.k; state.typeFilter.has(k) ? state.typeFilter.delete(k) : state.typeFilter.add(k); renderCalls(); });
  $('#callSearch').oninput = e => { state.search = e.target.value; renderTableBody(); };
  const base = `genoview_${res.chr}_${res.start + 1}-${res.end}`;
  $('#exVCF').onclick = () => download(base + '.vcf', callsToVCF(filteredResult()), 'text/plain');
  $('#exBED').onclick = () => { const f = filteredResult(); download(base + '.bed', callsToBED(f) + cnvBED(f), 'text/plain'); };
  $('#exBEDPE').onclick = () => download(base + '.bedpe', callsToBEDPE(filteredResult()), 'text/plain');
  $('#exTSV').onclick = () => download(base + '.tsv', callsToTSV(filteredResult()), 'text/tab-separated-values');
  $('#exShow').onclick = showCallsInBrowser;
  $('#exSpec').onclick = openSvSpectrum;

  const groups = res.groups;
  const head = [['type', 'Type'], ['start', 'Location'], ['svlen', 'Size', 'num'], ['score', 'Score', 'num'], ['ev', 'Evidence  SR / PE / CG / clip'],
    ...groups.map(g => [`g:${g}`, `${g}${res.control === g ? ' (ctl)' : ''}`, 'num']), ['hide', '']];
  $('#callsTable').innerHTML = `<table class="calls"><thead><tr>${head.map(([k, l, c]) => `<th data-k="${k}" class="${c || ''}">${esc(l)}${state.sort.key === k ? (state.sort.dir > 0 ? ' ▲' : ' ▼') : ''}</th>`).join('')}</tr></thead><tbody></tbody></table>
    <div class="hint" style="margin-top:6px">Per-sample cells: supporting reads · allele fraction · depth ratio inside/flanks (DEL ≈ 0.5 het / 0 hom; DUP ≈ 1.5 het). Click a row to jump there and draw its depth graph.</div>`;
  $$('#callsTable th').forEach(th => th.onclick = () => { const k = th.dataset.k; state.sort = { key: k, dir: state.sort.key === k ? -state.sort.dir : (k === 'score' || k.startsWith('g:') ? -1 : 1) }; renderCalls(); });
  renderTableBody();
}

const SV_SERIES = [['DEL', 'Deletions', '#e34948'], ['INS', 'Insertions', '#1baf7a'], ['DUP', 'Duplications', '#2a78d6'], ['INV', 'Inversions', '#6250d6'], ['LOSS', 'CNV loss (depth)', '#eda100'], ['GAIN', 'CNV gain (depth)', '#8a6d00']];
function svSpec(rows, name) {
  const series = SV_SERIES.map(([key, label, color]) => ({ key, label, color, map: new Map() }));
  let n = 0;
  for (const r of rows) {
    if (r.type === 'TRA') continue;
    const se = series.find(x => x.key === r.type); if (!se) continue;
    const len = Math.max(1, Math.abs(r.svlen ?? (r.end - r.start)));
    se.map.set(len, (se.map.get(len) || 0) + 1); n++;
  }
  return { name, unitLabel: 'calls', total: n, edited: n, events: n, series: series.filter(x => x.map.size || ['DEL', 'INS'].includes(x.key)) };
}
function openSvSpectrum() {
  const res = state.result; if (!res) return;
  const rows = visibleRows().filter(r => !state.typeFilter.has(r.klass));
  const inGroup = (r, g) => r.source === 'depth' ? r.group === g : (r.perGroup[g]?.support || 0) > 0;
  const groups = res.groups;
  const opt = [{ value: '__all', label: `All calls (${groups.length > 1 ? 'all samples/groups, each call once' : groups[0]})` }, ...groups.map(g => ({ value: g, label: g }))];
  if (groups.length > 1) opt.push({ value: '__mean', label: `Mean of ${groups.length} samples/groups ± SD` });
  const region = `${res.chr}:${fmt(res.start + 1)}-${fmt(res.end)}`;
  const getSpec = sel => {
    let sp;
    if (sel.src === '__all') { sp = svSpec(rows, 'all calls'); sp.subtitle = `${sp.total} calls · ${groups.length} sample/group${groups.length > 1 ? 's' : ''} · ${region}`; }
    else if (sel.src === '__mean') { const per = groups.map(g => svSpec(rows.filter(r => inGroup(r, g)), g)); sp = poolSpecs(per, 'mean', 'mean'); sp.subtitle = `Mean of ${groups.length} samples/groups · ${region}`; }
    else { sp = svSpec(rows.filter(r => inGroup(r, sel.src)), sel.src); sp.subtitle = `${sel.src} · ${sp.total} calls · ${region}`; }
    return sp;
  };
  openFigureEditor({ key: 'sv', title: 'Indel / SV size spectrum — figure editor', selects: [{ id: 'src', label: 'Calls from', options: opt, value: '__all' }], getSpec,
    defaults: { title: 'Indel and SV size spectrum', xLabel: 'Event length (bp)', measure: 'count', layout: 'grouped', showValues: true, bins: 'custom', customEdges: '20, 50, 100, 500, 1000, 5000, 10000, 100000' },
    filename: sel => `sv_spectrum_${sel.src === '__all' ? 'all' : sel.src === '__mean' ? 'mean' : sel.src}_${res.chr}_${res.start + 1}-${res.end}` });
}

function cnvBED(res) { return res.cnv.map((s, i) => `${s.chr}\t${s.start}\t${s.end}\tCNV_${s.type === 'DEL' ? 'loss' : 'gain'}_${s.group}_CN${s.cn}\t${Math.min(1000, Math.round(Math.abs(s.log2) * s.bins))}`).join('\n') + (res.cnv.length ? '\n' : ''); }

function sortVal(r, k) {
  if (k === 'type') return r.klass; if (k === 'start') return r.start; if (k === 'svlen') return Math.abs(r.svlen ?? 0); if (k === 'score') return r.score;
  if (k === 'ev') return (r.evidence?.split || 0) + (r.evidence?.disc || 0) + (r.evidence?.cigar || 0);
  if (k === 'hide') return hiddenReason(r) ? 1 : 0;
  if (k.startsWith('g:')) { const g = k.slice(2); return r.source === 'depth' ? (r.group === g ? Math.abs(r.log2) * 100 : 0) : (r.perGroup[g]?.support || 0); }
  return 0;
}

function renderTableBody() {
  const res = state.result; if (!res) return;
  const q = (state.search || '').toLowerCase();
  const showHidden = $('#fShowHidden').checked;
  let rows = allRows().filter(r => !state.typeFilter.has(r.klass) && (showHidden || !hiddenReason(r)));
  if (q) rows = rows.filter(r => `${r.klass} ${r.chr}:${r.start + 1}-${r.end} ${r.chr2 || ''} ${r.group || ''} ${r.id}`.toLowerCase().includes(q));
  const { key, dir } = state.sort;
  rows.sort((a, b) => { const x = sortVal(a, key), y = sortVal(b, key); return (x < y ? -1 : x > y ? 1 : 0) * dir || a.start - b.start; });
  state.visibleRows = rows;
  const maxS = Math.max(1, ...rows.map(r => r.score || 0));
  const tb = $('#callsTable tbody');
  tb.innerHTML = rows.slice(0, 3000).map(r => {
    const col = klassColor(r.klass);
    const loc = r.type === 'TRA' ? `${r.chr}:${fmt(r.start + 1)} → ${r.chr2}:${fmt(r.pos2 + 1)}` : r.small && r.type === 'INS' ? `${r.chr}:${fmt(r.start + 1)}` : `${r.chr}:${fmt(r.start + 1)}-${fmt(r.end)}`;
    const size = r.type === 'TRA' ? '—' : r.svlen == null ? '?' : (r.svlen > 0 && r.type !== 'DEL' && r.type !== 'LOSS' ? '+' : r.type === 'DEL' || r.type === 'LOSS' ? '−' : '') + fmtLen(r.svlen);
    const ev = r.source === 'depth' ? `<span class="ev">depth · log2 ${r.log2.toFixed(2)} · CN≈${r.cn}${r.supportedBy ? ' · ✓ breakpoints' : ''}</span>`
      : `<span class="ev">${r.evidence.split || 0} / ${r.evidence.disc || 0} / ${r.evidence.cigar || 0} / ${r.clipSupport ?? '–'}</span>${r.precise ? '' : ' <span class="muted" title="Breakpoints from discordant pairs only">≈</span>'}`;
    const cells = res.groups.map(g => {
      if (r.source === 'depth') return `<td class="num">${r.group === g ? `<b>${r.cn}</b> <span class="muted">CN</span>` : '<span class="muted">·</span>'}</td>`;
      const p = r.perGroup[g] || {};
      if (!p.support) return `<td class="num"><span class="muted">0</span>${p.ratio != null && !isNaN(p.ratio) ? ` <span class="muted">· ${p.ratio.toFixed(2)}</span>` : ''}</td>`;
      return `<td class="num"><b>${p.support}</b> <span class="muted">· ${(p.af ?? 0).toFixed(2)}${p.ratio != null && !isNaN(p.ratio) ? ' · ' + p.ratio.toFixed(2) : ''}</span></td>`;
    }).join('');
    const label = r.source === 'depth' ? `CNV ${r.type.toLowerCase()}` : r.klass;
    const why = hiddenReason(r);
    const hideBtn = why ? `<button class="ghost x unhide" data-k="${esc(rowKey(r))}" title="Hidden (${why}) — click to show${why === 'hidden by you' ? '' : '; change the filter to un-hide'}">↺</button>` : `<button class="ghost x hide" data-k="${esc(rowKey(r))}" title="Hide this call (polymorphism / not of interest)">⊘</button>`;
    return `<tr data-id="${r.id}" class="${state.selected?.id === r.id ? 'sel' : ''}" style="${why ? 'opacity:.45' : ''}"><td><span class="pill"><i style="background:${col}"></i>${label}</span></td><td class="loc">${loc}</td><td class="num">${size}</td>
      <td class="num"><span class="minibar" style="width:${Math.max(2, 36 * (r.score || 0) / maxS)}px;background:${col}"></span>${r.score}</td><td>${ev}</td>${cells}<td>${hideBtn}</td></tr>`;
  }).join('') || `<tr><td colspan="${6 + res.groups.length}" class="empty">No calls match.</td></tr>`;
  $$('button.hide, button.unhide', tb).forEach(b => b.onclick = ev => { ev.stopPropagation(); const k = b.dataset.k; state.hidden.has(k) ? state.hidden.delete(k) : state.hidden.add(k); renderCalls(); if ($('#pane-depth').classList.contains('on')) drawDepth(); });
  $$('tr[data-id]', tb).forEach(tr => tr.onclick = () => selectCall(rows.find(r => r.id === tr.dataset.id)));
}

async function selectCall(r) {
  state.selected = r;
  $$('#callsTable tr').forEach(tr => tr.classList.toggle('sel', tr.dataset.id === r.id));
  const w = r.type === 'TRA' ? 2000 : Math.max(r.end - r.start, 1);
  const pad = r.small ? 60 : Math.max(500, Math.round(w * 0.6));
  if (r.type === 'TRA') await goto(`${r.chr}:${fmt(Math.max(1, r.start - 1000))}-${fmt(r.start + 1000)} ${r.chr2}:${fmt(Math.max(1, r.pos2 - 1000))}-${fmt(r.pos2 + 1000)}`);
  else await goto(`${r.chr}:${Math.max(1, r.start - pad + 1)}-${r.end + pad}`);
  try { state.browser.clearROIs(); await state.browser.loadROI([{ name: 'Selected call', color: hexA(klassColor(r.klass), 0.18), features: [{ chr: r.chr, start: r.start, end: Math.max(r.end, r.start + 1), name: r.klass }] }]); } catch (e) { console.warn(e); }
  await drawDetail(r);
}

function hexA(hex, a) { const n = parseInt(hex.slice(1), 16); return `rgba(${n >> 16},${(n >> 8) & 255},${n & 255},${a})`; }

async function showCallsInBrowser() {
  const res = state.result; if (!res) return;
  const rows = visibleRows().filter(r => !state.typeFilter.has(r.klass));
  const bed = rows.filter(r => r.type !== 'TRA').map(r => `${r.chr}\t${r.start}\t${Math.max(r.end, r.start + 1)}\t${r.source === 'depth' ? 'CNV_' + r.type.toLowerCase() + '_' + r.group : r.klass.replace(' ', '_')}${r.svlen ? '_' + fmtLen(r.svlen).replace(' ', '') : ''}\t${Math.min(1000, r.score)}\t.\t${r.start}\t${Math.max(r.end, r.start + 1)}\t${rgb(klassColor(r.klass))}`).join('\n');
  const old = state.tracks.find(t => t.generated && t.name === 'GenoView calls'); if (old) removeTrackRecord(old, true);
  await addFileTrack(new File([bed + '\n'], 'genoview_calls.bed'), null, { name: 'GenoView calls', format: 'bed', generated: true, cfg: { itemRgb: true, displayMode: 'EXPANDED' } });
  const tra = res.calls.filter(c => c.type === 'TRA');
  if (tra.length) {
    const bedpe = tra.map(c => `${c.chr}\t${c.start}\t${c.start + 1}\t${c.chr2}\t${c.pos2}\t${c.pos2 + 1}\tTRA_${c.id}\t${c.score}`).join('\n');
    const o2 = state.tracks.find(t => t.generated && t.name === 'GenoView translocations'); if (o2) removeTrackRecord(o2, true);
    await addFileTrack(new File([bedpe + '\n'], 'genoview_tra.bedpe'), null, { name: 'GenoView translocations', format: 'bedpe', generated: true, color: TYPE_COLORS.TRA });
  }
  try { state.browser.clearROIs(); const byK = {}; rows.filter(r => r.type !== 'TRA' && !r.small).forEach(r => (byK[r.klass] ||= []).push({ chr: r.chr, start: r.start, end: r.end, name: r.source === 'depth' ? 'CNV ' + r.type.toLowerCase() : r.klass }));
    await state.browser.loadROI(Object.entries(byK).map(([k, f]) => ({ name: 'GenoView ' + k, color: hexA(klassColor(k), 0.12), features: f }))); } catch (e) { console.warn(e); }
  toast('Calls added as tracks; SV spans highlighted. Clear highlights from any track’s right-click menu → Regions of interest.');
}
function rgb(hex) { const n = parseInt(hex.slice(1), 16); return `${n >> 16},${(n >> 8) & 255},${n & 255}`; }

// ---------------- depth graphs ----------------
function drawDepth() {
  const res = state.result; const el = $('#depthChart');
  if (!res) return;
  const mode = $('#gMode').value;
  if (mode === 'log2' && !res.control) { $('#gMode').value = 'norm'; toast('Mark one BAM as control (left panel) to plot log2 ratios. Showing copy number instead.'); }
  $('#depthTitle').textContent = `Read depth · ${res.chr}:${fmt(res.start + 1)}-${fmt(res.end)} · ${fmt(res.binSize)} bp bins`;
  const shade = $('#gCalls').checked ? visibleRows().filter(r => !r.small && r.type !== 'TRA' && !state.typeFilter.has(r.klass)) : [];
  depthChart(el, res, { mode: $('#gMode').value, pooled: $('#gPooled').checked && res.groups.length > 1, shade, selected: state.selected,
    onClick: pos => goto(`${res.chr}:${Math.max(1, Math.round(pos - (res.end - res.start) / 40))}-${Math.round(pos + (res.end - res.start) / 40)}`) });
}
['gMode', 'gPooled', 'gCalls'].forEach(id => $('#' + id).onchange = () => { drawDepth(); if (state.selected) drawDetail(state.selected, true); });
window.addEventListener('resize', () => { if ($('#pane-depth').classList.contains('on')) drawDepth(); if (state.selected) drawDetail(state.selected, true); });

let detailCache = null;
async function drawDetail(r, reuse = false) {
  const res = state.result; if (!res) return;
  $('#detail').style.display = ''; $('.callsplit').classList.add('has');
  let det = reuse && detailCache?.id === r.id ? detailCache.res : null;
  if (!det) {
    const recs = state.tracks.filter(t => t.sample && t.sample.use && t.sample.bam);
    const samples = recs.map(t => ({ id: t.id, name: t.name, group: t.sample.group, isControl: t.sample.isControl, bam: t.sample.bam, chroms: t.sample.chroms, refNames: t.sample.refNames }));
    const span = r.type === 'TRA' ? 0 : r.end - r.start;
    const flank = r.small ? 150 : Math.max(800, Math.round(span * 0.75));
    const chrLen = state.browser.genome.getChromosome(r.chr)?.bpLength || Infinity;
    const s = Math.max(0, r.start - flank), e = Math.min(chrLen, (r.type === 'TRA' ? r.start : r.end) + flank);
    $('#detailTitle').textContent = 'Computing depth around the call…';
    const o = { ...res.options, binSize: Math.max(1, Math.ceil((e - s) / 400)) };
    det = await analyzeRegion({ samples, chr: r.chr, start: s, end: e, options: o });
    det.colors = res.colors;
    if (r.type === 'TRA') {
      const len2 = state.browser.genome.getChromosome(r.chr2)?.bpLength || Infinity;
      det.partner = await analyzeRegion({ samples, chr: r.chr2, start: Math.max(0, r.pos2 - flank), end: Math.min(len2, r.pos2 + flank), options: o });
      det.partner.colors = res.colors;
    }
    detailCache = { id: r.id, res: det };
  }
  const label = r.source === 'depth' ? `CNV ${r.type.toLowerCase()} in ${r.group} (CN≈${r.cn})` : `${r.klass} ${r.type === 'TRA' ? `${r.chr}:${fmt(r.start + 1)} → ${r.chr2}:${fmt(r.pos2 + 1)}` : `${r.chr}:${fmt(r.start + 1)}-${fmt(r.end)} (${fmtLen(r.svlen ?? 0)})`}`;
  $('#detailTitle').textContent = 'Depth around ' + label;
  const m = $('#gMode').value === 'log2' && !det.control ? 'norm' : $('#gMode').value;
  depthChart($('#detailChart'), det, { mode: m, pooled: $('#gPooled').checked && det.groups.length > 1, shade: [r], selected: r, compact: true, medians: res.medianDepth });
  if (det.partner) {
    const d2 = document.createElement('div'); d2.className = 'chart'; d2.style.marginTop = '6px'; $('#detailChart').appendChild(d2);
    depthChart(d2, det.partner, { mode: m, pooled: $('#gPooled').checked && det.groups.length > 1, shade: [{ ...r, chr: r.chr2, start: r.pos2, end: r.pos2 + 1 }], selected: r, compact: true, medians: res.medianDepth });
  }
  if (r.source === 'depth') { $('#evChart').innerHTML = `<div class="hint">Depth-only segment: ${r.bins} bins with mean log2 ratio ${r.log2.toFixed(2)} (≈ ${r.cn} copies) in <b>${esc(r.group)}</b>.${r.supportedBy ? ' Matches a breakpoint call (' + r.supportedBy + ').' : ' No split-read / discordant-pair call matches it: check mappability and GC before trusting it.'}</div>`; $('#evTable').innerHTML = ''; return; }
  evidenceChart($('#evChart'), r, res);
  $('#evTable').innerHTML = `<table class="calls" style="margin-top:6px"><thead><tr><th>Group</th><th class="num">Reads</th><th class="num">AF</th><th class="num">Depth in / flank</th></tr></thead><tbody>${res.groups.map(g => { const p = r.perGroup[g] || {}; return `<tr><td><span class="swatch" style="background:${res.colors[g]}"></span> ${esc(g)}</td><td class="num">${p.support || 0}</td><td class="num">${(p.af ?? 0).toFixed(2)}</td><td class="num">${p.depthInside != null ? `${p.depthInside.toFixed(1)} / ${p.depthFlank.toFixed(1)} = <b>${p.ratio.toFixed(2)}</b>` : '–'}</td></tr>`; }).join('')}</tbody></table>`;
}

$('#btnDepthAnn').onclick = () => { const svg = $('#depthChart svg'); if (!svg) return toast('Run an analysis first.', true); const r = state.result; openAnnotator({ svg: serializeSVG(svg).replace(/^<\?xml[^>]*>\s*/, ''), key: 'depth', title: 'Annotate & export · depth graph', filename: `depth_${r.chr}_${r.start + 1}-${r.end}` }); };
$('#btnDetAnn').onclick = () => { const svg = $('#detailChart svg'); if (!svg) return; openAnnotator({ svg: serializeSVG(svg).replace(/^<\?xml[^>]*>\s*/, ''), key: 'detail', title: 'Annotate & export · call depth', filename: `call_${state.selected?.id || ''}_depth` }); };
$('#btnDetSVG').onclick = () => { const svg = $('#detailChart svg'); if (svg) download(`call_${state.selected?.id || ''}_depth.svg`, serializeSVG(svg), 'image/svg+xml'); };
$('#btnDetPNG').onclick = () => { const svg = $('#detailChart svg'); if (svg) svgToPNG(svg, `call_${state.selected?.id || ''}_depth.png`); };
$('#btnDepthSVG').onclick = () => { const svg = $('#depthChart svg'); if (!svg) return; download(`depth_${state.result.chr}_${state.result.start + 1}-${state.result.end}.svg`, serializeSVG(svg), 'image/svg+xml'); };
$('#btnDepthPNG').onclick = () => { const svg = $('#depthChart svg'); if (!svg) return; svgToPNG(svg, `depth_${state.result.chr}_${state.result.start + 1}-${state.result.end}.png`); };
$('#btnDepthBG').onclick = () => {
  const res = state.result; if (!res) return;
  for (const g of [...res.groups, ...(res.groups.length > 1 ? ['*'] : [])]) download(`depth_${g === '*' ? 'pooled' : g}_${res.chr}_${res.start + 1}-${res.end}.bedGraph`, depthToBedGraph(res, g), 'text/plain');
};
$('#btnDepthTrack').onclick = async () => {
  const res = state.result; if (!res) return toast('Run an analysis first.', true);
  const gs = [...res.groups, ...(res.groups.length > 1 && $('#gPooled').checked ? ['*'] : [])];
  const tracks = gs.map(g => ({ url: new File([depthToBedGraph(res, g)], `${g === '*' ? 'pooled' : g}.bedgraph`), format: 'bedgraph', name: g === '*' ? 'pooled' : g, color: g === '*' ? '#52514e' : res.colors[g] }));
  const name = `Depth ${res.chr}:${fmt(res.start + 1)}-${fmt(res.end)}`;
  const cfg = { type: 'merged', name, tracks, height: 110, autoscale: true, order: state.nextId * 10 };
  const t = await state.browser.loadTrack(cfg);
  state.tracks.push({ id: state.nextId++, name, format: 'merged', igv: t, color: '#52514e', height: 110, mode: 'full', generated: true, cfg });
  renderTrackList(); toast('Depth overlay added to the browser.');
};

// ---------------- junctions ----------------
function renderJunctions() {
  const res = state.result; const j = res?.junctions || [];
  $('#nJunc').textContent = j.length;
  if (!j.length) { $('#juncTable').innerHTML = res?.options.library === 'rna' ? '<div class="empty">No spliced reads in this region.</div>' : ''; return; }
  $('#juncInfo').textContent = `${j.length} junctions (reads with N in CIGAR) in ${res.chr}:${fmt(res.start + 1)}-${fmt(res.end)}`;
  const max = Math.max(...j.map(x => x.total));
  $('#juncTable').innerHTML = `<table class="calls"><thead><tr><th>Junction</th><th class="num">Intron</th><th class="num">Total reads</th>${res.groups.map(g => `<th class="num">${esc(g)}</th>`).join('')}</tr></thead><tbody>${j.map((x, i) => `<tr data-i="${i}"><td class="loc">${x.chr}:${fmt(x.start)}-${fmt(x.end + 1)}</td><td class="num">${fmtLen(x.end - x.start)}</td><td class="num"><span class="minibar" style="width:${Math.max(2, 50 * x.total / max)}px"></span>${x.total}</td>${res.groups.map(g => `<td class="num">${x.support[g] || 0}</td>`).join('')}</tr>`).join('')}</tbody></table>`;
  $$('#juncTable tr[data-i]').forEach(tr => tr.onclick = () => { const x = j[+tr.dataset.i]; goto(`${x.chr}:${x.start - 200}-${x.end + 200}`); });
}
$('#btnJuncBED').onclick = () => { const r = state.result; if (!r?.junctions.length) return toast('No junctions yet.', true); download(`junctions_${r.chr}_${r.start + 1}-${r.end}.bed`, junctionsToBED(r), 'text/plain'); };
$('#btnJuncTrack').onclick = async () => {
  const r = state.result; if (!r?.junctions.length) return toast('No junctions yet.', true);
  const jb = r.junctions.map(j => `${j.chr}\t${j.start}\t${j.end}\tdepth=${j.total};uniquely_mapped=${j.total};multi_mapped=0;annotated_junction=false;${Object.entries(j.support).map(([g, n]) => `${g.replace(/[;= ]/g, '_')}=${n}`).join(';')}\t${j.total}\t+`).join('\n') + '\n';
  const old = state.tracks.find(t => t.generated && t.name === 'Splice junctions'); if (old) removeTrackRecord(old, true);
  await addFileTrack(new File([jb], 'junctions.bed'), null, { name: 'Splice junctions', format: 'bed', generated: true, color: '#6250d6', cfg: { type: 'junction' } });
};

// ---------------- downloads ----------------
function download(name, text, type) {
  const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([text], { type })); a.download = name;
  document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500);
}
function serializeSVG(svg) {
  const c = svg.cloneNode(true); c.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  c.querySelectorAll('.hover').forEach(n => n.remove());
  return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(c);
}
function svgToPNG(svg, name, scale = 3) {
  const vb = svg.viewBox.baseVal; const img = new Image();
  img.onload = () => { const cv = document.createElement('canvas'); cv.width = vb.width * scale; cv.height = vb.height * scale; const ctx = cv.getContext('2d'); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, cv.width, cv.height); ctx.drawImage(img, 0, 0, cv.width, cv.height); cv.toBlob(b => { const a = document.createElement('a'); a.href = URL.createObjectURL(b); a.download = name; a.click(); }); };
  img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(serializeSVG(svg));
}

const crCtx = { openAnnotator, handleFiles, state, $, $$, toast, download, serializeSVG, svgToPNG, goto, currentLocus, esc, fmt };
initCrispr(crCtx);

// ---------------- stay active in background tabs ----------------
// Chrome freezes / throttles hidden tabs; a held Web Lock exempts the page from freezing,
// and a beforeunload prompt protects loaded work from an accidental close or tab discard.
let keepAlive = false;
function ensureKeepAlive() {
  if (keepAlive) return; keepAlive = true;
  try { navigator.locks?.request('genoview-keep-active', () => new Promise(() => {})); } catch (e) { console.warn(e); }
  window.addEventListener('beforeunload', e => { if (state.tracks.length) { e.preventDefault(); e.returnValue = ''; } });
}
document.addEventListener('visibilitychange', () => {
  if (document.visibilityState !== 'visible' || !state.browser) return;
  try { state.browser.visibilityChange?.(); } catch {}
  for (const r of state.tracks) r.igv?.trackView?.repaintViews?.();
});
window.genoviewKeepAlive = ensureKeepAlive;

initSession({ state, $, $$, toast, download, esc, fmt, goto, currentLocus, renderTrackList, renderCalls, drawDepth, applyHideIndel, runAnalysis,
  handleFiles, loadGenomeFromFiles, genomeFromBam, setGenome, addFileTrack, findIndex, isIndex, formatOf, setMode, setColor, setHeight, setCoverageHeight, setBamPart,
  setFeatureOpt, setCompact, sortReads, applyOrder, rebuildOverviewTrack, getAnnotationStore, setAnnotationStore, fmeta });

// expose for debugging / automation
window.genoview = { state, handleFiles, runAnalysis, goto, selectCall, allRows };
