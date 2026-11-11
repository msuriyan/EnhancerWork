import { yieldNow } from './yield.js';
// Fast "where are the reads?" overview for a BAM, usable at any zoom level.
// Small BAMs: scan the reads once (exact depth in bins). Large BAMs: estimate density from the index.

const F_UNMAPPED = 4, F_SECONDARY = 256, F_SUPP = 2048;

export async function computeOverview(bam, { fileSize, maxReadsForScan = 3_000_000, maxBytesForScan = 600e6, onProgress = () => {} } = {}) {
  await bam.getHeader();
  const refs = bam.indexToChr || [];
  const counts = [];
  let total = 0;
  for (let i = 0; i < refs.length; i++) {
    let n = 0; try { n = await bam.index.lineCount(i); } catch { n = -1; }
    counts.push(n); if (n > 0) total += n;
  }
  const scan = fileSize <= maxBytesForScan && total <= maxReadsForScan;
  const perChrom = {};
  let done = 0;
  for (let i = 0; i < refs.length; i++) {
    const { refName: chr, length: len } = refs[i];
    if (counts[i] === 0) continue;
    if (scan) {
      const bin = Math.max(25, Math.ceil(len / 40000));
      const nb = Math.ceil(len / bin);
      const cov = new Float32Array(nb);
      let reads = 0;
      const W = 5_000_000;
      for (let ws = 0; ws < len; ws += W) {
        const recs = await bam.getRecordsForRange(chr, ws, Math.min(len, ws + W));
        for (const r of recs) {
          if (r.start < ws && ws > 0) continue;
          const f = r.flags; if (f & (F_UNMAPPED | F_SECONDARY | F_SUPP)) continue;
          reads++;
          const cig = r.NUMERIC_CIGAR; let p = r.start;
          for (let k = 0; k < cig.length; k++) {
            const op = cig[k] & 0xf, l = cig[k] >>> 4;
            if (op === 0 || op === 7 || op === 8 || op === 2) { addSpan(cov, p, p + l, bin); p += l; }
            else if (op === 3) p += l;
          }
        }
        done += Math.min(W, len - ws);
        onProgress(done);
        await yieldNow();
      }
      if (reads) perChrom[chr] = { bin, cov, reads, len, exact: true };
    } else {
      // index-based density (compressed bytes per 16 kb window) — relative only
      try {
        const ic = await bam.index.indexCov(i);
        if (ic.length) {
          const bin = ic[0].end - ic[0].start || 16384;
          const cov = new Float32Array(Math.ceil(len / bin));
          for (const e of ic) { const b = Math.floor(e.start / bin); if (b < cov.length) cov[b] = e.score; }
          perChrom[chr] = { bin, cov, reads: counts[i], len, exact: false };
        }
      } catch { /* CSI or old index: no overview */ }
    }
  }
  const hotspots = findHotspots(perChrom);
  return { scan, total, perChrom, hotspots };
}

function addSpan(cov, a, b, bin) {
  let i = Math.floor(a / bin);
  while (a < b && i < cov.length) {
    const e = Math.min(b, (i + 1) * bin);
    cov[i] += (e - a) / bin; a = e; i++;
  }
}

/** contiguous regions with reads (≥ 2% of the chromosome's peak), merged across small gaps */
export function findHotspots(perChrom) {
  const out = [];
  for (const [chr, d] of Object.entries(perChrom)) {
    let max = 0; for (const v of d.cov) if (v > max) max = v;
    if (!max) continue;
    const thr = Math.max(d.exact ? 0.5 : 0, max * 0.02);
    const gap = Math.max(1, Math.ceil(1000 / d.bin));
    let i = 0;
    while (i < d.cov.length) {
      if (d.cov[i] <= thr) { i++; continue; }
      let j = i, last = i, sum = 0, pk = 0;
      while (j < d.cov.length && j - last <= gap) { if (d.cov[j] > thr) { last = j; } sum += d.cov[j]; pk = Math.max(pk, d.cov[j]); j++; }
      out.push({ chr, start: i * d.bin, end: Math.min(d.len, (last + 1) * d.bin), peak: pk, mass: sum * d.bin, exact: d.exact });
      i = last + 1;
    }
  }
  out.sort((a, b) => b.mass - a.mass);
  return out;
}

export function overviewToBedGraph(ov, name) {
  const lines = [`track type=bedGraph name="${name}"`];
  for (const [chr, d] of Object.entries(ov.perChrom)) {
    let i = 0;
    while (i < d.cov.length) {
      const v = d.cov[i]; let j = i + 1;
      while (j < d.cov.length && Math.abs(d.cov[j] - v) < 1e-6) j++;
      if (v > 0) lines.push(`${chr}\t${i * d.bin}\t${Math.min(d.len, j * d.bin)}\t${+v.toFixed(2)}`);
      i = j;
    }
  }
  return lines.join('\n') + '\n';
}
