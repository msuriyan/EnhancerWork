import { yieldNow } from './yield.js';
// A BAM without an index: read every record once into memory, then answer range queries
// with the same interface GenoView uses from @gmod/bam's BamFile.
import { streamBamRecords } from '@gmod/bam';
import { BlobFile } from 'generic-filehandle2';

export class InMemoryBam {
  constructor(file) { this.file = file; this.byRef = []; this.maxLen = []; this.header = ''; this.indexToChr = []; this.loaded = null; }

  load(onProgress = () => {}) {
    if (this.loaded) return this.loaded;
    this.loaded = (async () => {
      let n = 0;
      for await (const batch of streamBamRecords({
        bamFilehandle: new BlobFile(this.file),
        onHeader: h => { this.header = h.headerText; this.indexToChr = h.indexToChr; this.byRef = h.indexToChr.map(() => []); this.maxLen = h.indexToChr.map(() => 0); },
      })) {
        for (const r of batch) {
          const id = r.ref_id; if (id < 0 || !this.byRef[id]) continue;
          this.byRef[id].push(r);
          const L = r.end - r.start; if (L > this.maxLen[id]) this.maxLen[id] = L;
        }
        n += batch.length; onProgress(n);
        await yieldNow();
      }
      for (const arr of this.byRef) {
        let sorted = true; for (let i = 1; i < arr.length; i++) if (arr[i].start < arr[i - 1].start) { sorted = false; break; }
        if (!sorted) arr.sort((a, b) => a.start - b.start);
      }
      this.chrToIndex = Object.fromEntries(this.indexToChr.map((c, i) => [c.refName, i]));
      this.index = { lineCount: async i => this.byRef[i]?.length || 0, indexCov: async () => [] };
      return this;
    })();
    return this.loaded;
  }

  async getHeader() { await this.load(); return this.header; }

  async getRecordsForRange(chr, start, end) {
    await this.load();
    const id = this.chrToIndex[chr]; if (id == null) return [];
    const arr = this.byRef[id];
    const from = start - this.maxLen[id];
    let lo = 0, hi = arr.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (arr[mid].start < from) lo = mid + 1; else hi = mid; }
    const out = [];
    for (let i = lo; i < arr.length && arr[i].start < end; i++) if (arr[i].end > start) out.push(arr[i]);
    return out;
  }
}
