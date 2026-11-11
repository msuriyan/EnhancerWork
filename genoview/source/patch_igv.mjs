// esbuild plugin: speed up igv.js CoverageMap.incCounts for long reads
// (skip negative offsets that turn the coverage array into a dictionary; avoid dynamic string keys)
export const igvPatch = {
  name: 'igv-patch',
  setup(build) {
    build.onLoad({ filter: /igv\.esm\.js$/ }, async args => {
      const fs = await import('fs');
      let src = fs.readFileSync(args.path, 'utf8');
      const start = src.indexOf('        function incBlockCount(block) {');
      const end = src.indexOf('    getPosCount(pos, base) {', start);
      if (start < 0 || end < 0) throw new Error('igv patch anchor not found');
      const fast = `        function incBlockCount(block) {
            if ('S' === block.type) return
            const seq = alignment.seq, qual = alignment.qual, seqOffset = block.seqOffset;
            const pos = alignment.strand, cov = self.coverage, thr = self.threshold;
            const hasSeq = seq != undefined && seq !== '*';
            let max = self.maximum;
            let i = block.start - self.bpStart, j = 0;
            if (i < 0) { j = -i; i = 0; }
            for (; j < block.len; i++, j++) {
                let c = cov[i];
                if (!c) c = cov[i] = new Coverage(thr);
                const k = seqOffset + j;
                const q = qual && k < qual.length ? qual[k] : 30;
                const code = hasSeq ? seq.charCodeAt(k) : 78;
                switch (code) {
                    case 65: case 97: if (pos) c.posA++; else c.negA++; c.qualA += q; break;
                    case 67: case 99: if (pos) c.posC++; else c.negC++; c.qualC += q; break;
                    case 71: case 103: if (pos) c.posG++; else c.negG++; c.qualG += q; break;
                    case 84: case 116: if (pos) c.posT++; else c.negT++; c.qualT += q; break;
                    default: if (pos) c.posN++; else c.negN++; c.qualN += q;
                }
                c.total++; c.qual += q;
                if (c.total > max) max = c.total;
            }
            self.maximum = max;
        }
    }

`;
      src = src.slice(0, start) + fast + src.slice(end);
      // --- hide indels / mismatches in reads by size or user regions (GenoView hooks) ---
      const rep = (a, b) => { if (!src.includes(a)) throw new Error('igv patch anchor missing: ' + a.slice(0, 60)); src = src.replace(a, b); };
      rep(`                const yStrokedLine = y + alignmentHeight / 2;
                for (let gap of alignment.gaps) {`, `                const yStrokedLine = y + alignmentHeight / 2;
                for (let gap of alignment.gaps) {
                    if (gap.type === 'D' && globalThis.__gvHide && globalThis.__gvHide(alignment.chr, gap.start, gap.len, 'D')) {
                        IGVGraphics.fillRect(ctx, (gap.start - bpStart) / bpPerPixel, y, Math.max(1, gap.len / bpPerPixel), alignmentHeight, {fillStyle: alignmentColor});
                        continue
                    }`);
      rep(`                for (let insertionBlock of alignment.insertions) {
                    if (this.hideSmallIndels && insertionBlock.len <= this.indelSizeThreshold) {`, `                for (let insertionBlock of alignment.insertions) {
                    if (globalThis.__gvHide && globalThis.__gvHide(alignment.chr, insertionBlock.start, insertionBlock.len, 'I')) continue
                    if (this.hideSmallIndels && insertionBlock.len <= this.indelSizeThreshold) {`);
      rep(`                        let readChar = seq ? seq.charAt(seqOffset + i) : '';
                        const refChar = offsetBP + i >= 0 ? referenceSequence.charAt(offsetBP + i) : '';`, `                        if (!isSoftClip && !showAllBases && globalThis.__gvHideMM && globalThis.__gvHideMM(alignment.chr, block.start + i)) continue
                        let readChar = seq ? seq.charAt(seqOffset + i) : '';
                        const refChar = offsetBP + i >= 0 ? referenceSequence.charAt(offsetBP + i) : '';`);
      // --- region-wide read sorting: largest/smallest deletion, largest indel, random (GenoView) ---
      rep(`    getSortValue({ position, option, tag, sortAsPairs }, alignmentContainer) {

        if (!option) option = "BASE";
`, `    getSortValue({ position, option, tag, sortAsPairs, start, end }, alignmentContainer) {

        if (!option) option = "BASE";
        if (option === "GV_DEL_MAX" || option === "GV_DEL_MIN" || option === "GV_INDEL_MAX" || option === "GV_RANDOM") {
            if (option === "GV_RANDOM") return this.__rnd ?? (this.__rnd = Math.random())
            let best = 0, any = false;
            for (const a of this.alignments) {
                if (a.start >= end || a.start + a.lengthOnRef <= start) continue
                if (a.gaps) for (const g of a.gaps) if (g.type === 'D' && g.start < end && g.start + g.len > start) { best = Math.max(best, g.len); any = true; }
                if (option === "GV_INDEL_MAX" && a.insertions) for (const ins of a.insertions) if (ins.start >= start && ins.start <= end) { best = Math.max(best, ins.len); any = true; }
                // split read on the same chromosome/strand: the unaligned gap is a large deletion
                try {
                    const sa = a.getTag && a.getTag('SA');
                    if (sa) for (const part of String(sa).split(';')) {
                        const f = part.split(','); if (f.length < 4 || f[0] !== a.chr) continue
                        const ps = +f[1] - 1; let rl = 0; for (const m of f[3].matchAll(/(\\d+)([MDN=X])/g)) rl += +m[1];
                        const ae = a.start + a.lengthOnRef;
                        const gap = ps >= ae ? ps - ae : (ps + rl <= a.start ? a.start - (ps + rl) : 0);
                        const g0 = ps >= ae ? ae : ps + rl, g1 = g0 + gap;
                        if (gap > 0 && g0 < end && g1 > start) { best = Math.max(best, gap); any = true; }
                    }
                } catch (e) {}
            }
            if (option === "GV_DEL_MIN") return any ? best : Number.MAX_VALUE
            return any ? -best : 0
        }
`);
      rep(`        for (let row of this.rows) {
            const alignment = row.findAlignment(options.position, options.sortAsPairs);`, `        for (let row of this.rows) {
            const alignment = (options.start !== undefined && String(options.option).startsWith('GV_'))
                ? row.alignments.find(a => a.start < options.end && a.start + a.lengthOnRef > options.start)
                : row.findAlignment(options.position, options.sortAsPairs);`);
      // --- per-track: hide feature names / strand arrows (GenoView) ---
      rep(`options.drawLabel = options.labelAllFeatures || pixelsPerFeature > 10;`, `options.drawLabel = !this.gvHideNames && (options.labelAllFeatures || pixelsPerFeature > 10);`);
      rep(`                options.drawLabel = true;
                this.render.call(this, feature, bpStart, bpPerPixel, pixelHeight, context, options);`, `                options.drawLabel = !this.gvHideNames;
                this.render.call(this, feature, bpStart, bpPerPixel, pixelHeight, context, options);`);
      rep(`const direction = feature.strand === '+' ? 1 : feature.strand === '-' ? -1 : 0;`, `const direction = this.gvHideArrows ? 0 : feature.strand === '+' ? 1 : feature.strand === '-' ? -1 : 0;`);
      rep(`            for (let x = xLeft; x < xRight; x += step) {`, `            for (let x = xLeft; direction !== 0 && x < xRight; x += step) {`);
      // --- ideogram red region box can be switched off; hidden tracks are skipped in SVG export (GenoView) ---
      rep(`            // Draw red box
            context.strokeStyle = "red";`, `            // Draw red box
            context.strokeStyle = globalThis.__gvNoIdeoBox ? "rgba(0,0,0,0)" : "red";`);
      rep(`    renderSVGContext(context, {deltaX, deltaY}) {

        renderSVGAxis(context, this.track, this.axisCanvas, deltaX, deltaY);`, `    renderSVGContext(context, {deltaX, deltaY}) {

        if (this.track && this.track.gvHidden) return
        renderSVGAxis(context, this.track, this.axisCanvas, deltaX, deltaY);`);
      return { contents: src, loader: 'js' };
    });
  },
};
