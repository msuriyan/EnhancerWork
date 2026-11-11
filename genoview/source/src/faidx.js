import { yieldNow } from './yield.js';
// Build a .fai index for a plain (uncompressed) FASTA File in the browser, streaming.
export async function buildFai(file, onProgress = () => {}) {
  const reader = file.stream().getReader();
  const out = [];
  let offset = 0;            // absolute byte offset of current chunk start
  let state = 'start';       // 'header' | 'seq'
  let name = '', nameBytes = [];
  let seqLen = 0, seqStart = 0, lineBases = 0, lineBytes = 0, curLineBases = 0, curLineBytes = 0, linesDone = 0, irregular = false;
  let lastProgress = 0, shortSeen = false;
  const finish = () => {
    if (!name) return;
    if (curLineBases > 0) { closeLine(true); }
    out.push(`${name}\t${seqLen}\t${seqStart}\t${lineBases || seqLen}\t${lineBytes || seqLen + 1}`);
  };
  const closeLine = (last) => {
    if (curLineBases === 0 && curLineBytes <= 1) { curLineBytes = 0; return; }
    if (shortSeen) irregular = true;            // a line after a short line
    if (linesDone === 0) { lineBases = curLineBases; lineBytes = curLineBytes; }
    else if (curLineBases > lineBases) irregular = true;
    else if (curLineBases < lineBases) shortSeen = true;
    linesDone++; curLineBases = 0; curLineBytes = 0;
  };
  for (;;) {
    const { value, done } = await reader.read();
    if (done) break;
    for (let i = 0; i < value.length; i++) {
      const b = value[i];
      if (state === 'header') {
        if (b === 10) {
          name = new TextDecoder().decode(new Uint8Array(nameBytes)).trim().split(/\s+/)[0];
          state = 'seq'; seqStart = offset + i + 1; seqLen = 0; linesDone = 0; shortSeen = false; lineBases = 0; lineBytes = 0; curLineBases = 0; curLineBytes = 0;
        } else nameBytes.push(b);
        continue;
      }
      if (b === 62 && curLineBytes === 0) { // '>' at line start
        if (state === 'seq') {
          if (irregular) throw new Error(`FASTA sequence "${name}" has uneven line lengths; reformat it (e.g. seqkit seq -w 60) or supply a .fai`);
          finish();
        }
        state = 'header'; nameBytes = []; continue;
      }
      if (state !== 'seq') continue;
      curLineBytes++;
      if (b === 10) { closeLine(false); }
      else if (b !== 13) { curLineBases++; seqLen++; }
    }
    offset += value.length;
    if (offset - lastProgress > 50e6) { lastProgress = offset; onProgress(offset / file.size); await yieldNow(true); }
  }
  if (state === 'seq') {
    if (curLineBases > 0) { // final line without newline
      if (linesDone === 0) { lineBases = curLineBases; lineBytes = curLineBases + 1; }
      curLineBases = 0;
    }
    if (irregular) throw new Error(`FASTA sequence "${name}" has uneven line lengths`);
    out.push(`${name}\t${seqLen}\t${seqStart}\t${lineBases || seqLen}\t${lineBytes || seqLen + 1}`);
  }
  onProgress(1);
  return out.join('\n') + '\n';
}
