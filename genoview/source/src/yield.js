// Yield to the event loop without timers: background tabs throttle setTimeout to >=1 s
// (and to once a minute after a while), which made long analyses look frozen.
// MessageChannel tasks are not throttled. Time-sliced: only actually yields every ~25 ms.
let last = 0;
const ch = typeof MessageChannel !== 'undefined' ? new MessageChannel() : null;
const queue = [];
if (ch) ch.port1.onmessage = () => { const r = queue.shift(); if (r) r(); };
export function yieldNow(force = false) {
  const now = typeof performance !== 'undefined' ? performance.now() : Date.now();
  if (!force && now - last < 25) return Promise.resolve();
  last = now;
  if (!ch) return new Promise(r => setTimeout(r, 0));
  return new Promise(r => { queue.push(r); ch.port2.postMessage(0); });
}
