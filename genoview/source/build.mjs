import { build } from 'esbuild';
import { readFileSync, writeFileSync, mkdirSync } from 'fs';
import { igvPatch } from './patch_igv.mjs';

const r = await build({
  entryPoints: ['src/app.js'], bundle: true, minify: true, format: 'iife', write: false,
  target: ['chrome100', 'firefox100', 'safari15'], legalComments: 'none', charset: 'utf8', plugins: [igvPatch],
});
let js = r.outputFiles[0].text.replace(/<\/script/gi, '<\\/script');
const tpl = readFileSync('src/template.html', 'utf8');
const html = tpl.replace('<script>/*__BUNDLE__*/</script>', () => `<script>\n${js}\n</script>`);
mkdirSync('dist', { recursive: true });
writeFileSync('dist/GenoView.html', html);
console.log('dist/GenoView.html', (html.length / 1e6).toFixed(2), 'MB');
