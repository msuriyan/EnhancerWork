# EnhancerWork

Four single-file, browser-based tools for making **publication figures** from comparative and functional genomics data. Each tool is one HTML file. There is nothing to install and no server, and your data never leave your computer.

**▶ Open the tools online:** https://msuriyan.github.io/EnhancerWork/

**⬇ Download for offline use:** [latest release](https://github.com/msuriyan/EnhancerWork/releases/latest). Save any `.html` file and double-click it.

| # | Tool | What it does | Open |
|---|------|--------------|------|
| 1 | [Cactus Alignment Viewer](cactus-alignment-viewer/) | Cactus whole-genome alignments (bigMaf/MAF) shown with the phylogeny, ATAC-seq signal and peaks, BED and GFF3/GTF tracks | [launch](https://msuriyan.github.io/EnhancerWork/cactus-alignment-viewer/) |
| 2 | [ATAC Comparative Panels](atac-comparative-panels/) | Multi-species ATAC-seq panels with tissue overlays, and homologous-peak tracing through a Cactus alignment | [launch](https://msuriyan.github.io/EnhancerWork/atac-comparative-panels/atac_comparative_panels.html) |
| 3 | [TFBS Enhancer Map](tfbs-enhancer-map/) | Strand-aware, TRANSFAC-style maps of predicted TF binding sites on enhancers/CREs | [launch](https://msuriyan.github.io/EnhancerWork/tfbs-enhancer-map/) |
| 4 | [GenoView](genoview/) | IGV/UCSC-style browser for BAM/bigWig/BED with indel/SV detection and CRISPR crispant edit analysis | [launch](https://msuriyan.github.io/EnhancerWork/genoview/) |

Each folder has its own README with full instructions and example data.

## Running the tools

- **Online:** use the launch links above. The page is served by GitHub, but the files you open are read only by your own browser.
- **Offline:** download an HTML file from the latest release (or clone this repository) and open it in Chrome, Edge, Firefox or Safari.

## Repository layout

```
index.html                  landing page (GitHub Pages)
cactus-alignment-viewer/    Tool 1 (+ tests/)
atac-comparative-panels/    Tool 2
tfbs-enhancer-map/          Tool 3 (+ example/)
genoview/                   Tool 4 (+ example/, scripts/ Python helpers, source/ build files)
```

## For maintainers

- **GitHub Pages:** *Settings → Pages → Deploy from a branch → `main` / `(root)`*.
- **New release:** *Actions → release → Run workflow*, or push a tag such as `v1.0.0`. The workflow attaches all four tools as standalone HTML files.

## License

MIT. See [LICENSE](LICENSE).
