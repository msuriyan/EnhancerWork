# ATAC Comparative Panels

A single-file, browser-based viewer for making publication figures of chromatin accessibility (ATAC-seq) across several species. You can also use it to find which peaks in a reference species are conserved and open in the others, using a Progressive Cactus whole-genome alignment.

- **One HTML file, nothing to install.** Download it and open it in Chrome, Edge or Firefox.
- **Your data stay on your computer.** Files are read locally by the browser and never uploaded.
- **No external code libraries.** All file readers and drawing are written in plain JavaScript.

## Quick start

1. Download [`atac_comparative_panels.html`](atac_comparative_panels.html): open it and use **Download raw file**.
2. Open it in your browser.
3. Click **Load example** to see a 5-species figure made from synthetic data, or **New figure** to start your own.

You can also run it online from GitHub Pages (see below). Even then, files are only read by your own browser.

## What it does

- **Any number of species**, each with its own genome, region, chromosome and orientation.
- **Per-species data:** genome FASTA (+ `.fai`), gene annotation (GFF3 / GTF / BED12), peaks (BED / narrowPeak / bigBed) and ATAC signal (bigWig / bedGraph), one file per tissue or sample. Every file is optional.
- **Tissue layout:** draw each sample above or below a baseline, overlay samples in one track or give each its own track. Replicates can be averaged.
- **Readable peaks:** robust scaling that ignores outlier spikes, plus gain, square-root/log transforms, smoothing and per-sample boost.
- **Gene models:** all transcripts, or longest / shortest per gene. Search by symbol, LOC/gene ID, transcript ID or description. Hide overlapping genes to keep only the gene of interest.
- **Navigation:** zoom by dragging on the figure, double-click, or the ◀ + − ▶ buttons.
- **Phylogeny:** load a Newick tree, e.g. the Cactus guide tree. It is pruned to your species, and branch lengths are optional.
- **Homology and conserved peaks** from a Cactus **bigMaf** (see below).
- **Marks:** draw boxes over elements and link them across species with ribbons. Color, opacity, border and height are adjustable.
- **Publication control:** width in mm, sizes in points, and font, size, color and visibility for every label. Editable legend and free text boxes.
- **Export** to SVG (fully editable in Illustrator or Inkscape) or PNG at 300, 600 or 1200 dpi.
- **Sessions:** save and reopen a figure. In Chrome and Edge, a session's data files reopen automatically; elsewhere, pick the data folder once.

## Finding conserved peaks with a Cactus alignment

You need:

- A **bigMaf** (`.bb`) exported from your HAL file, with your reference species as the reference genome.
- A **peak file** for the reference species, and for each species you want to test.

In the tool, open **7 · Homology**:

1. Load the bigMaf. Species are matched to the alignment's genome names automatically: `Jcoe`, `J_coenia` and `Junonia_coenia` all match *Junonia coenia*. You can correct any match.
2. Choose the reference species and which peaks to test: the current window, chromosome or whole genome.
3. Click **Find in other species**. For every reference peak and species you get one call:
   - **Open:** the peak aligns, and that species has a peak within 500 bp.
   - **Aligned, no peak:** the sequence is there but not accessible; the nearest peak's distance is reported.
   - **Partly aligned / not aligned.**
4. Click any peak to move every species to its homologous region. Each species gets a short explanation of its call.
5. Export a per-peak table or a table of all reference × homologous peak pairs (TSV).

Turn on **Follow reference** to move all species to the homologous region whenever you zoom or pan the reference.

### Making a bigMaf from a HAL file

```bash
hal2maf --refGenome Jcoe --noAncestors alignment.hal Jcoe.maf
mafToBigMaf Jcoe Jcoe.maf stdout | sort -k1,1 -k2,2n > Jcoe.bigMaf.txt
# bigMaf.as is available from the UCSC Genome Browser documentation
bedToBigBed -type=bed3+1 -as=bigMaf.as -tab Jcoe.bigMaf.txt Jcoe.chrom.sizes Jcoe.bigMaf.bb
```

A bigMaf only works from its reference genome. To start from another species, make a bigMaf with that species as reference.

### RefSeq vs GenBank names

If a species' files use RefSeq names (`NC_…`/`NW_…`) and the alignment uses GenBank names (`CM…`), or the reverse, load that assembly's NCBI `*_assembly_report.txt` under **1 · Species & data → Chromosome name aliases**. Both names are then treated as the same sequence everywhere. This applies to the same assembly only; data on a different assembly version must be lifted over first.

## Supported files

| Data | Formats |
| --- | --- |
| ATAC signal | bigWig, bedGraph (`.gz` ok) |
| Peaks | BED, narrowPeak, broadPeak (`.gz` ok), bigBed |
| Genes | GFF3, GTF, BED12, bigBed12 |
| Genome | FASTA, optional `.fai` |
| Alignment | Cactus bigMaf (`.bb`) |
| Tree | Newick |
| Name aliases | NCBI assembly report, or any table of names |

bigWig, bigBed and FASTA files are read in pieces, so multi-GB files open quickly.

## Browser support

Chrome and Edge support every feature. Firefox and Safari support everything except automatically reopening a saved session's files; there you choose the data folder once instead.

## Running it on GitHub Pages

In this repository: **Settings → Pages → Deploy from a branch → `main` / root**. The tool is then served at `https://<user>.github.io/atac-comparative-panels/`. Files you load are still only read by your own browser.

## How it works (short)

- bigWig/bigBed headers, chromosome B+ trees and R-tree indexes are parsed directly. Only the data blocks for the current window are fetched and decompressed (`DecompressionStream`).
- Signal is summarised per display bin (maximum by default), using bigWig zoom levels where they are fine enough.
- Each reference peak is projected base by base through the aligned columns of the bigMaf blocks. Projected positions are split into copies at gaps larger than max(2 kb, 2 × peak length), and copies with at least 20% of the peak aligned are tested for nearby peaks.
- The figure is drawn as SVG in points; PNG is rendered from the same SVG.

## Citing

If you use this tool, please cite this repository. Please also cite the tools it builds on:

- Armstrong J. et al. Progressive Cactus is a multiple-genome aligner for the thousand-genome era. *Nature* 587, 246–251 (2020). https://doi.org/10.1038/s41586-020-2871-y
- Hickey G. et al. HAL: a hierarchical format for storing and analyzing multiple genome alignments. *Bioinformatics* 29, 1341–1342 (2013). https://doi.org/10.1093/bioinformatics/btt128
- Kent W.J. et al. BigWig and BigBed: enabling browsing of large distributed datasets. *Bioinformatics* 26, 2204–2207 (2010). https://doi.org/10.1093/bioinformatics/btq351

## License

MIT, see [LICENSE](LICENSE).
