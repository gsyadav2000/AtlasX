# AtlasX

AtlasX is an open-source Python toolkit for analyzing single-cell ATAC-seq data: scoring which genes are statistically enriched near each cell's accessible chromatin, clustering cells by that enrichment signal, and matching cells against reference profiles with a proper statistical null model.

It's built as a from-scratch, transparent reimplementation of the core ideas behind [scEpiSearch](https://genome.cshlp.org/content/33/2/218) (Mishra et al., *Genome Research* 2023) — gene-enrichment scoring in place of raw peak comparison, and null-model-calibrated matching in place of fixed thresholds — validated step by step against real PBMC data rather than assumed to work.

---

## What it actually does right now

**Gene enrichment scoring.** For each cell, ranks which genes are statistically enriched near that cell's most accessible peaks (TF-IDF-weighted peak selection, binomial/hypergeometric significance testing with automatic exact-test fallback for small samples).

**Cell clustering.** Groups cells by the overlap in their enriched-gene sets (vectorized Jaccard similarity, hierarchical clustering), then tests each cluster against real marker gene panels.

**Marker gene enrichment.** Tests whether a cluster's characteristic genes are statistically enriched for known cell-type markers, using panels sourced from the [Seurat PBMC3k reference tutorial](https://satijalab.org/seurat/archive/v4.3/pbmc3k_tutorial) rather than an ad hoc gene list.

**Reference matching.** Matches a query cell against reference cluster profiles using a synthetic random-gene null model, so p-values aren't biased by how common a reference's cell type happens to be in the training data.

**Peak annotation and file I/O.** Loads 10x-style HDF5 and BED files, parses GENCODE GTF annotation, maps peaks to nearby genes by chromosome-indexed binary search.

### Validated results (10x Genomics 10k PBMC dataset, all 8,728 cells)

- Unsupervised clustering, with no cell-type labels given, separates T cells (p = 1.1e-4), CD14+ monocytes and FCGR3A+ monocytes as distinct subclusters (p = 0.024 and 0.024), and NK cells (p = 0.012) — each confirmed against the Seurat PBMC3k marker panel.
- Reference matching, evaluated on cells held out from reference-building, correctly identifies a cell's true cluster 82.1% of the time (36 points above the 46.1% majority-class baseline), with match confidence tracking correctness (98% of matches statistically significant, 81.9% both correct and significant).

### Important limitation

Reference matching so far is **self-referential**: reference profiles are built from a training split of the *same* dataset being queried, not from an independent external atlas. This was the right first step — it validates that the matching architecture (synthetic null model, held-out evaluation, no circularity) works correctly before investing in sourcing real cross-dataset reference data — but it is not yet the cross-dataset, cross-study matching that is scEpiSearch's actual contribution. That remains the next major piece of work.

---

## Installation

```bash
python -m pip install git+https://github.com/gsyadav2000/AtlasX.git
```

Or clone and install in editable mode for development:

```bash
git clone https://github.com/gsyadav2000/AtlasX.git
cd AtlasX
python -m pip install -e .
```

You'll also need a GENCODE annotation matching your data's genome build (see note below on build-matching) — download from [gencodegenes.org](https://www.gencodegenes.org/human/).

---

## Quick start

The full worked pipeline — load data, score gene enrichment, cluster cells, test against marker panels — lives in `examples/cell_clustering_demo.py`. Reference matching is in `examples/reference_matching_demo.py`. A minimal enrichment example:

```python
from atlasx.loader.atac_loader import ATACLoader
from atlasx.database.gene_database import GeneDatabase
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.scoring.gene_enrichment import GeneEnrichmentScorer

dataset = ATACLoader("data/raw/your_dataset.h5").load()

genes = GeneDatabase(
    "data/reference/your_gencode_annotation.gtf",
    protein_coding_only=True
).load()

index = ChromosomeIndex(genes).build()
finder = NearbyGeneFinder(index)

scorer = GeneEnrichmentScorer(dataset.peaks, finder, dataset.matrix)

top_genes = scorer.enrich(cell_index=0, top_n=2000, top_genes=20)

for gene, p_value, p_adjusted in top_genes:
    print(f"{gene:15} p = {p_value:.2e}")
```

**Important: match your genome build.** Peak coordinates and your GTF annotation must be in the same genome build (both hg19/GRCh37, or both hg38/GRCh38). A mismatch produces plausible-looking but biologically meaningless gene assignments with no error or warning — this exact bug cost significant debugging time during development. Older 10x Genomics scATAC-seq datasets (the "v1" chemistry era) are typically hg19; GENCODE's default releases are hg38, so you likely want a [GRCh37-mapped GENCODE release](https://www.gencodegenes.org/human/grch37_mapped_releases.html) (filename contains `lift37`) unless you know your peaks are hg38.

---

## Command-line interface

```bash
atlasx summary path/to/dataset.h5
atlasx peaks path/to/dataset.h5 -n 20
atlasx genes path/to/annotation.gtf -n 20
atlasx version
```

---

## Project structure
AtlasX/
│
├── atlasx/
│ ├── core/ # Peak, Gene, Genome, Dataset objects
│ ├── loader/ # HDF5 (10x-style) dataset loading
│ ├── io/ # BED file reading
│ ├── database/ # GTF parsing, chromosome indexing, nearby-gene search
│ ├── annotation/ # Peak-to-gene annotation, CSV/DataFrame export
│ ├── scoring/ # Gene enrichment, clustering, marker panels, reference matching
│ └── cli/ # Command-line interface
│
├── examples/ # Runnable demo scripts for every module
├── tests/ # pytest test suite (fast, no large data dependencies)
├── data/
│ ├── example/ # Small tracked example files
│ ├── raw/ # Large datasets (gitignored)
│ ├── processed/ # Gitignored
│ └── reference/ # GTF annotation files (gitignored)
│
├── README.md
├── pyproject.toml
└── requirements.txt
---

## Testing

```bash
python -m pytest
```

The test suite is fast and self-contained — no large dataset downloads required. Demo scripts that need real data live in `examples/`, separate from `tests/`, so `pytest` never depends on files that aren't in the repository.

---

## Roadmap

### Completed

- HDF5 and BED dataset loading
- GTF parsing with protein-coding filtering
- Chromosome-indexed nearby-gene search
- Peak-to-gene annotation, CSV/DataFrame export
- Command-line interface
- Gene enrichment scoring (TF-IDF weighting, binomial/hypergeometric significance)
- Batch enrichment across cell populations
- Marker gene set enrichment (Seurat PBMC3k-sourced panels)
- Cell-to-cell similarity and hierarchical clustering
- Reference matching with bias-free synthetic null model, validated on held-out data

### Planned

- Cross-dataset reference matching against an independent external atlas
- Automatically updating reference database (scheduled ingestion from public repositories, with QC gating before any dataset is added)
- Approximate nearest-neighbor search for scaling matching to large reference pools
- PyPI release

---

## License

MIT License

---

## Author

Ghanshyam Yadav

CSIR-NET JRF (AIR 68)

---

## Citation

If AtlasX contributes to your research, please cite the GitHub repository until a formal publication is available.