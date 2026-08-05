# EpiMatch

EpiMatch is an open-source Python toolkit for analyzing single-cell ATAC-seq data: scoring which genes are statistically enriched near each cell's accessible chromatin, clustering cells by that enrichment signal, and matching cells against reference profiles with a proper statistical null model.

It's built as a from-scratch, transparent reimplementation of the core ideas behind [scEpiSearch](https://genome.cshlp.org/content/33/2/218) (Mishra et al., *Genome Research* 2023) — gene-enrichment scoring in place of raw peak comparison, and null-model-calibrated matching in place of fixed thresholds — validated step by step against real public data, and checked directly against scEpiSearch's own source code, rather than assumed to work.

---

## What it actually does right now

**Gene enrichment scoring.** For each cell, ranks which genes are statistically enriched near that cell's most accessible peaks (TF-IDF-weighted peak selection, binomial significance testing with an exact hypergeometric fallback for small samples).

**Cell clustering.** Groups cells by their enrichment profiles — either Jaccard similarity over top-enriched gene sets, or Spearman correlation over the full continuous enrichment vector (scEpiSearch's own matching approach) — then hierarchical clustering, tested against real marker gene panels.

**Marker gene enrichment.** Tests whether a cluster's characteristic genes are statistically enriched for known cell-type markers, using panels sourced from the [Seurat PBMC3k reference tutorial](https://satijalab.org/seurat/archive/v4.3/pbmc3k_tutorial), extended with a small number of markers observed directly in this project's own validated reference profiles (documented in `scoring/marker_panels.py`).

**Reference matching.** Matches a query cell against reference cluster profiles using a synthetic random-gene null model, so p-values aren't biased by how common a reference's cell type happens to be in the training data. Validated both within one dataset (held-out cells) and across two independent public datasets.

**Auto-updating ingestion pipeline.** Searches NCBI GEO for candidate single-cell datasets, downloads and safely extracts them (path-traversal and size guarded), tries known loaders, and runs automated genome-build and quality checks before anything enters the reference pool — every decision logged to a versioned manifest.

**Peak annotation and file I/O.** Loads 10x-style HDF5 (including multiome/feature-barcode matrices) and BED/BED.gz files, parses GENCODE GTF annotation, maps peaks to nearby genes by chromosome-indexed binary search.

### Validated results

- **Unsupervised clustering** (10x Genomics 10k PBMC dataset, full 8,728 cells), with no cell-type labels given, correctly separates T cells, two monocyte subtypes (CD14+ and FCGR3A+), NK cells, B cells, and plasmacytoid dendritic cells — each confirmed against real marker genes, several found directly in the data during validation rather than assumed from the reference panel alone.
- **Reference matching, held-out evaluation** (same dataset, train/holdout split): 82.1% accuracy on held-out cells, 36 points above the 46.1% majority-class baseline, with match confidence tracking correctness (98% of matches statistically significant).
- **Cross-dataset matching**: reference profiles built from the 10k-PBMC dataset, tested against held-out cells from an independent 5k-PBMC dataset (different donor). NK and T-cell signal replicated independently; monocyte replicated after expanding the marker panel with genes observed directly in the reference profile.
- **GEO ingestion pipeline**: tested live against 15 real accessions discovered by keyword search; 3 genuinely usable datasets ingested end-to-end (real peak data, correct genome build, passing QC), the rest correctly and specifically rejected (unsupported format, oversized download, or confirmed non-ATAC data) rather than silently failing.
- **scEpiSearch parameter comparison**: obtained and read scEpiSearch's actual source code, confirmed its real parameter values (1,000,000 bp gene-search window, 10,000 peaks per cell, Spearman correlation matching) differ substantially from this project's own defaults. Under a matched, controlled comparison — same cells, same matching method, only these parameters changed — EpiMatch's own parameter values produced a statistically significant dominant cluster where scEpiSearch's confirmed values did not, with the mechanism (dilution toward generic, less cell-type-specific genes) measured directly rather than inferred.

### Honest scope

- The scEpiSearch parameter comparison above used a 3,000-cell subsample, not the full dataset, due to memory constraints on the development machine — disclosed, not hidden.
- Cross-dataset matching so far uses two related PBMC datasets, not scEpiSearch's full external reference pool (millions of cells).
- scEpiSearch's actual, complete pipeline (their real null model built from merged real-cell queries, their full reference atlas) has not been run end-to-end for a head-to-head benchmark — their code requires a legacy Python/R environment and a 34GB reference data download, neither attempted yet.
- The ingestion pipeline currently recognizes HDF5 and BED formats only; a meaningful fraction of real GEO submissions use other formats (R objects, sparse matrix triplets, bigWig) not yet supported.

---

## Installation

```bash
python -m pip install git+https://github.com/gsyadav2000/EpiMatch.git
```

Or clone and install in editable mode for development:

```bash
git clone https://github.com/gsyadav2000/EpiMatch.git
cd EpiMatch
python -m pip install -e .
```

You'll also need a GENCODE annotation matching your data's genome build (see note below on build-matching) — download from [gencodegenes.org](https://www.gencodegenes.org/human/).

---

## Quick start

The full worked pipeline — load data, score gene enrichment, cluster cells, test against marker panels — lives in `examples/cell_clustering_demo.py`. Reference matching is in `examples/reference_matching_demo.py`. Cross-dataset matching is in `examples/build_reference_from_dataset_a.py` and `examples/match_dataset_b_against_saved_reference.py`. The GEO ingestion pipeline is in `examples/run_ingestion_pipeline_demo.py`. A minimal enrichment example:

```python
from epimatch.loader.atac_loader import ATACLoader
from epimatch.database.gene_database import GeneDatabase
from epimatch.database.chromosome_index import ChromosomeIndex
from epimatch.database.nearby_gene_finder import NearbyGeneFinder
from epimatch.scoring.gene_enrichment import GeneEnrichmentScorer

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
epimatch summary path/to/dataset.h5
epimatch peaks path/to/dataset.h5 -n 20
epimatch genes path/to/annotation.gtf -n 20
epimatch version
```

---

## Project structure
EpiMatch/
│
├── epimatch/
│ ├── core/ # Peak, Gene, Genome, Dataset objects
│ ├── loader/ # HDF5 (10x-style, incl. multiome) dataset loading
│ ├── io/ # BED file reading
│ ├── database/ # GTF parsing, chromosome indexing, nearby-gene search
│ ├── annotation/ # Peak-to-gene annotation, CSV/DataFrame export
│ ├── scoring/ # Gene enrichment, clustering, marker panels, reference matching
│ ├── ingestion/ # GEO discovery, download, QC, ingestion pipeline
│ └── cli/ # Command-line interface
│
├── examples/ # Runnable demo scripts for every module
├── tests/ # pytest test suite (fast, no large data dependencies)
├── data/
│ ├── example/ # Small tracked example files
│ ├── raw/ # Large datasets (gitignored)
│ ├── processed/ # Gitignored
│ ├── reference/ # GTF annotation files (gitignored)
│ └── manifests/ # Ingestion pipeline manifest (tracked - small, human-readable)
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

- HDF5 (including multiome/feature-barcode) and BED dataset loading
- GTF parsing with protein-coding filtering
- Chromosome-indexed nearby-gene search
- Peak-to-gene annotation, CSV/DataFrame export
- Command-line interface
- Gene enrichment scoring (TF-IDF weighting, binomial significance with exact fallback)
- Batch enrichment across cell populations
- Marker gene set enrichment (Seurat PBMC3k-sourced panels, extended with validated additions)
- Cell-to-cell similarity and hierarchical clustering (Jaccard and Spearman-correlation methods)
- Reference matching with bias-free synthetic null model, validated on held-out data
- Cross-dataset reference matching, validated against an independent public dataset
- Automatically updating reference database (GEO discovery, download, safe extraction, QC gating, versioned manifest) — tested live against 15 real accessions
- Direct comparison against scEpiSearch's real source code and confirmed parameter values

### Planned

- Full-scale (non-subsampled) run of the scEpiSearch parameter comparison, pending more available memory
- Head-to-head benchmark against scEpiSearch's actual, complete pipeline (their real reference pool and null model)
- Broader ingestion format support (sparse matrix triplets, R objects, bigWig)
- Approximate nearest-neighbor search for scaling matching to large reference pools
- Scheduled/automatic execution of the ingestion pipeline (currently run manually)
- PyPI release

---

## License

MIT License

---

## Author

Ghanshyam Yadav

CSIR-NET JRF (AIR 68)

PhD, IIIT Delhi, under the guidance of Dr. Vibhor Kumar

---

## Citation

If EpiMatch contributes to your research, please cite the GitHub repository until a formal publication is available.