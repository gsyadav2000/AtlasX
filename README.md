# AtlasX

AtlasX is an open-source Python toolkit for genomic peak annotation and downstream analysis.

It is designed to provide a modular framework for working with genomic datasets, including scATAC-seq and other peak-based experiments.

---

## Features

Current features include:

- Load ATAC-seq datasets
- Parse GENCODE GTF annotation files
- Chromosome-based gene indexing
- Efficient nearby gene searching
- Peak-to-gene annotation
- Distance calculation between peaks and genes
- CSV export
- Pandas DataFrame export
- Distance-based annotation filtering

---

## Installation

Install directly from GitHub:

```bash
python -m pip install git+https://github.com/gsyadav2000/AtlasX.git
```

---

## Quick Start

```python
from atlasx.loader.atac_loader import ATACLoader
from atlasx.database.gene_database import GeneDatabase
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.annotation.gene_annotator import GeneAnnotator

loader = ATACLoader(
    "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
)

dataset = loader.load()

database = GeneDatabase(
    "data/reference/gencode.v47.basic.annotation.gtf"
)

genes = database.load()

index = ChromosomeIndex(genes).build()

finder = NearbyGeneFinder(index)

annotator = GeneAnnotator(finder)

annotations = annotator.annotate_dataset(dataset.peaks[:100])

annotator.to_csv(
    annotations,
    "annotations.csv"
)
```

---

## Project Structure

```
AtlasX/
│
├── atlasx/
│   ├── annotation/
│   ├── core/
│   ├── database/
│   ├── loader/
│   └── utils/
│
├── tests/
├── data/
├── README.md
├── pyproject.toml
└── requirements.txt
```

---

## Roadmap

### Completed

- Dataset loading
- GTF parser
- Gene database
- Chromosome indexing
- Nearby gene search
- Peak annotation
- CSV export
- DataFrame export
- Distance filtering
- Installable Python package

### Planned

- BED file support
- Command-line interface
- Gene Ontology analysis
- Pathway analysis
- Peak visualization
- PyPI release

---

## License

MIT License

---

## Author

Ghanshyam Yadav

CSIR-UGC NET JRF

IIIT-Delhi

---

## Citation

If AtlasX contributes to your research, please cite the GitHub repository until a formal publication is available.