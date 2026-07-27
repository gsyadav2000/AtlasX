"""
AtlasX Command Line Interface
"""

import argparse

from atlasx.loader.atac_loader import ATACLoader
from atlasx.database.gene_database import GeneDatabase

VERSION = "0.14.0"


def dataset_summary(args):
    """Display a summary of an ATAC dataset."""

    loader = ATACLoader(args.file)
    dataset = loader.load()

    print()
    print("=" * 50)
    print("AtlasX Dataset Summary")
    print("=" * 50)
    print(f"Cells : {dataset.n_cells:,}")
    print(f"Peaks : {dataset.n_peaks:,}")
    print(f"Shape : {dataset.matrix.shape}")
    print("=" * 50)


def show_version(args):
    """Display AtlasX version."""

    print()
    print("=" * 40)
    print("AtlasX")
    print("=" * 40)
    print(f"Version : {VERSION}")
    print("Author  : Ghanshyam Yadav")
    print("=" * 40)


def show_peaks(args):
    """Display the first N peaks."""

    loader = ATACLoader(args.file)
    dataset = loader.load()

    print()
    print("=" * 50)
    print(f"First {args.number} Peaks")
    print("=" * 50)

    limit = min(args.number, len(dataset.peaks))

    for i, peak in enumerate(dataset.peaks[:limit], start=1):
        print(
            f"{i:>3}. "
            f"{peak.chromosome}:{peak.start}-{peak.end}"
        )

    print("=" * 50)


def show_genes(args):
    """Display the first N genes."""

    db = GeneDatabase(args.file)
    genes = db.load()

    print()
    print("=" * 50)
    print(f"First {args.number} Genes")
    print("=" * 50)

    limit = min(args.number, len(genes))

    for i, gene in enumerate(genes[:limit], start=1):

        print(f"{i:>3}. {gene.name}")
        print(
            f"     {gene.chromosome}:{gene.start}-{gene.end} ({gene.strand})"
        )
        print()

    print("=" * 50)


def main():

    parser = argparse.ArgumentParser(
        prog="atlasx",
        description="AtlasX : Single-cell ATAC Analysis Toolkit"
    )

    subparsers = parser.add_subparsers(dest="command")

    # ---------------- Summary ---------------- #

    summary_parser = subparsers.add_parser(
        "summary",
        help="Display dataset summary"
    )

    summary_parser.add_argument(
        "file",
        help="Path to HDF5 dataset"
    )

    summary_parser.set_defaults(
        func=dataset_summary
    )

    # ---------------- Version ---------------- #

    version_parser = subparsers.add_parser(
        "version",
        help="Display AtlasX version"
    )

    version_parser.set_defaults(
        func=show_version
    )

    # ---------------- Peaks ---------------- #

    peaks_parser = subparsers.add_parser(
        "peaks",
        help="Display the first N peaks"
    )

    peaks_parser.add_argument(
        "file",
        help="Path to HDF5 dataset"
    )

    peaks_parser.add_argument(
        "-n",
        "--number",
        type=int,
        default=10,
        help="Number of peaks to display"
    )

    peaks_parser.set_defaults(
        func=show_peaks
    )

    # ---------------- Genes ---------------- #

    genes_parser = subparsers.add_parser(
        "genes",
        help="Display the first N genes"
    )

    genes_parser.add_argument(
        "file",
        help="Path to GTF annotation"
    )

    genes_parser.add_argument(
        "-n",
        "--number",
        type=int,
        default=10,
        help="Number of genes to display"
    )

    genes_parser.set_defaults(
        func=show_genes
    )

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()