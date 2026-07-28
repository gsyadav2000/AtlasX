"""
AtlasX Marker Gene Panels

Marker genes are taken from the official Seurat PBMC3k guided
clustering tutorial (Satija Lab), the standard reference marker set
used throughout the single-cell field for identifying PBMC cell
types from clustering results:
https://satijalab.org/seurat/archive/v4.3/pbmc3k_tutorial

This replaces an earlier version of this file that used a larger,
hand-curated list assembled from general immunology knowledge rather
than a single traceable source. The Seurat panels are smaller
per-lineage (1-3 genes each, versus the earlier list's 4-17), which
is an honest tradeoff: less statistical power per test (a smaller
panel needs a larger fraction of it to show up before a hypergeometric
test calls it significant), but every gene's inclusion is citable to
a specific, authoritative, PBMC-specific source rather than to this
project's own judgment.

Two lineages new to this project - Dendritic Cells and
Platelets/Megakaryocytes - are included here because they're part of
the Seurat reference set, even though neither has been tested for in
this project before. Whether either shows up as its own cluster in
this dataset is itself informative, not assumed.
"""

CD4_T_CELL_MARKERS = frozenset({"IL7R", "CCR7"})
CD8_T_CELL_MARKERS = frozenset({"CD8A"})

# Combined T cell panel: this project has been testing at the
# lineage level (T cell vs monocyte vs NK, etc.), not the CD4/CD8
# subtype level, so CD4 and CD8 markers are combined here for
# consistency with earlier results. All three genes are still from
# the same Seurat source, just regrouped.
T_CELL_MARKERS = CD4_T_CELL_MARKERS | CD8_T_CELL_MARKERS

B_CELL_MARKERS = frozenset({"MS4A1"})

CD14_MONOCYTE_MARKERS = frozenset({"CD14", "LYZ"})
FCGR3A_MONOCYTE_MARKERS = frozenset({"FCGR3A", "MS4A7"})
MONOCYTE_MARKERS = CD14_MONOCYTE_MARKERS | FCGR3A_MONOCYTE_MARKERS

NK_CELL_MARKERS = frozenset({"GNLY", "NKG7"})

DENDRITIC_CELL_MARKERS = frozenset({"FCER1A", "CST3"})

PLATELET_MARKERS = frozenset({"PPBP"})

# Union of every panel above - kept for backward compatibility with
# earlier code that tests against "any known PBMC marker" rather
# than a specific lineage.
PBMC_IMMUNE_MARKERS = frozenset().union(
    T_CELL_MARKERS,
    B_CELL_MARKERS,
    MONOCYTE_MARKERS,
    NK_CELL_MARKERS,
    DENDRITIC_CELL_MARKERS,
    PLATELET_MARKERS,
)

LINEAGE_PANELS = {
    "T cell": T_CELL_MARKERS,
    "B cell": B_CELL_MARKERS,
    "Monocyte": MONOCYTE_MARKERS,
    "NK cell": NK_CELL_MARKERS,
    "Dendritic cell": DENDRITIC_CELL_MARKERS,
    "Platelet": PLATELET_MARKERS,
}