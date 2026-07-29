"""
AtlasX Marker Gene Panels

Marker genes are taken from the official Seurat PBMC3k guided
clustering tutorial (Satija Lab), the standard reference marker set
used throughout the single-cell field for identifying PBMC cell
types from clustering results:
https://satijalab.org/seurat/archive/v4.3/pbmc3k_tutorial

CD300E, C5AR2, and CD93 were added to the monocyte panel afterward,
not from the Seurat tutorial - they were observed as prominent,
well-established monocyte/myeloid marker genes in this project's own
dataset-A reference profile for the monocyte cluster, during
cross-dataset validation. The original 4-gene Seurat monocyte panel
(CD14, LYZ, FCGR3A, MS4A7) was too narrow to detect real monocyte
signal that was demonstrably present in that reference profile. This
is a real, defensible expansion (all three are established monocyte
markers in the literature), but it is a project-derived addition, not
part of the original Seurat source - flagged here so that distinction
stays visible rather than getting silently absorbed into "the Seurat
panel."

Two lineages - Dendritic Cells and Platelets/Megakaryocytes - are
included here because they're part of the Seurat reference set, even
though neither has been tested for extensively in this project. Both
have very small panels (1-2 genes), so their statistical power to
detect real signal is correspondingly limited.
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

# CD300E, C5AR2, CD93 added from this project's own data - see
# module docstring for why.
ADDITIONAL_MONOCYTE_MARKERS = frozenset({"CD300E", "C5AR2", "CD93"})

MONOCYTE_MARKERS = (
    CD14_MONOCYTE_MARKERS
    | FCGR3A_MONOCYTE_MARKERS
    | ADDITIONAL_MONOCYTE_MARKERS
)

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