"""
EpiMatch Marker Gene Panels

Marker genes are taken from the official Seurat PBMC3k guided
clustering tutorial (Satija Lab):
https://satijalab.org/seurat/archive/v4.3/pbmc3k_tutorial

Several genes below were added afterward, not from the Seurat
tutorial - observed as prominent, specific marker genes in this
project's own dataset-A reference profiles during validation, then
checked against known literature before inclusion:

- CD300E, C5AR2, CD93: monocyte, found in cluster 3's reference
  profile (see project history) - established monocyte/myeloid genes.
- POU2AF1: B cell, found in cluster 6's reference profile. NOTE: the
  gene BCR ALSO appeared in this profile and looked like an obvious
  B-cell-receptor marker by name, but the gene symbol "BCR" refers to
  Breakpoint Cluster Region (the BCR-ABL1 fusion gene in CML), a
  broadly-expressed gene unrelated to the B-cell receptor complex -
  deliberately NOT included here despite the misleading name.
- SPIB: plasmacytoid dendritic cell (pDC), found in cluster 7's
  reference profile. This is a different DC subtype than the
  original Seurat conventional-DC panel (FCER1A, CST3) detects, which
  is why that panel found nothing for a real pDC-like cluster.

As with the monocyte panel expansion, these additions were chosen
after seeing which genes were present in the data being tested - real
and literature-supported, but not a blind test. See project history
for the fuller caveat.

Cluster 8 in dataset A shows a strong immediate-early-gene signature
(FOS, FOSB, JUNB, EGR1, NR4A1, NR4A2) that looks like a dissociation-
stress technical artifact rather than a distinct cell type (a
documented phenomenon in single-cell literature), so no panel was
created to "confirm" it - manufacturing a marker set for a likely
artifact would be circular, not real validation.
"""

CD4_T_CELL_MARKERS = frozenset({"IL7R", "CCR7"})
CD8_T_CELL_MARKERS = frozenset({"CD8A"})

T_CELL_MARKERS = CD4_T_CELL_MARKERS | CD8_T_CELL_MARKERS

# POU2AF1 added - see module docstring. BCR deliberately excluded
# despite the misleading gene symbol name.
B_CELL_MARKERS = frozenset({"MS4A1", "POU2AF1"})

CD14_MONOCYTE_MARKERS = frozenset({"CD14", "LYZ"})
FCGR3A_MONOCYTE_MARKERS = frozenset({"FCGR3A", "MS4A7"})
ADDITIONAL_MONOCYTE_MARKERS = frozenset({"CD300E", "C5AR2", "CD93"})

MONOCYTE_MARKERS = (
    CD14_MONOCYTE_MARKERS
    | FCGR3A_MONOCYTE_MARKERS
    | ADDITIONAL_MONOCYTE_MARKERS
)

NK_CELL_MARKERS = frozenset({"GNLY", "NKG7"})

# Conventional dendritic cells (original Seurat panel).
DENDRITIC_CELL_MARKERS = frozenset({"FCER1A", "CST3"})

# Plasmacytoid dendritic cells - biologically distinct subtype from
# conventional DCs, added separately rather than merged into
# DENDRITIC_CELL_MARKERS, since they use a different marker set.
PLASMACYTOID_DC_MARKERS = frozenset({"SPIB"})

PLATELET_MARKERS = frozenset({"PPBP"})

PBMC_IMMUNE_MARKERS = frozenset().union(
    T_CELL_MARKERS,
    B_CELL_MARKERS,
    MONOCYTE_MARKERS,
    NK_CELL_MARKERS,
    DENDRITIC_CELL_MARKERS,
    PLASMACYTOID_DC_MARKERS,
    PLATELET_MARKERS,
)

LINEAGE_PANELS = {
    "T cell": T_CELL_MARKERS,
    "B cell": B_CELL_MARKERS,
    "Monocyte": MONOCYTE_MARKERS,
    "NK cell": NK_CELL_MARKERS,
    "Dendritic cell (conventional)": DENDRITIC_CELL_MARKERS,
    "Dendritic cell (plasmacytoid)": PLASMACYTOID_DC_MARKERS,
    "Platelet": PLATELET_MARKERS,
}