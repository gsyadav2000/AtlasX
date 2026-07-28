"""
AtlasX Marker Gene Panels

IMPORTANT LIMITATION: these are small, hand-curated lists of
well-established, textbook PBMC/immune marker genes, not proper
downloaded gene sets from a curated database (e.g. MSigDB, GO, or
Enrichr's Human_Gene_Atlas). They are useful as a quick sanity
check, but every gene here was selected from general knowledge of
immune cell biology rather than from a verifiable, versioned source
file. For any result meant to support a real claim or publication,
these panels should be replaced with gene sets downloaded directly
from a proper source (e.g. a GO Biological Process term such as
"immune system process" via AmiGO/QuickGO, or an MSigDB immune
signature), so exact gene membership is traceable and citable.
"""

T_CELL_MARKERS = frozenset({
    "CD3D", "CD3E", "CD3G", "CD8A", "CD4", "IL7R", "CCR7",
    "TCF7", "LEF1", "SELL", "FOXP3", "IL2RA", "BCL11B",
    "GZMB", "GZMK", "PRF1", "NR4A1",
})

B_CELL_MARKERS = frozenset({
    "CD19", "MS4A1", "CD79A", "CD79B", "VPREB3", "IGHM",
})

MONOCYTE_MARKERS = frozenset({
    "CD14", "LYZ", "FCN1", "S100A8", "S100A9", "ITGAM",
    "CSF1R", "FCGR3A", "CD80", "CD86", "TLR4", "PILRA",
})

NK_CELL_MARKERS = frozenset({
    "NKG7", "GNLY", "KLRD1", "CX3CR1",
})

GENERAL_IMMUNE_MARKERS = frozenset({
    "PTPRC", "PTPN6", "TNF", "IFNG", "CCL5", "HLA-DRA", "HLA-DRB1",
    "MPO", "ELANE",
})

# Union of every panel above - kept for backward compatibility with
# earlier code that tests against "any immune marker" rather than a
# specific lineage.
PBMC_IMMUNE_MARKERS = frozenset().union(
    T_CELL_MARKERS,
    B_CELL_MARKERS,
    MONOCYTE_MARKERS,
    NK_CELL_MARKERS,
    GENERAL_IMMUNE_MARKERS,
)

LINEAGE_PANELS = {
    "T cell": T_CELL_MARKERS,
    "B cell": B_CELL_MARKERS,
    "Monocyte": MONOCYTE_MARKERS,
    "NK cell": NK_CELL_MARKERS,
}