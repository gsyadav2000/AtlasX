"""
AtlasX Marker Gene Panels

IMPORTANT LIMITATION: this is a small, hand-curated list of
well-established, textbook PBMC/immune marker genes, not a proper
downloaded gene set from a curated database (e.g. MSigDB, GO, or
Enrichr's Human_Gene_Atlas). It is useful as a quick sanity check,
but every gene here was selected from general knowledge of immune
cell biology rather than from a verifiable, versioned source file.
For any result meant to support a real claim or publication, this
panel should be replaced with a gene set downloaded directly from
a proper source (e.g. a GO Biological Process term such as
"immune system process" via AmiGO/QuickGO, or an MSigDB immune
signature), so the exact gene membership is traceable and citable.
"""

PBMC_IMMUNE_MARKERS = frozenset({
    # T cells
    "CD3D", "CD3E", "CD3G", "CD8A", "CD4", "IL7R", "CCR7",
    "TCF7", "LEF1", "SELL", "FOXP3", "IL2RA", "BCL11B",
    "GZMB", "GZMK", "PRF1", "NR4A1",
    # B cells
    "CD19", "MS4A1", "CD79A", "CD79B", "VPREB3", "IGHM",
    # Monocytes / myeloid
    "CD14", "LYZ", "FCN1", "S100A8", "S100A9", "ITGAM",
    "CSF1R", "FCGR3A", "CD80", "CD86", "TLR4", "PILRA",
    # NK cells
    "NKG7", "GNLY", "KLRD1", "CX3CR1",
    # General immune / signaling
    "PTPRC", "PTPN6", "TNF", "IFNG", "CCL5", "HLA-DRA", "HLA-DRB1",
    # Granulocyte
    "MPO", "ELANE",
})