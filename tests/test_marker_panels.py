from atlasx.scoring.marker_panels import (
    T_CELL_MARKERS,
    B_CELL_MARKERS,
    MONOCYTE_MARKERS,
    NK_CELL_MARKERS,
    DENDRITIC_CELL_MARKERS,
    PLATELET_MARKERS,
    PBMC_IMMUNE_MARKERS,
    LINEAGE_PANELS,
)


def test_expected_marker_genes_present():

    assert "IL7R" in T_CELL_MARKERS
    assert "CD8A" in T_CELL_MARKERS
    assert "MS4A1" in B_CELL_MARKERS
    assert "CD14" in MONOCYTE_MARKERS
    assert "FCGR3A" in MONOCYTE_MARKERS
    assert "NKG7" in NK_CELL_MARKERS
    assert "FCER1A" in DENDRITIC_CELL_MARKERS
    assert "PPBP" in PLATELET_MARKERS


def test_pbmc_immune_markers_is_union_of_all_panels():

    assert PBMC_IMMUNE_MARKERS == T_CELL_MARKERS | B_CELL_MARKERS | MONOCYTE_MARKERS | NK_CELL_MARKERS | DENDRITIC_CELL_MARKERS | PLATELET_MARKERS


def test_lineage_panels_covers_six_lineages():

    assert set(LINEAGE_PANELS.keys()) == {
        "T cell", "B cell", "Monocyte", "NK cell", "Dendritic cell", "Platelet"
    }