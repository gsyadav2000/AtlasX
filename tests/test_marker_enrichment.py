from atlasx.scoring.marker_enrichment import marker_set_enrichment, genes_above_frequency


def test_marker_enrichment_detects_overrepresentation():

    background = {f"Gene{i}" for i in range(1000)}
    marker_panel = {f"Gene{i}" for i in range(20)}

    # Observed set is heavily loaded with marker genes - far more
    # than the 20/1000 base rate would predict by chance.
    observed_genes = set(list(marker_panel)[:15]) | {f"Gene{i}" for i in range(500, 510)}

    result = marker_set_enrichment(observed_genes, marker_panel, background)

    assert result["panel_hit_count"] == 15
    assert result["p_value"] < 0.001


def test_genes_above_frequency_filters_by_fraction():

    gene_hit_counts = {
        "AlwaysPresent": 10,
        "HalfPresent": 5,
        "RarelyPresent": 1,
    }

    frequent = genes_above_frequency(gene_hit_counts, num_cells=10, min_fraction=0.3)

    assert frequent == {"AlwaysPresent", "HalfPresent"}
    assert "RarelyPresent" not in frequent


def test_genes_above_frequency_handles_zero_cells():

    result = genes_above_frequency({}, num_cells=0, min_fraction=0.1)

    assert result == set()