from collections import Counter

from atlasx.scoring.marker_enrichment import marker_set_enrichment


def test_marker_enrichment_detects_overrepresentation():

    background = {f"Gene{i}" for i in range(1000)}
    marker_panel = {f"Gene{i}" for i in range(20)}

    # Observed set is heavily loaded with marker genes - far more
    # than the 20/1000 base rate would predict by chance.
    marker_hits = {gene: 1 for gene in list(marker_panel)[:15]}
    non_marker_hits = {f"Gene{i}": 1 for i in range(500, 510)}
    gene_hit_counts = Counter({**marker_hits, **non_marker_hits})

    result = marker_set_enrichment(gene_hit_counts, marker_panel, background)

    assert result["panel_hit_count"] == 15
    assert result["p_value"] < 0.001