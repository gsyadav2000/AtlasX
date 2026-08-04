from epimatch.core.peak import Peak


def test_peak_creation():
    peak = Peak(
        chromosome="chr1",
        start=100,
        end=500,
    )

    assert peak.chromosome == "chr1"
    assert peak.start == 100
    assert peak.end == 500
    assert peak.length == 400