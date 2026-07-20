from atlasx.core.peak import Peak

peak1 = Peak(
    chromosome="chr1",
    start=1000,
    end=1500,
    accessibility=7.5
)

peak2 = Peak(
    chromosome="chr1",
    start=1400,
    end=1700,
    accessibility=6.2
)

peak1.summary()

print("Overlap :", peak1.overlaps(peak2))

print("Distance to Gene :", peak1.distance_to(2000))