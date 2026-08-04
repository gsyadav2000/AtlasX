from epimatch.core.peak import Peak

peak1 = Peak(
    chromosome="chr1",
    start=1000,
    end=1500
)

peak2 = Peak(
    chromosome="chr1",
    start=1400,
    end=1700
)

print("=" * 50)
print("Peak 1")
print("=" * 50)
peak1.summary()

print()

print("=" * 50)
print("Peak 2")
print("=" * 50)
peak2.summary()

print()

print("Overlap :", peak1.overlaps(peak2))

print("Distance to Position :", peak1.distance_to(2000))