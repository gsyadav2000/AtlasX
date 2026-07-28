from atlasx.core.peak import Peak
from atlasx.core.gene import Gene

peak = Peak.from_string(
    "chr17:7670000-7670500"
)

gene = Gene(
    name="TP53",
    chromosome="chr17",
    start=7661779,
    end=7687550,
    strand="-"
)

distance = peak.distance_to_gene(gene)

print("=" * 50)
print("Peak → Gene Distance")
print("=" * 50)
print(f"Peak : {peak.chromosome}:{peak.start}-{peak.end}")
print(f"Gene : {gene.name}")
print(f"TSS  : {gene.tss}")
print(f"Distance : {distance} bp")
print("=" * 50)