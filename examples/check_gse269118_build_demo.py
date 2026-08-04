from epimatch.io.bed_reader import BEDReader

peaks = BEDReader("data/raw/ingested/GSE269118_extracted/GSM8306617_Patient1_techrep1_summits.bed.gz").load()

print("First 5 peaks from GSE269118:")
for peak in peaks[:5]:
    print(f"  {peak.chromosome}:{peak.start}-{peak.end}")