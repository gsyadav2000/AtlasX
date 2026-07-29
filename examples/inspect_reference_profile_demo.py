import pickle

with open("data/processed/dataset_a_reference.pkl", "rb") as f:
    saved = pickle.load(f)

print("Cluster 3 (monocyte) reference profile, top 30 genes:")
for gene in sorted(saved["reference_profiles"][3])[:30]:
    print(f"  {gene}")