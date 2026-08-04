import pickle

from epimatch.scoring.marker_panels import LINEAGE_PANELS

with open("data/processed/dataset_a_reference.pkl", "rb") as f:
    saved = pickle.load(f)

reference_profiles = saved["reference_profiles"]
cluster_sizes = saved["cluster_sizes"]

print("For each cluster, checking whether ANY marker gene from each")
print("lineage panel appears anywhere in the top-100 reference profile -")
print("not just whether it's significant, since a gene can be present")
print("without being frequent enough to register as statistically")
print("significant against a small panel.\n")

for cluster_id in sorted(reference_profiles.keys()):

    profile = reference_profiles[cluster_id]
    size = cluster_sizes.get(cluster_id, "?")

    print(f"--- Cluster {cluster_id} ({size} cells) ---")

    for lineage_name, marker_panel in LINEAGE_PANELS.items():

        present = sorted(marker_panel & profile)

        if present:
            print(f"  {lineage_name:10} FOUND: {', '.join(present)}")
        else:
            print(f"  {lineage_name:10} none of {sorted(marker_panel)} present")

    print()

print("=" * 60)
print("Full reference profile for clusters with no lineage hits so far")
print("=" * 60)
print("(clusters 2, 6, 7, 8 have never shown a significant lineage match)\n")

for cluster_id in [2, 6, 7, 8]:
    if cluster_id not in reference_profiles:
        continue
    print(f"--- Cluster {cluster_id} full profile ---")
    for gene in sorted(reference_profiles[cluster_id]):
        print(f"  {gene}")
    print()