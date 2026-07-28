from atlasx.scoring.regulatory_score import RegulatoryScore

score = RegulatoryScore(
    distance_bp=17300,
    promoter=False,
    enhancer=True
)

print("=" * 50)
print("AtlasX Regulatory Score")
print("=" * 50)

print("Distance Score :", score.distance_score())
print("Promoter Score :", score.promoter_score())
print("Enhancer Score :", score.enhancer_score())

print()

print("Final Score :", score.final_score())

print("=" * 50)