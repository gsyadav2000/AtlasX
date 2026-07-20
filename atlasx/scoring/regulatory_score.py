"""
AtlasX Regulatory Confidence Score
Version 0.1
"""


class RegulatoryScore:

    def __init__(
        self,
        distance_bp,
        promoter=False,
        enhancer=False
    ):

        self.distance_bp = distance_bp
        self.promoter = promoter
        self.enhancer = enhancer

    def distance_score(self):
        """
        Convert genomic distance into a score.
        Smaller distance = higher score.
        """

        if self.distance_bp < 1000:
            return 1.0

        elif self.distance_bp < 10000:
            return 0.8

        elif self.distance_bp < 50000:
            return 0.5

        elif self.distance_bp < 100000:
            return 0.2

        return 0.05

    def promoter_score(self):

        return 1.0 if self.promoter else 0.0

    def enhancer_score(self):

        return 1.0 if self.enhancer else 0.0

    def final_score(self):

        score = (
            0.5 * self.distance_score()
            + 0.3 * self.promoter_score()
            + 0.2 * self.enhancer_score()
        )

        return round(score, 3)