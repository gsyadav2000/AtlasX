"""
EpiMatch Gene Database
"""

from epimatch.core.gene import Gene


class GeneDatabase:

    def __init__(self, gtf_file, protein_coding_only=False):

        self.gtf_file = gtf_file
        self.genes = []
        self.protein_coding_only = protein_coding_only

    def load(self):

        print("Loading GENCODE annotation...")

        with open(self.gtf_file, "r", encoding="utf-8") as f:

            for line in f:

                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")

                if len(fields) != 9:
                    continue

                chromosome = fields[0]
                feature = fields[2]

                if feature != "gene":
                    continue

                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]

                attributes = fields[8]

                gene_name = None
                gene_type = None

                for item in attributes.split(";"):

                    item = item.strip()

                    if item.startswith("gene_name"):
                        gene_name = item.split('"')[1]

                    elif item.startswith("gene_type"):
                        gene_type = item.split('"')[1]

                if gene_name is None:
                    continue

                if self.protein_coding_only and gene_type != "protein_coding":
                    continue

                gene = Gene(
                    name=gene_name,
                    chromosome=chromosome,
                    start=start,
                    end=end,
                    strand=strand
                )

                self.genes.append(gene)

        print(f"Loaded {len(self.genes):,} genes.")

        return self.genes