import numpy as np
import logging
from .graph import VariantNotFoundException

class MostSimilarVariantLookup:
    def __init__(self, lookup_array, prob_same_genotype):
        self.lookup_array = lookup_array
        self.prob_same_genotype = prob_same_genotype

    def get_most_similar_variant(self, variant_id):
        return self.lookup_array[variant_id]

    def prob_of_having_the_same_genotype_as_most_similar(self, variant_id):
        return self.prob_same_genotype[variant_id]

    def to_file(self, file_name):
        np.savez(file_name, lookup_array=self.lookup_array, prob_same_genotype=self.prob_same_genotype)

    @classmethod
    def from_file(cls, file_name):
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npy")

        return cls(data["lookup_array"], data["prob_same_genotype"])


class GenotypeMatrixAnalyser:
    def __init__(self, genotype_matrix):
        self.matrix = genotype_matrix.matrix

    def get_most_similar_previous_variant(self, variant_id):
        matrix = self.matrix
        variant_genotypes = matrix[:,variant_id]
        #print("Variant genotypes: %s" % variant_genotypes)
        submatrix_start = 0
        # Only look at 1000 variants back
        if variant_id > 1000:
            submatrix_start = variant_id - 1000

        submatrix = matrix[:,submatrix_start:variant_id]
        #print("Submatrix: %s" % submatrix)
        similarity_scores = np.sum(
                submatrix.transpose() == matrix[:,variant_id],
                axis=1
        )
        #print(similarity_scores)

        most_similar = np.argmax(similarity_scores) + submatrix_start
        value = similarity_scores[most_similar - submatrix_start]
        return most_similar, value

    def analyse(self):
        n_variants = self.matrix.shape[1]
        n_individuals = self.matrix.shape[0]
        most_similar_lookup = np.zeros(n_variants, dtype=np.uint32)
        prob_same_genotype = np.zeros(n_variants, dtype=np.float)

        for variant_id in range(1, self.matrix.shape[1]):
            if variant_id % 1000 == 0:
                logging.info("%d variants analysed" % variant_id)

            most_similar, score = self.get_most_similar_previous_variant(variant_id)
            #logging.info("Most similar to %d is %d with score %d. Genotype distribution: %s" % (variant_id, most_similar, score, np.unique(self.matrix[:,variant_id], return_counts=True)))
            most_similar_lookup[variant_id] = most_similar
            prob_same_genotype[variant_id] = score / n_individuals

        return MostSimilarVariantLookup(most_similar_lookup, prob_same_genotype)


class GenotypeMatrix:
    def __init__(self, matrix):
        self.matrix = matrix

    @classmethod
    def from_nodes_to_haplotypes_and_variants(cls, nodes_to_haplotypes, variants, graph, n_individuals):
        n_variants = len(variants)

        matrix = np.zeros((n_individuals, n_variants), dtype=np.uint8)

        for variant_number, variant in enumerate(variants):
            if variant_number % 100 == 0:
                logging.info("%d variants processeed" % variant_number)

            try:
                reference_node, variant_node = graph.get_variant_nodes(variant)
            except VariantNotFoundException:
                continue

            for individual_id in nodes_to_haplotypes.get_individuals_having_node_pair(reference_node, reference_node):
                matrix[individual_id, variant_number] = 1

            for individual_id in nodes_to_haplotypes.get_individuals_having_node_pair(variant_node, variant_node):
                matrix[individual_id, variant_number] = 2

            for individual_id in nodes_to_haplotypes.get_individuals_having_node_pair(reference_node, variant_node):
                matrix[individual_id, variant_number] = 3

        return cls(matrix)

    def to_file(self, file_name):
        np.save(file_name, self.matrix)

    @classmethod
    def from_file(cls, file_name):
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npy")

        return cls(data)
