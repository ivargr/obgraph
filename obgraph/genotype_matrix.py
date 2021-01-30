import numpy as np
import logging
import time
from .graph import VariantNotFoundException
from graph_kmer_index.shared_mem import to_shared_memory, from_shared_memory
from multiprocessing import Pool, Process
from .variant_to_nodes import VariantToNodes

class MostSimilarVariantLookup:
    properties = {"lookup_array", "prob_same_genotype"}
    def __init__(self, lookup_array=None, prob_same_genotype=None):
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


class GenotypeFrequencies:
    properties = {"homo_ref", "homo_alt", "hetero"}
    def __init__(self, homo_ref=None, homo_alt=None, hetero=None):
        self.homo_ref = homo_ref
        self.homo_alt = homo_alt
        self.hetero = hetero

    @staticmethod
    def create_using_shared_memory(variant_interval):
        from_variant, to_variant = variant_interval
        logging.info("Creating on interval %d-%d" % (from_variant, to_variant))

        genotype_matrix = from_shared_memory(GenotypeMatrix, "genotype_matrix_shared_for_frequencies")
        genotype_frequencies = from_shared_memory(GenotypeFrequencies, "genotype_frequencies_shared")

        n_variants = genotype_matrix.matrix.shape[1]
        n_individuals = len(np.where(genotype_matrix.matrix[:,0])[0] != 0)  # can be zeros for non-individuals, so all non-zero is an individual
        # Less memory hungry, but slower

        for numeric_genotype, array in zip([1, 2, 3], [genotype_frequencies.homo_ref, genotype_frequencies.homo_alt, genotype_frequencies.hetero]):
            logging.info("Finding for genotype %d" % numeric_genotype)
            prev_time = time.time()
            for variant_id in range(from_variant, to_variant):
                if variant_id % 100000 == 0:
                    logging.info("%d/%d variants processed (genotype now is %d). Prev 100k processed in %.3f s" % (variant_id-from_variant, to_variant-from_variant, numeric_genotype, time.time()-prev_time))
                    prev_time = time.time()

                array[variant_id] = len(np.where(genotype_matrix.matrix[:, variant_id] == numeric_genotype)[0]) / n_individuals

    @classmethod
    def from_genotype_matrix(cls, genotype_matrix, n_threads=10):
        to_shared_memory(genotype_matrix, "genotype_matrix_shared_for_frequencies")

        n_variants = genotype_matrix.matrix.shape[1]
        n_individuals = len(np.where(genotype_matrix.matrix[:,0])[0] != 0)  # can be zeros for non-individuals, so all non-zero is an individual
        logging.info("Assumes there are %d individuals and %d variants" % (n_individuals, n_variants))
        data = {1: np.zeros(n_variants, dtype=float), 2: np.zeros(n_variants, dtype=float), 3: np.zeros(n_variants, dtype=float)}
        genotype_frequences = cls(data[1], data[2], data[3])
        to_shared_memory(genotype_frequences, "genotype_frequencies_shared")

        intervals = [int(i) for i in np.linspace(0, n_variants, n_threads)]
        variant_intervals = [(from_id, to_id) for from_id, to_id in zip(intervals[0:-1], intervals[1:])]
        logging.info("Will analyse intervals: %s" % variant_intervals)

        pool = Pool(n_threads)

        for result in pool.imap(GenotypeFrequencies.create_using_shared_memory, variant_intervals):
            logging.info("Done with one job")

        """
        for numeric_genotype, array in data.items():
            logging.info("Finding for genotype %d" % numeric_genotype)
            # the second index from np where gives the columns that have a hit, every column 1 time for each hit
            column_hits = np.where(genotype_matrix.matrix == numeric_genotype)[1]
            logging.info("Making frequencies")
            unique_columns, n_hits_per_column = np.unique(column_hits, return_counts=True)
            data[numeric_genotype][unique_columns] = n_hits_per_column / n_individuals
        """
        """
        # Less memory hungry, but slower
        for numeric_genotype, array in data.items():
            logging.info("Finding for genotype %d" % numeric_genotype)
            for variant_id in range(n_variants):
                if variant_id % 10000 == 0:
                    logging.info("%d variants processed" % variant_id)

                array[variant_id] = len(np.where(genotype_matrix.matrix[:,variant_id] == numeric_genotype)[0]) / n_individuals
        """
        return from_shared_memory(GenotypeFrequencies, "genotype_frequencies_shared")
        #return cls(data[1], data[2], data[3])

    def get_frequencies_for_variant(self, variant_id):
        return self.homo_ref[variant_id], self.homo_alt[variant_id], self.hetero[variant_id]

    @classmethod
    def from_file(cls, file_name):
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npz")

        return cls(data["homo_ref"], data["homo_alt"], data["hetero"])

    def to_file(self, file_name):
        np.savez(file_name, homo_ref=self.homo_ref, homo_alt=self.homo_alt, hetero=self.hetero)


class GenotypeMatrixAnalyser:
    def __init__(self, genotype_matrix):
        self.matrix = genotype_matrix

    @staticmethod
    def analyse_variants_on_shared_memody(variant_interval):
        from_id, to_id = variant_interval
        if from_id == 0:
            from_id = 1
        logging.info("Analysing variant %d to %d in one job" % (from_id, to_id))
        matrix = from_shared_memory(GenotypeMatrix, "genotype_matrix")
        lookup = from_shared_memory(MostSimilarVariantLookup, "most_similar_variant_lookup")
        n_individuals = matrix.matrix.shape[0]
        prev_time = time.time()
        for i, variant_id in enumerate(range(from_id, to_id)):
            if i % 5000 == 0 and i > 0:
                logging.info("%d/%d variants analysed (last 5k analyse in %.3f s)" % (i, to_id-from_id, time.time()-prev_time))
                prev_time = time.time()

            most_similar, score = matrix.get_most_similar_previous_variant(variant_id)
            #logging.info("Most similar to %d is %d with score %d. Genotype distribution: %s" % (variant_id, most_similar, score, np.unique(self.matrix[:,variant_id], return_counts=True)))
            lookup.lookup_array[variant_id] = most_similar
            lookup.prob_same_genotype[variant_id] = score / n_individuals


    def analyse(self, n_threads=10):
        n_variants = self.matrix.matrix.shape[1]
        n_individuals = self.matrix.matrix.shape[0]

        most_similar_lookup = np.zeros(n_variants, dtype=np.uint32)
        prob_same_genotype = np.zeros(n_variants, dtype=np.float)

        lookup = MostSimilarVariantLookup(most_similar_lookup, prob_same_genotype)
        to_shared_memory(self.matrix, "genotype_matrix")
        to_shared_memory(lookup, "most_similar_variant_lookup")

        intervals = [int(i) for i in np.linspace(0, n_variants, n_threads)]
        variant_intervals = [(from_id, to_id) for from_id, to_id in zip(intervals[0:-1], intervals[1:])]
        logging.info("Will analyse intervals: %s" % variant_intervals)

        pool = Pool(n_threads)

        for result in pool.imap(GenotypeMatrixAnalyser.analyse_variants_on_shared_memody, variant_intervals):
            logging.info("Done with one job")

        lookup = from_shared_memory(MostSimilarVariantLookup, "most_similar_variant_lookup")

        return lookup


class GenotypeMatrix:
    properties = {"matrix"}

    def __init__(self, matrix=None):
        self.matrix = matrix

    @classmethod
    def from_variants(cls, variants, n_individuals, n_variants, n_threads=10, chunk_size=10000):
        matrix = np.zeros((n_individuals, n_variants), dtype=np.uint8)
        matrix = cls(matrix)
        logging.info("Putting genotype matrix in shared memory")
        to_shared_memory(matrix, "genotype_matrix")

        logging.info("Getting variant chunks")
        variant_chunks = variants.get_chunks(chunk_size=chunk_size)

        pool = Pool(n_threads)

        i = 0
        for result in pool.imap(GenotypeMatrix.fill_shared_memory_matrix_with_variants, variant_chunks):
            i += 1
            logging.info("Done with %d variant chunks" % i)

        logging.info("Done with all variant chunks")
        matrix = from_shared_memory(GenotypeMatrix, "genotype_matrix")
        return cls(matrix.matrix)

    def get_most_similar_previous_variant(self, variant_id):
        matrix = self.matrix
        #variant_genotypes = matrix[:,variant_id]
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

    @staticmethod
    def fill_shared_memory_matrix_with_variants(variants):
        matrix = from_shared_memory(GenotypeMatrix, "genotype_matrix")

        for variant in variants:
            variant_number = variant.vcf_line_number
            if variant_number % 10000 == 0:
                logging.info("%d variants processeed" % variant_number)

            for individual_id, genotype in variant.get_individuals_and_numeric_genotypes():
                matrix.matrix[individual_id, variant_number] = genotype

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
