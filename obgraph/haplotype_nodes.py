import numpy as np
from collections import defaultdict
import logging
from .graph import VariantNotFoundException, Graph
import pickle
from multiprocessing import Pool, Process
from graph_kmer_index.shared_mem import to_shared_memory, from_shared_memory
from itertools import repeat
from .traversing import traverse_graph_by_following_nodes
import random


class GenotypeToNodes:
    def __init__(self):
        # Could this just be a map from genotype to two haplotypes?
        pass

    @classmethod
    def make_from_n_random_haplotypes(cls, graph, variants, n_haplotypes=10):
        pass

class HaplotypeToNodes:
    properties = {"_haplotype_to_index", "_haplotype_to_n_nodes", "_nodes"}
    def __init__(self, haplotype_to_index=None, haplotype_to_n_nodes=None, nodes=None):
        self._haplotype_to_index = haplotype_to_index
        self._haplotype_to_n_nodes = haplotype_to_n_nodes
        self._nodes = nodes

    def to_file(self, file_name):
        np.savez(file_name, index=self._haplotype_to_index, n=self._haplotype_to_n_nodes, haplotypes=self._nodes)

    def get_nodes(self, haplotype):
        assert type(haplotype) == int
        index = self._haplotype_to_index[haplotype]
        n = self._haplotype_to_n_nodes[haplotype]

        if n == 0:
            return np.array([])

        return self._nodes[index:index+n]

    def __getitem__(self, item):
        return self.get_nodes(item)

    def get_n_haplotypes_on_nodes_array(self, n_nodes):
        counts = np.zeros(n_nodes)
        for haplotype in range(len(self._haplotype_to_index)):
            nodes = self.get_nodes(haplotype)
            counts[nodes] += 1

        return counts

    @classmethod
    def from_file(cls, file_name):
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npz")

        return cls(data["index"], data["n"], data["haplotypes"])

    def get_new_by_traversing_graph(self, graph, n_haplotypes, store_only_variant_nodes=False):
        haplotype_to_index = []
        haplotype_to_n_nodes = []
        nodes = []

        index = 0
        for haplotype in range(n_haplotypes):
            haplotype_to_index.append(index)
            nodes_in_haplotype = self.get_nodes(haplotype)

            # Traverse graph by following these nodes,
            new_nodes = traverse_graph_by_following_nodes(graph, set(nodes_in_haplotype))
            logging.info("Got %d new nodes by traversing graph for haplotype %d" % (len(new_nodes), haplotype))

            nodes.extend(new_nodes)
            index += len(new_nodes)
            haplotype_to_n_nodes.append(len(new_nodes))

        new = HaplotypeToNodes(np.array(haplotype_to_index, dtype=np.uint32), np.array(haplotype_to_n_nodes, dtype=np.uint32), np.array(nodes, dtype=np.uint32))
        print("N nodes: %s" % new._haplotype_to_n_nodes)

        return new

    @classmethod
    def from_flat_haplotypes_and_nodes(cls, haplotypes, nodes):
        assert len(haplotypes) == len(nodes)

        logging.info("Creating numpy arrays from %d nodes" % len(nodes))
        nodes = np.array(nodes, dtype=np.uint32)
        haplotypes = np.array(haplotypes, np.uint16)

        logging.info("Sorting haplotypes and nodes")
        sorting = np.argsort(haplotypes)
        nodes = nodes[sorting]
        haplotypes = haplotypes[sorting]

        # Find positions where nodes change (these are our index entries)
        logging.info("Making index")
        diffs = np.ediff1d(haplotypes, to_begin=1)
        unique_entry_positions = np.nonzero(diffs)[0]
        print(unique_entry_positions)
        unique_haplotypes = haplotypes[unique_entry_positions]

        lookup_size = int(np.max(haplotypes)) + 1
        lookup = np.zeros(lookup_size, dtype=np.uint32)
        lookup[unique_haplotypes] = unique_entry_positions
        n_entries = np.ediff1d(unique_entry_positions, to_end=len(haplotypes) - unique_entry_positions[-1])
        print("N entries: %s" % n_entries)
        n_nodes = np.zeros(lookup_size, dtype=np.uint32)
        n_nodes[unique_haplotypes] = n_entries

        return cls(lookup, n_nodes, nodes)

    @staticmethod
    def _multiprocess_wrapper(shared_memory_graph_name, variants, limit_to_n_haplotypes=10):
        graph = from_shared_memory(Graph, shared_memory_graph_name)
        return HaplotypeToNodes.get_flat_haplotypes_and_nodes_from_graph_and_variants(graph, variants, limit_to_n_haplotypes)

    @classmethod
    def make_from_n_random_haplotypes(cls, graph, variants, n_haplotypes=10, weight_by_allele_frequency=True):
        # Simple way of making "arbitrary" haplotypes, just give every nth variant to every haplotype
        if not weight_by_allele_frequency:
            logging.info("Will not weight by allele frequency, will divide haplotypes equally between ref and var nodes")
        current_haplotype = 0

        haplotype_ids = list(range(n_haplotypes))

        flat_haplotypes = []
        flat_nodes = []
        for i, variant in enumerate(variants):
            if i % 100000 == 0:
                logging.info("%d variants processed" % i)

            try:
                reference_node, variant_node = graph.get_variant_nodes(variant)
            except VariantNotFoundException:
                continue

            # Select number of haplotypes on variant node by allele frequency, always minimum 1
            # never more than n_haplotypes-1 (guaranteeing min 1 on ref and min 1 on alt)
            if weight_by_allele_frequency:
                n_haplotypes_on_variant_node = int(round(min(max(1, graph.get_node_allele_frequency(variant_node) * n_haplotypes), n_haplotypes-1)))
            else:
                assert n_haplotypes % 2 == 0, "Number of haplotypes most be divisible by 2"
                n_haplotypes_on_variant_node = n_haplotypes // 2

            haplotypes_on_variant_node = set(random.sample(haplotype_ids, n_haplotypes_on_variant_node))

            for haplotype in haplotype_ids:
                flat_haplotypes.append(haplotype)
                if haplotype in haplotypes_on_variant_node:
                    flat_nodes.append(variant_node)
                else:
                    flat_nodes.append(reference_node)

            """
            # Give variant node to current haplotype
            flat_haplotypes.append(current_haplotype % n_haplotypes)
            flat_nodes.append(variant_node)

            # Give ref node to all others
            for haplotype in range(n_haplotypes):
                if haplotype != current_haplotype  % n_haplotypes:
                    flat_haplotypes.append(haplotype)
                    flat_nodes.append(reference_node)

            current_haplotype += 1
            """

        return cls.from_flat_haplotypes_and_nodes(flat_haplotypes, flat_nodes)

    @staticmethod
    def get_flat_haplotypes_and_nodes_from_graph_and_variants(graph, variants, limit_to_n_haplotypes=10):
        logging.info("Processing %d variants" % len(variants))
        flat_haplotypes = []
        flat_nodes = []
        haplotypes = list(range(0, limit_to_n_haplotypes))
        for i, variant in enumerate(variants):
            if i % 1000000 == 0:
                logging.info("%d variants processed" % i)

            try:
                reference_node, variant_node = graph.get_variant_nodes(variant)
            except VariantNotFoundException:
                continue

            genotypes = variant.vcf_line.split()[9:]
            for haplotype in haplotypes:
                individual_number = haplotype // 2
                haplotype_number = haplotype - individual_number * 2
                haplotype_string = genotypes[individual_number].replace("/", "|").split("|")[haplotype_number]
                if haplotype_string == "1":
                    # Follows the variant, add variant node here. Do not store reference node in order to svae space
                    flat_haplotypes.append(haplotype)
                    flat_nodes.append(variant_node)
                #else:
                #flat_nodes.append(reference_node)

        return flat_haplotypes, flat_nodes

    @classmethod
    def from_graph_and_variants(cls, graph, variants, limit_to_n_haplotypes=10, n_threads=10):
        # Flat structures used to make the index later
        flat_nodes = []
        flat_haplotypes = []

        pool = Pool(n_threads)
        shared_memory_graph_name = "graph_shared"
        to_shared_memory(graph, shared_memory_graph_name)

        for haplotypes, nodes in pool.starmap(HaplotypeToNodes._multiprocess_wrapper, zip(repeat(shared_memory_graph_name), variants.get_chunks(chunk_size=1000), repeat(limit_to_n_haplotypes))):
            logging.info("Done with 1 iteration")
            flat_haplotypes.extend(haplotypes)
            flat_nodes.extend(nodes)
            logging.info("Added nodes and haplotypes")

        logging.info("Done processing all variants")

        return cls.from_flat_haplotypes_and_nodes(flat_haplotypes, flat_nodes)


# Simple placeholder class for representing a matrix
# rows are haplotypes
# columns contain all nodes covered by that haplotype
class NodeToHaplotypes:
    def __init__(self, nodes_to_index, nodes_to_n_haplotypes, haplotypes):
        self._nodes_to_index = nodes_to_index
        self._nodes_to_n_haplotypes = nodes_to_n_haplotypes
        self._haplotypes = haplotypes
        logging.info("Finding n haplotypes (could be slow?)")
        self._n_haplotypes = np.max(haplotypes)
        logging.info("Done finding n haplotypes")

    def get_haplotypes_on_node(self, node):
        if node >= len(self._nodes_to_index):
            return np.array([])

        index_pos = self._nodes_to_index[node]
        n = self._nodes_to_n_haplotypes[node]
        if n == 0:
            return np.array([])

        return self._haplotypes[index_pos:index_pos+n]

    def get_individuals_having_node_pair(self, node1, node2):
        node1_haplotypes = set(self.get_haplotypes_on_node(node1))
        node2_haplotypes = set(self.get_haplotypes_on_node(node2))

        individuals = set()
        for individual_id in range(self._n_haplotypes // 2):
            haplotype1 = individual_id * 2
            haplotype2 = haplotype1 + 1

            if (haplotype1 in node1_haplotypes and haplotype2 in node2_haplotypes) or (haplotype2 in node1_haplotypes and haplotype1 in node2_haplotypes):
                individuals.add(individual_id)

        return individuals


    @classmethod
    def from_haplotype_nodes(cls, haplotype_nodes):

        # "flat" lists of nodes and corresponding haplotypes having those nodes
        nodes = []
        haplotypes = []

        n_haplotypes = haplotype_nodes.nodes.shape[0]
        for haplotype in range(n_haplotypes):
            logging.info("Processing haplotype %d" % haplotype)
            for node in haplotype_nodes.nodes[haplotype, :]:
                if node > 0:
                    nodes.append(node)
                    haplotypes.append(haplotype)

        nodes = np.array(nodes, dtype=np.uint32)
        haplotypes = np.array(haplotypes, np.uint16)

        sorting = np.argsort(nodes)
        nodes = nodes[sorting]
        haplotypes = haplotypes[sorting]

        # Find positions where nodes change (these are our index entries)
        diffs = np.ediff1d(nodes, to_begin=1)
        unique_entry_positions = np.nonzero(diffs)[0]
        unique_nodes = nodes[unique_entry_positions]

        lookup_size = int(np.max(nodes)) + 1
        lookup = np.zeros(lookup_size, dtype=np.int)
        lookup[unique_nodes] = unique_entry_positions
        n_entries = np.ediff1d(unique_entry_positions, to_end=len(nodes) - unique_entry_positions[-1])
        n_haplotypes = np.zeros(lookup_size, dtype=np.uint16)
        n_haplotypes[unique_nodes] = n_entries

        return cls(lookup, n_haplotypes, haplotypes)

    def to_file(self, file_name):
        np.savez(file_name, index=self._nodes_to_index, n=self._nodes_to_n_haplotypes, haplotypes=self._haplotypes)

    @classmethod
    def from_file(cls, file_name):
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npz")

        return cls(data["index"], data["n"], data["haplotypes"])




class HaplotypeNodes:
    def __init__(self, nodes, n_haplotypes_on_node):
        self.nodes = nodes
        self.n_haplotypes_on_node = n_haplotypes_on_node


    def to_file(self, file_name):
        np.savez(file_name, nodes=self.nodes, n_haplotypes_on_node=self.n_haplotypes_on_node)

    def __get__(self, item):
        return self.nodes[item]

    @classmethod
    def from_file(cls, file_name):
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npz")

        return cls(data["nodes"], data["n_haplotypes_on_node"])

    @classmethod
    def from_graph_and_variants(cls, graph, variants, limit_to_n_haplotypes=10):

        # First find all variant nodes that the haplotype has
        haplotypes = list(range(0, limit_to_n_haplotypes))
        variant_nodes_in_haplotype = defaultdict(set)
        for i, variant in enumerate(variants):
            if i % 1000 == 0:
                logging.info("%d variants processed" % i)

            try:
                reference_node, variant_node = graph.get_variant_nodes(variant)
            except VariantNotFoundException:
                continue

            if variant.position == 4871514:
                logging.info("Variant 4871514 has nodes %d/%d" % (reference_node, variant_node))

            genotypes = variant.vcf_line.split()[9:]
            for haplotype in haplotypes:
                individual_number = haplotype // 2
                haplotype_number = haplotype - individual_number * 2
                haplotype_string = genotypes[individual_number].replace("/", "|").split("|")[haplotype_number]
                if haplotype_string == "1":
                    # Follows the variant, add variant node here
                    variant_nodes_in_haplotype[haplotype].add(variant_node)
                else:
                    variant_nodes_in_haplotype[haplotype].add(reference_node)

        # Iterate graph
        logging.info("Iterating graph for each haplotype")
        nodes = np.zeros((len(haplotypes), len(graph.nodes)), dtype=np.uint32)
        n_haplotypes_on_node = np.zeros(len(graph.nodes) + 1, dtype=np.uint32)

        for haplotype in haplotypes:
            logging.info("Handling haplotype %d" % haplotype)
            current_node = graph.get_first_node()
            i = 0
            while True:
                nodes[haplotype, i] = current_node
                n_haplotypes_on_node[current_node] += 1

                next_nodes = graph.get_edges(current_node)
                if len(next_nodes) == 0:
                    break

                next_node = None
                if len(next_nodes) == 1:
                    next_node = next_nodes[0]
                else:
                    for potential_next in next_nodes:
                        if potential_next in variant_nodes_in_haplotype[haplotype]:
                            next_node = potential_next

                if next_node is None:
                    logging.error("Could not find next node from node %d" % current_node)
                    logging.error("Possible next nodes are %s" % next_nodes)
                    raise Exception("")

                current_node = next_node
                i += 1

            nodes[haplotype, i] = current_node
            n_haplotypes_on_node[current_node] += 1

        return cls(nodes, n_haplotypes_on_node)
