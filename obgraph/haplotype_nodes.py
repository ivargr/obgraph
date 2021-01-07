import numpy as np
from collections import defaultdict
import logging
from .graph import VariantNotFoundException

# Simple placeholder class for representing a matrix
# rows are haplotypes
# columns contain all nodes covered by that haplotype

class HaplotypeNodes:
    def __init__(self, nodes):
        self.nodes = nodes

    def to_file(self, file_name):
        np.save(file_name, self.nodes)

    @classmethod
    def from_file(cls, file_name):
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npz")

        return cls(data)

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
                    if variant.position == 4871514:
                        logging.info("Haplotype %d has variant" % (haplotype))
                    # Follows the variant, add variant node here
                    variant_nodes_in_haplotype[haplotype].add(variant_node)
                else:
                    if variant.position == 4871514:
                        logging.info("Haplotype %d does not have variant" % (haplotype))
                    variant_nodes_in_haplotype[haplotype].add(reference_node)

        # Iterate graph
        logging.info("Iterating graph for each haplotype")
        nodes = np.zeros((len(haplotypes), len(graph.nodes)), dtype=np.uint32)
        for haplotype in haplotypes:
            logging.info("Handling haplotype %d" % haplotype)
            current_node = graph.get_first_node()
            i = 0
            while True:
                nodes[haplotype, i] = current_node

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

        return cls(nodes)
