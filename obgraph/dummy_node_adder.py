from alignment_free_graph_genotyper.variants import GenotypeCalls
# Modifies a graph so that deletions/insertions has parallel dummy nodes
import logging
from itertools import product
from .graph import Graph


class DummyNodeAdder:
    def __init__(self, graph,  variants):
        self.graph = graph
        self.variants = variants
        self.mutable_graph = None
        self.current_new_node_id = len(self.graph.nodes) + 1

    def create_new_graph_with_dummy_nodes(self):
        logging.info("Creating a mutable graph that can be changed")
        self.mutable_graph = self.graph.to_mutable_graph()

        logging.info("Adding dummy nodes")
        for i, variant in enumerate(self.variants):
            if i % 100 == 0:
                logging.info("%d variants processed" % i)

            if variant.type == "SNP":
                continue

            self._add_dummy_edges_around_indel(variant)

        logging.info("Creating a new immutable graph from the mutable graph")
        print(self.mutable_graph.get_all_nodes())
        print(self.mutable_graph.node_sequences)
        return Graph.from_mutable_graph(self.mutable_graph)

    def get_nodes_for_inserted_sequence_at_ref_pos(self, inserted_sequence, ref_pos):
        node_before_inserted_nodes = self.graph.get_node_at_ref_offset(ref_pos)
        inserted_nodes = self.mutable_graph.find_nodes_from_node_that_matches_sequence(node_before_inserted_nodes, inserted_sequence)
        if not inserted_nodes:
            logging.error("Could not find inserted nodes for sequence %s at ref pos %d. Node before is %d" % (inserted_sequence, ref_pos, node_before_inserted_nodes))
            raise Exception("Could not find inserted nodes")

        return inserted_nodes

    def _add_dummy_edges_around_indel(self, variant):
        if variant.type == "DELETION":
            inserted_sequence = variant.ref_sequence[1:]
        elif variant.type == "INSERTION":
            inserted_sequence = variant.variant_sequence[1:]
        else:
            raise Exception("Unsupported variant")

        inserted_nodes = self.get_nodes_for_inserted_sequence_at_ref_pos(inserted_sequence, variant.position)

        # Find edges going from any node that goes into the inserted nodes and that ends at any node going out from the end of the inserted nodes
        # These edges are bypassing the inserted nodes
        in_nodes = self.mutable_graph.get_nodes_before(inserted_nodes[0])
        out_nodes = self.mutable_graph.get_edges(inserted_nodes[-1])
        bypass_edges = [(from_node, to_node) for from_node, to_node in product(in_nodes, out_nodes) if to_node in self.mutable_graph.get_edges(from_node)]
        print("BYpass edges: %s" % bypass_edges)
        self._add_new_dummy_node_for_edges(bypass_edges)

    def _add_new_dummy_node_for_edges(self, edges):
        dummy_node = self.current_new_node_id
        self.mutable_graph.add_node(dummy_node)

        # Remove all edges
        for from_node, to_node in edges:
            self.mutable_graph.remove_edge(from_node, to_node)

        # For all unique from_nodes add an edge to the dummy node
        for from_node in set((from_node for from_node, to_node in edges)):
            self.mutable_graph.add_edge(from_node, dummy_node)

        # For all unique to_nodes, add an edge from dummy node to to node
        for to_node in set((to_node for from_node, to_node in edges)):
            print("Adding node from dummy node %d to node %d" % (dummy_node, to_node))
            self.mutable_graph.add_edge(dummy_node, to_node)

        self.current_new_node_id += 1

