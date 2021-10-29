from .variants import VcfVariants
import logging
from collections import defaultdict
from .mutable_graph import MutableGraph
from .graph import Graph
from .dummy_node_adder import DummyNodeAdder
import numpy as np

class GraphConstructor:
    def __init__(self, reference_sequence, variants: VcfVariants):
        self.reference_sequence = reference_sequence
        self.variants = variants
        self.breakpoints = []
        self.make_sorted_breakpoints()

        self._ref_pos_to_node_after = defaultdict(list)
        self._node_to_ref_pos_after = defaultdict(list)
        self._ref_pos_to_node_start = defaultdict(list)

        self._node_ids = []
        self._node_sizes = []
        self._node_sequences = []
        self._reference_nodes = []
        self._current_node_id = 1
        self._edges_from = []
        self._edges_to = []

        self._edges_from_ref_pos_to_node = []
        self._deletions = defaultdict(list)  # A lookup from a ref pos to another ref pos representing a deltion
        self._edges_added = set()

        self._mutable_graph = MutableGraph()


        self.make_nodes()
        self.make_edges()
        #self._mutable_graph.linear_ref_nodes = self._reference_nodes
        self._graph = None
        self._graph_with_dummy_nodes = None
        logging.info("Making a graph")
        self._graph = self.get_graph()
        self._graph.to_file("tmpgraph")
        logging.info("Wrote a tmpgraph")
        logging.info("Adding dummy nodes for indels")
        self.add_dummy_nodes()



    def get_graph(self):
        return Graph.from_mutable_graph(self._mutable_graph)
        """
        return Graph.from_flat_nodes_and_edges(np.array(self._node_ids), np.array(self._node_sequences, dtype=object),
                                               np.array(self._node_sizes), np.array(self._edges_from),
                                               np.array(self._edges_to), np.array(self._reference_nodes), np.array([1]))
        """

    def get_graph_with_dummy_nodes(self):
        return self._graph_with_dummy_nodes


    def add_dummy_nodes(self):
        dummy_node_adder = DummyNodeAdder(self._graph, self.variants)
        new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes(self._mutable_graph)
        self._graph_with_dummy_nodes = new_graph

    def _make_edge(self, from_node, to_node, ref_pos_before_to_node):
        if (from_node, to_node) in self._edges_added:
            return

        assert to_node != from_node, "Trying to make edge from same node to same node %d to %d" % (from_node, to_node)

        self._edges_from.append(from_node)
        self._edges_to.append(to_node)
        self._edges_added.add((from_node, to_node))
        self._mutable_graph.add_edge(from_node, to_node)

        # Check if there is a deletion going from start of to_node, if so we should also have an edge past to node
        # Check if ref_pos_after comes in at a deletion, if so also add edges to after deletion
        for next_ref_pos in self._deletions[ref_pos_before_to_node+1]:
            for node_after_deletion in self._ref_pos_to_node_after[next_ref_pos - 1]:
                #logging.info("    Adding deletion edge from %d to %d" % (from_node, node_after_deletion))
                self._make_edge(from_node, node_after_deletion, next_ref_pos-1)


    def _make_node(self, ref_position_before_node, ref_position_after_node, sequence, is_ref_node=False, is_deletion=False):
        assert sequence != "", "Empty sequence for node"
        #size = len(sequence)
        #self._node_ids.append(self._current_node_id)
        #self._node_sequences.append(np.array(list(sequence)))
        #self._node_sizes.append(size)

        #if is_ref_node and len(sequence) > 0:  # We never want empty dummy nodes as reference node (by definition)
        #    self._reference_nodes.append(self._current_node_id)

        self._ref_pos_to_node_after[ref_position_before_node].append(self._current_node_id)
        self._ref_pos_to_node_start[ref_position_before_node+1].append(self._current_node_id)
        #logging.info("         Addding node %d to ref_pos_to_node_after for ref pos %d" % (self._current_node_id, ref_position_before_node))
        self._node_to_ref_pos_after[self._current_node_id].append(ref_position_after_node)
        node_id = self._current_node_id

        self._mutable_graph.add_node(node_id, sequence, is_ref_node=is_ref_node)

        self._current_node_id += 1
        return node_id

    def make_nodes(self):
        logging.info("Making nodes")
        prev_ref_node_end = -1

        i = 0
        for breakpoint_position, variant in self.breakpoints:
            i += 1

            if i % 100000 == 0:
                logging.info("%d/%d breakpoints processed" % (i, len(self.breakpoints)))

            #logging.info("At breakpoint %s, %s" % (breakpoint_position, variant))
            # Always make a ref variant from prev_ref_node_end to this breakpoints
            if breakpoint_position > prev_ref_node_end:
                #logging.info("   Making a reference node %d. Pos before/after: %d/%d." % (
                #self._current_node_id, prev_ref_node_end, breakpoint_position + 1))
                prev_ref_node_id = self._make_node(prev_ref_node_end, breakpoint_position+1, self.reference_sequence[prev_ref_node_end+1:breakpoint_position+1], is_ref_node=True)
                prev_ref_node_end = breakpoint_position
                #logging.info("    Setting prev ref node end to %d" % prev_ref_node_end)

            # If breakpoint is a new variant, make a variant node
            if variant is not None:
                if variant.type != "DELETION":
                    #logging.info("   Making a variant node %d. Pos before/after: %d/%d. Variant node sequence: %s" % (
                    #self._current_node_id, variant.get_reference_position_before_variant(),
                    #variant.get_reference_position_after_variant(), variant.get_variant_sequence()))
                    self._make_node(variant.get_reference_position_before_variant(), variant.get_reference_position_after_variant(), variant.get_variant_sequence())
                else:
                    #logging.info("Adding deletion info from ref pos %d to ref pos %d" % (variant.get_reference_position_before_variant()+1, variant.get_reference_position_after_variant()))
                    self._deletions[variant.get_reference_position_before_variant()+1].append(variant.get_reference_position_after_variant())
                    #self._node_to_ref_pos_after[prev_ref_node_id].append(variant.get_reference_position_after_variant())

        # Add one reference node for the rest of the reference sequence that is left
        self._make_node(prev_ref_node_end, len(self.reference_sequence), self.reference_sequence[prev_ref_node_end+1:len(self.reference_sequence)], is_ref_node=True)

    def make_edges(self):
        #logging.info("Ref pos to node after: %s" % self._ref_pos_to_node_after)
        #logging.info("Ref pos to node start: %s" % self._ref_pos_to_node_start)
        #logging.info("Node to ref pos after: %s" % self._node_to_ref_pos_after)
        #logging.info("Adding edges")
        i = 0
        for from_node, ref_positions_after in self._node_to_ref_pos_after.items():
            if i % 500000 == 0:
                logging.info("%d edge nodes processed" % i)
            i += 1
            for ref_pos_after in ref_positions_after:
                # Find all nodes where the ref pos before ref_pos_after goes to the node
                for to_node in self._ref_pos_to_node_after[ref_pos_after-1]:
                    if to_node == from_node:
                        continue
                    self._make_edge(from_node, to_node, ref_pos_after-1)


    def make_sorted_breakpoints(self):
        # Assume variants are already sorted
        logging.info("Traversing variants")
        for i, variant in enumerate(self.variants):
            if i % 100 == 0:
                logging.info("%d variants processed" % i)

            # Breakpoints are last base pair in a reference node
            self.breakpoints.append((variant.get_reference_position_before_variant(), variant))
            self.breakpoints.append((variant.get_reference_position_after_variant()-1, None))

        logging.info("Sorting breakpoints")
        self.breakpoints = sorted(self.breakpoints, key=lambda b: b[0])
        logging.info("Done sorting breakpoints")

        #logging.info("Breakpoints: %s" % self.breakpoints)
