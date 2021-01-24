from alignment_free_graph_genotyper.variants import GenotypeCalls
import logging
from collections import defaultdict
from .graph import Graph
import numpy as np



class GraphConstructor___:
    def __init__(self, reference_sequence, variants: GenotypeCalls):
        self.reference_sequence = reference_sequence
        self.variants = variants
        self.breakpoints = defaultdict(list)
        self.unique_breakpoint_positions = set()
        self.make_sorted_breakpoints()

        self._node_ids = []
        self._node_sizes = []
        self._node_sequences = []
        self._reference_nodes = []
        self._current_node_id = 1
        self._edges_from = []
        self._edges_to = []
        self._added_variant_nodes = defaultdict(int)
        self._variant_nodes = {}

        self.make_graph()

    def get_graph(self):
        logging.info("Node ids: %s" % self._node_ids)
        logging.info("Node sizes: %s" % self._node_sizes)
        logging.info("Node sequences: %s" % self._node_sequences)
        logging.info("Edges from: %s" % self._edges_from)
        logging.info("Edges to  : %s" % self._edges_to)
        logging.info("Reference nodes: %s" % self._reference_nodes)
        return Graph.from_flat_nodes_and_edges(np.array(self._node_ids), np.array(self._node_sequences, dtype=object), np.array(self._node_sizes), np.array(self._edges_from), np.array(self._edges_to), np.array(self._reference_nodes), np.array([1]))

    def _add_edge(self, from_node, to_node):
        self._edges_from.append(from_node)
        self._edges_to.append(to_node)

    def _add_node(self, node_sequence, is_reference_node=False):
        node_id = self._current_node_id
        self._node_ids.append(node_id)
        self._node_sizes.append(len(node_sequence))
        self._node_sequences.append(np.array(list(node_sequence)))

        if is_reference_node:
            self._reference_nodes.append(node_id)

        self._current_node_id += 1
        return node_id

    def make_graph(self):
        prev_ref_node_end = 0

        for breakpoint_position in self.unique_breakpoint_positions:
            logging.info("--")
            logging.info("Handling breakpoint position %d" % breakpoint_position)
            # Always make a ref node from previous ref_node_end
            logging.info("Making ref node from %d to %d" % (prev_ref_node_end, breakpoint_position))
            ref_node_id = self._add_node(self.reference_sequence[prev_ref_node_end:breakpoint_position], is_reference_node=True)
            prev_ref_node_end = breakpoint_position

            variants_out = [variant for type, variant in self.breakpoints[breakpoint_position] if type == "out"]
            variants_in = [variant for type, variant in self.breakpoints[breakpoint_position] if type == "in"]

            assert len(variants_out) > 0 or len(variants_in) > 0, "Found no variants out or in on breakpoint position %d" % breakpoint_position
            # Make a variant-node for all variants going out from here
            var_nodes_added = []
            for variant in variants_out:
                logging.info("  Adding variant node for variant %s with sequence %s" % (variant, variant.get_variant_sequence()))
                if variant.type == "DELETION":
                    # Don't add a variant node. The variant node is the previous ref we are going from
                    self._variant_nodes[variant.id()] = ref_node_id
                else:
                    node_id = self._add_node(variant.get_variant_sequence())
                    var_nodes_added.append(node_id)
                    self._variant_nodes[variant.id()] = node_id

            # Add edge from the ref node we just created to all variant-nodes going out from here
            for variant, variant_node_id in zip(variants_out, var_nodes_added):
                logging.info("Adding edge from previous ref node %d to variant node %d" % (ref_node_id, variant_node_id))
                self._add_edge(ref_node_id, variant_node_id)

            next_ref_node_id = self._current_node_id
            # Add edges from all variants ending at this position to the variants going out AND the next ref node
            # BUT: If a variant coming in is an insertion, we know it started after the previous base pair, and
            # we don't want to add it to insertions starting at this base pair (because then it could get an edge to itself
            # or to other insertions)
            for variant_in in variants_in:
                variant_in_node = self._variant_nodes[variant_in.id()]
                for variant_out in variants_out:
                    variant_out_node = self._variant_nodes[variant_out.id()]
                    if variant_in.type == "INSERTION" and variant_out.type == "INSERTION":
                        logging.info("SKipping edge between variants %s and %s" % (variant_in, variant_out))
                    elif variant_out.type == "DELETION":
                        # Deletions don't have variant nodes
                        logging.info("SKipping edge from variant %s to %s because the last one is deletion" % (variant_in, variant_out))
                    else:
                        logging.info("   Making edge from variant %s to variant %s" % (variant_in, variant_out))
                        self._add_edge(variant_in_node, variant_out_node)

                # All variant_in should also have an ege to the next ref node
                logging.info("Adding variant node from variant coming in %s (node=%s) to next ref node %d" % (variant_in, variant_in_node, next_ref_node_id))
                # If one variant is an insertion, the edges in should go to the dummy reference nod, not to the next reference node (except for insertions themselves, which should go to next ref node)
                self._add_edge(variant_in_node, next_ref_node_id)

            # Also always add an edge from previous ref-node to this
            self._add_edge(ref_node_id, next_ref_node_id)


        # Make one ref node at the end for the rest of the graph
        logging.info("Making ref node from %d" % (prev_ref_node_end))
        self._add_node(self.reference_sequence[prev_ref_node_end:], is_reference_node=True)


    def make_sorted_breakpoints(self):
        # Assume variants are already sorted
        logging.info("Traversing variants")
        for i, variant in enumerate(self.variants):
            if i % 100 == 0:
                logging.info("%d variants processed" % i)

            # Breakpoints are base pairs where variants went out the base pair before or variants coming in to this base pair
            self.breakpoints[variant.get_reference_position_before_variant()+1].append(("out", variant))
            self.breakpoints[variant.get_reference_position_after_variant()].append(("in", variant))
            self.unique_breakpoint_positions.add(variant.get_reference_position_before_variant()+1)
            self.unique_breakpoint_positions.add(variant.get_reference_position_after_variant())

        logging.info("Sorting breakpoints")
        self.unique_breakpoint_positions = sorted(list(self.unique_breakpoint_positions))


        logging.info("Breakpoints: %s" % self.breakpoints)
        logging.info("Breakpoints positions: %s" % self.unique_breakpoint_positions)

class GraphConstructor__:
    def __init__(self, reference_sequence, variants: GenotypeCalls):
        self.reference_sequence = reference_sequence
        self.variants = variants
        self.breakpoints = defaultdict(list)
        self.unique_breakpoint_positions = set()
        self.make_sorted_breakpoints()

        self._node_ids = []
        self._node_sizes = []
        self._node_sequences = []
        self._reference_nodes = []
        self._current_node_id = 1
        self._edges_from = []
        self._edges_to = []
        self._added_variant_nodes = defaultdict(int)
        self._variant_nodes = {}

        self.make_graph()

    def get_graph(self):
        logging.info("Node ids: %s" % self._node_ids)
        logging.info("Node sizes: %s" % self._node_sizes)
        logging.info("Node sequences: %s" % self._node_sequences)
        logging.info("Edges from: %s" % self._edges_from)
        logging.info("Edges to  : %s" % self._edges_to)
        logging.info("Reference nodes: %s" % self._reference_nodes)
        return Graph.from_flat_nodes_and_edges(np.array(self._node_ids), np.array(self._node_sequences, dtype=object), np.array(self._node_sizes), np.array(self._edges_from), np.array(self._edges_to), np.array(self._reference_nodes), np.array([1]))

    def _add_edge(self, from_node, to_node):
        self._edges_from.append(from_node)
        self._edges_to.append(to_node)

    def _add_node(self, node_sequence, is_reference_node=False):
        node_id = self._current_node_id
        self._node_ids.append(node_id)
        self._node_sizes.append(len(node_sequence))
        self._node_sequences.append(np.array(list(node_sequence)))

        if is_reference_node:
            self._reference_nodes.append(node_id)

        self._current_node_id += 1
        return node_id

    def make_graph(self):
        prev_ref_node_end = 0

        for breakpoint_position in self.unique_breakpoint_positions:
            logging.info("--")
            logging.info("Handling breakpoint position %d" % breakpoint_position)
            # Always make a ref node from previous ref_node_end
            logging.info("Making ref node from %d to %d" % (prev_ref_node_end, breakpoint_position))
            ref_node_id = self._add_node(self.reference_sequence[prev_ref_node_end:breakpoint_position], is_reference_node=True)
            prev_ref_node_end = breakpoint_position

            variants_out = [variant for type, variant in self.breakpoints[breakpoint_position] if type == "out"]
            variants_in = [variant for type, variant in self.breakpoints[breakpoint_position] if type == "in"]

            assert len(variants_out) > 0 or len(variants_in) > 0, "Found no variants out or in on breakpoint position %d" % breakpoint_position
            # Make a variant-node for all variants going out from here
            var_nodes_added = []
            insertion_exists = False
            for variant in variants_out:
                logging.info("  Adding variant node for variant %s with sequence %s" % (variant, variant.get_variant_sequence()))
                node_id = self._add_node(variant.get_variant_sequence())
                var_nodes_added.append(node_id)
                self._variant_nodes[variant.id()] = node_id
                if variant.type == "INSERTION":
                    insertion_exists = True

            # Add edge from the ref node we just created to all variant-nodes going out from here
            for variant, variant_node_id in zip(variants_out, var_nodes_added):
                logging.info("Adding edge from previous ref node %d to variant node %d" % (ref_node_id, variant_node_id))
                self._add_edge(ref_node_id, variant_node_id)

            # If any of the variants was an insertion, we need to add one insertion dummy node on the reference
            # All variants coming in on this event should then go to this dummy node instead of the next ref node after the insertion
            dummy_reference_node = None
            if insertion_exists:
                dummy_reference_node = self._current_node_id
                logging.info("!! An insertion exists. Add dummy node with id %d" % dummy_reference_node)
                self._add_edge(ref_node_id, dummy_reference_node)
                self._add_node("")

            next_ref_node_id = self._current_node_id
            # Add edges from all variants ending at this position to the variants going out AND the next ref node
            # BUT: If a variant coming in is an insertion, we know it started after the previous base pair, and
            # we don't want to add it to insertions starting at this base pair (because then it could get an edge to itself
            # or to other insertions)
            for variant_in in variants_in:
                variant_in_node = self._variant_nodes[variant_in.id()]
                for variant_out in variants_out:
                    variant_out_node = self._variant_nodes[variant_out.id()]
                    if variant_in.type == "INSERTION" and variant_out.type == "INSERTION":
                        logging.info("SKipping edge between variants %s and %s" % (variant_in, variant_out))
                    else:
                        logging.info("   Making edge from variant %s to variant %s" % (variant_in, variant_out))
                        self._add_edge(variant_in_node, variant_out_node)

                # All variant_in should also have an ege to the next ref node
                logging.info("Adding variant node from variant coming in %s to next ref node %d" % (variant_in_node, next_ref_node_id))
                # If one variant is an insertion, the edges in should go to the dummy reference nod, not to the next reference node (except for insertions themselves, which should go to next ref node)
                if insertion_exists and variant_in.type != "INSERTION":
                    self._add_edge(variant_in_node, dummy_reference_node)
                else:
                    self._add_edge(variant_in_node, next_ref_node_id)

            # Also always add an edge from previous ref-node to this
            if dummy_reference_node is not None:
                self._add_edge(dummy_reference_node, next_ref_node_id)
            else:
                self._add_edge(ref_node_id, next_ref_node_id)


        # Make one ref node at the end for the rest of the graph
        logging.info("Making ref node from %d" % (prev_ref_node_end))
        self._add_node(self.reference_sequence[prev_ref_node_end:], is_reference_node=True)


    def make_sorted_breakpoints(self):
        # Assume variants are already sorted
        logging.info("Traversing variants")
        for i, variant in enumerate(self.variants):
            if i % 100 == 0:
                logging.info("%d variants processed" % i)

            # Breakpoints are base pairs where variants went out the base pair before or variants coming in to this base pair
            self.breakpoints[variant.get_reference_position_before_variant()+1].append(("out", variant))
            self.breakpoints[variant.get_reference_position_after_variant()].append(("in", variant))
            self.unique_breakpoint_positions.add(variant.get_reference_position_before_variant()+1)
            self.unique_breakpoint_positions.add(variant.get_reference_position_after_variant())

        logging.info("Sorting breakpoints")
        self.unique_breakpoint_positions = sorted(list(self.unique_breakpoint_positions))


        logging.info("Breakpoints: %s" % self.breakpoints)
        logging.info("Breakpoints positions: %s" % self.unique_breakpoint_positions)

class GraphConstructor:
    def __init__(self, reference_sequence, variants: GenotypeCalls):
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
        self._deletions = defaultdict(list)


        self.make_nodes()
        self.make_edges()

    def get_graph(self):
        return Graph.from_flat_nodes_and_edges(np.array(self._node_ids), np.array(self._node_sequences, dtype=object), np.array(self._node_sizes), np.array(self._edges_from), np.array(self._edges_to), np.array(self._reference_nodes), np.array([1]))

    def _make_edge(self, from_node, to_node, ref_pos_before_to_node):
        self._edges_from.append(from_node)
        self._edges_to.append(to_node)

        # Check if there is a deletion going from start of to_node, if so we should also have an edge past to node
        # Check if ref_pos_after comes in at a deletion, if so also add edges to after deletion
        for next_ref_pos in self._deletions[ref_pos_before_to_node+1]:
            for node_after_deletion in self._ref_pos_to_node_after[next_ref_pos - 1]:
                logging.info("    Adding deletion edge from %d to %d" % (from_node, node_after_deletion))
                self._make_edge(from_node, node_after_deletion, next_ref_pos-1)


    def _make_node(self, ref_position_before_node, ref_position_after_node, sequence, is_ref_node=False, is_deletion=False):
        assert sequence != "", "Empty sequence for node"
        size = len(sequence)
        self._node_ids.append(self._current_node_id)
        self._node_sequences.append(np.array(list(sequence)))
        self._node_sizes.append(size)

        if is_ref_node and len(sequence) > 0:  # We never want empty dummy nodes as reference node (by definition)
            self._reference_nodes.append(self._current_node_id)

        self._ref_pos_to_node_after[ref_position_before_node].append(self._current_node_id)
        self._ref_pos_to_node_start[ref_position_before_node+1].append(self._current_node_id)
        logging.info("         Addding node %d to ref_pos_to_node_after for ref pos %d" % (self._current_node_id, ref_position_before_node))
        self._node_to_ref_pos_after[self._current_node_id].append(ref_position_after_node)
        node_id = self._current_node_id
        self._current_node_id += 1
        return node_id

    def make_nodes(self):
        prev_ref_node_end = -1
        prev_ref_node_id = None

        for breakpoint_position, variant in self.breakpoints:
            logging.info("At breakpoint %s, %s" % (breakpoint_position, variant))
            # Always make a ref variant from prev_ref_node_end to this breakpoints
            if breakpoint_position > prev_ref_node_end:
                logging.info("   Making a reference node %d. Pos before/after: %d/%d." % (
                    self._current_node_id, prev_ref_node_end, breakpoint_position + 1))
                prev_ref_node_id = self._make_node(prev_ref_node_end, breakpoint_position+1, self.reference_sequence[prev_ref_node_end+1:breakpoint_position+1], is_ref_node=True)
                prev_ref_node_end = breakpoint_position
                logging.info("    Setting prev ref node end to %d" % prev_ref_node_end)

            # If breakpoint is a new variant, make a variant node
            if variant is not None:
                if variant.type != "DELETION":
                    logging.info("   Making a variant node %d. Pos before/after: %d/%d. Variant node sequence: %s" % (
                        self._current_node_id, variant.get_reference_position_before_variant(),
                        variant.get_reference_position_after_variant(), variant.get_variant_sequence()))
                    self._make_node(variant.get_reference_position_before_variant(), variant.get_reference_position_after_variant(), variant.get_variant_sequence())
                else:
                    logging.info("Adding deletion info from ref pos %d to ref pos %d" % (variant.get_reference_position_before_variant()+1, variant.get_reference_position_after_variant()))
                    self._deletions[variant.get_reference_position_before_variant()+1].append(variant.get_reference_position_after_variant())
                    #self._node_to_ref_pos_after[prev_ref_node_id].append(variant.get_reference_position_after_variant())

        # Add one reference node for the rest of the reference sequence that is left
        self._make_node(prev_ref_node_end, len(self.reference_sequence), self.reference_sequence[prev_ref_node_end+1:len(self.reference_sequence)], is_ref_node=True)

    def make_edges(self):
        logging.info("Ref pos to node after: %s" % self._ref_pos_to_node_after)
        logging.info("Ref pos to node start: %s" % self._ref_pos_to_node_start)
        logging.info("Node to ref pos after: %s" % self._node_to_ref_pos_after)
        for from_node, ref_positions_after in self._node_to_ref_pos_after.items():
            for ref_pos_after in ref_positions_after:
                # Find all nodes where the ref pos before ref_pos_after goes to the node
                for to_node in self._ref_pos_to_node_after[ref_pos_after-1]:
                    if to_node == from_node:
                        continue
                    logging.info("Making edge from %d to %d" % (from_node, to_node))
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

        logging.info("Breakpoints: %s" % self.breakpoints)