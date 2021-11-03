import itertools
import json
import logging
import re
import numpy as np
from .variants import VcfVariants
from collections import defaultdict
from .mutable_graph import MutableGraph

class VariantNotFoundException(Exception):
    pass

class Graph:
    properties = {
        "nodes", "node_to_sequence_index", "node_sequences", "node_to_edge_index", "node_to_n_edges", "edges",
        "node_to_ref_offset", "ref_offset_to_node", "chromosome_start_nodes", "linear_ref_nodes_index", "allele_frequencies", "numeric_node_sequences"
    }

    def __init__(self, nodes=None, node_to_sequence_index=None, node_sequences=None, node_to_edge_index=None, node_to_n_edges=None, edges=None,
                 node_to_ref_offset=None, ref_offset_to_node=None, chromosome_start_nodes=None, allele_frequencies=None, linear_ref_nodes_index=None,
                 numeric_node_sequences=None, linear_ref_nodes_and_dummy_nodes_index=None):
        self.nodes = nodes
        self.node_to_sequence_index = node_to_sequence_index
        self.node_sequences = node_sequences
        self.node_to_edge_index = node_to_edge_index
        self.node_to_n_edges = node_to_n_edges
        self.edges = edges
        self.node_to_ref_offset = node_to_ref_offset
        self.ref_offset_to_node = ref_offset_to_node
        self._linear_ref_nodes_cache = None
        self.chromosome_start_nodes = chromosome_start_nodes
        self.allele_frequencies = allele_frequencies
        self.numeric_node_sequences = numeric_node_sequences

        if self.nodes is not None and linear_ref_nodes_index is None:
            logging.info("Makng a linear ref lookup, since it is not provided")
            self.linear_ref_nodes_index = np.zeros(len(self.nodes), dtype=np.uint8)
            self.linear_ref_nodes_index[np.array(list(self.linear_ref_nodes()))] = 1
            logging.info("Done making linear ref nodes index")
        else:
            self.linear_ref_nodes_index = linear_ref_nodes_index
        
        if self.nodes is not None and linear_ref_nodes_and_dummy_nodes_index is None:
            self.linear_ref_nodes_and_dummy_nodes_index = None
            logging.info("Making linear ref nodes and dummy nodes index")
            self.make_linear_ref_node_and_ref_dummy_node_index()
        else:
            self.linear_ref_nodes_and_dummy_nodes_index = linear_ref_nodes_and_dummy_nodes_index

    def node_has_edges(self, node, edges):
        e = self.get_edges(node)
        if len(edges) != len(e):
            return False

        for edge in edges:
            if edge not in edges:
                return False

        return True

    def get_all_nodes(self):
        return np.union1d(np.where(self.nodes != 0)[0], np.where(self.node_to_n_edges > 0)[0])
        # This does not return nodes of size 0
        #return np.where(self.nodes != 0)[0]

    def __str__(self):
        mutable_graph = self.to_mutable_graph()
        return str(mutable_graph)

    def __repr__(self):
        return self.__str__()

    def to_mutable_graph(self):
        all_nodes = self.get_all_nodes()
        nodes = {node: self.get_node_size(node) for node in all_nodes}
        node_sequences = {node: self.get_node_sequence(node) for node in all_nodes}
        edges = {node: list(self.get_edges(node)) for node in all_nodes}
        return MutableGraph(nodes, node_sequences, edges, list(self.linear_ref_nodes()), self.node_to_ref_offset, self.ref_offset_to_node, self.chromosome_start_nodes, self.allele_frequencies)

    @classmethod
    def from_mutable_graph(cls, mutable_graph):
        return cls.from_dicts(mutable_graph.node_sequences, mutable_graph.edges, mutable_graph.linear_ref_nodes)

    def get_node_size(self, node):
        return self.nodes[node]

    def get_node_allele_frequency(self, node):
        if self.allele_frequencies is None:
            return 1.0
        return self.allele_frequencies[node]

    def get_node_subsequence(self, node, start, end):
        index_position = self.node_to_sequence_index[node]
        return ''.join(self.node_sequences[int(index_position+start):int(index_position+end)])

    def get_node_sequence(self, node):
        index_position = self.node_to_sequence_index[node]
        return ''.join(self.node_sequences[index_position:index_position+self.nodes[node]])

    def get_nodes_sequence_index_positions(self, nodes):
        # First find the number of indexes we ned
        sequence_length = np.sum(self.nodes[nodes])
        indexes = np.zeros(sequence_length, dtype=np.uint64)
        i = 0
        for node in nodes:
            node_size = self.nodes[node]
            index_position = self.node_to_sequence_index[node]
            indexes[i:i+node_size] = np.arange(index_position, index_position+node_size)
            i += self.nodes[node]

        return indexes

    def max_node_id(self):
        return len(self.nodes)-1

    def get_numeric_node_sequences(self, nodes):
        logging.info("Getting sequence length")
        sequence_length = np.sum(self.nodes[nodes])
        numeric_sequence = np.zeros(sequence_length, dtype=np.uint8)

        i = 0
        for node in nodes:
            node_size = self.nodes[node]
            index_position = self.node_to_sequence_index[node]
            numeric_sequence[i:i + node_size] = self.numeric_node_sequences[index_position:index_position+node_size]
            i += self.nodes[node]

        return numeric_sequence

    def get_nodes_sequences2(self, nodes):
        return ''.join(self.node_sequences[self.get_nodes_sequence_index_positions(nodes)])

    def get_nodes_sequence(self, nodes):
        sequences = []

        for node in nodes:
            index_position = self.node_to_sequence_index[node]
            sequences.extend(self.node_sequences[index_position:index_position+self.nodes[node]])

        #for node in nodes:
        #    sequences.append(self.get_node_sequence(node))
        return ''.join(sequences)

    def get_node_offset_at_chromosome_and_chromosome_offset(self, chromosome, offset):
        chromosome_position = chromosome - 1
        chromosome_start_node = self.chromosome_start_nodes[chromosome_position]
        chromosome_offset = self.get_ref_offset_at_node(chromosome_start_node)
        real_offset = int(chromosome_offset + offset)
        node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, offset)
        node_offset = self.get_ref_offset_at_node(node)
        return real_offset - node_offset

    def get_node_at_chromosome_and_chromosome_offset(self, chromosome, offset):
        chromosome_position = chromosome - 1
        try:
            chromosome_start_node = self.chromosome_start_nodes[chromosome_position]
        except IndexError:
            logging.error("Could not find chromosome start position for chromosome %d. Chromosome start nodes are %s" % (chromosome, self.chromosome_start_nodes))
            raise


        # Shift offset with chromosome start offset
        chromosome_offset = self.get_ref_offset_at_node(chromosome_start_node)
        real_offset = int(chromosome_offset + offset)
        return self.get_node_at_ref_offset(real_offset)

    def get_edges(self, node):
        if node >= len(self.node_to_edge_index):
            return []

        index = self.node_to_edge_index[node]
        n_edges = self.node_to_n_edges[node]

        if n_edges == 0:
            return []

        return self.edges[index:index+n_edges]

    def linear_ref_nodes(self):
        if self._linear_ref_nodes_cache is not None:
            return self._linear_ref_nodes_cache
        else:
            nodes = set(np.unique(self.ref_offset_to_node))
            self._linear_ref_nodes_cache = nodes
            return nodes

    def make_linear_ref_node_and_ref_dummy_node_index(self):
        linear_ref_nodes_and_dummy_nodes_index = np.zeros(len(self.linear_ref_nodes_index), dtype=np.uint8)
        for node in range(len(linear_ref_nodes_and_dummy_nodes_index)):
            if node % 100000 == 0:
                logging.info("%d nodes processed" % node)

            if self.is_linear_ref_node_or_linear_ref_dummy_node(node):
                linear_ref_nodes_and_dummy_nodes_index[node] = 1

        self.linear_ref_nodes_and_dummy_nodes_index = linear_ref_nodes_and_dummy_nodes_index

    def is_linear_ref_node_or_linear_ref_dummy_node(self, node):
        if self.linear_ref_nodes_and_dummy_nodes_index is not None:
            if self.linear_ref_nodes_and_dummy_nodes_index[node] == 1:
                return True
            else:
                return False

        if self.is_linear_ref_node(node):
            return True

        if self.get_node_size(node) == 0:
            # is a linear if it has a linear ref node as next node and the previous node on the linear ref from this linear ref node has an edge to node
            if len([next_node for next_node in self.get_edges(node) if self.is_linear_ref_node(next_node) and node in self.get_edges(self.get_node_at_ref_offset(int(self.get_ref_offset_at_node(next_node)-1)))]) > 0:
                return True

        return False

    def is_linear_ref_node(self, node):
        if self.linear_ref_nodes_index[node] != 0:
            return True
        return False

    def linear_ref_length(self):
        return len(self.ref_offset_to_node)

    def get_node_at_ref_offset(self, ref_offset):
        try:
            return self.ref_offset_to_node[ref_offset]
        except IndexError:
            logging.error("Invalid ref offset %s" % str(ref_offset))
            raise

    def get_ref_offset_at_node(self, node):
        return self.node_to_ref_offset[node]

    def get_node_offset_at_ref_offset(self, ref_offset):
        node = self.get_node_at_ref_offset(ref_offset)
        offset_at_node = self.get_ref_offset_at_node(node)
        return ref_offset - offset_at_node

    def to_file(self, file_name):
        allele_frequencies = self.allele_frequencies
        if allele_frequencies is None:
            allele_frequencies = np.zeros(0)

        numeric_node_sequences = self.numeric_node_sequences
        if numeric_node_sequences is None:
            numeric_node_sequences = np.zeros(0)

        np.savez(file_name,
                 nodes=self.nodes,
                 node_to_sequence_index=self.node_to_sequence_index,
                 node_sequences=self.node_sequences,
                 node_to_edge_index=self.node_to_edge_index,
                 node_to_n_edges=self.node_to_n_edges,
                 edges=self.edges,
                 node_to_ref_offset=self.node_to_ref_offset,
                 ref_offset_to_node=self.ref_offset_to_node,
                 chromosome_start_nodes=self.chromosome_start_nodes,
                 allele_frequencies=allele_frequencies,
                 linear_ref_nodes_index=self.linear_ref_nodes_index,
                 numeric_node_sequences=numeric_node_sequences,
                 linear_ref_nodes_and_dummy_nodes_index=self.linear_ref_nodes_and_dummy_nodes_index
                 )

        logging.info("Saved to file %s" % file_name)

    @classmethod
    def from_file(cls, file_name):
        logging.info("Reading file from %s" % file_name)
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npz")

        allele_frequencies = None
        if "allele_frequencies" in data:
            allele_frequencies = data["allele_frequencies"]

        linear_ref_nodes_index = None
        if "linear_ref_nodes_index" in data:
            linear_ref_nodes_index = data["linear_ref_nodes_index"]

        linear_ref_nodes_and_dummy_nodes_index = None
        if "linear_ref_nodes_and_dummy_nodes_index" in data:
            linear_ref_nodes_and_dummy_nodes_index = data["linear_ref_nodes_and_dummy_nodes_index"]

        numeric_node_sequences = None
        if "numeric_node_sequences" in data:
            numeric_node_sequences = data["numeric_node_sequences"]

        return cls(data["nodes"],
                   data["node_to_sequence_index"],
                   data["node_sequences"],
                   data["node_to_edge_index"],
                   data["node_to_n_edges"],
                   data["edges"],
                   data["node_to_ref_offset"],
                   data["ref_offset_to_node"],
                   data["chromosome_start_nodes"],
                   allele_frequencies,
                   linear_ref_nodes_index,
                   numeric_node_sequences,
                   linear_ref_nodes_and_dummy_nodes_index)


    def get_flat_graph(self):
        node_ids = list(np.where(self.nodes > 0)[0])
        node_sizes = list(self.nodes[node_ids])
        node_sequences = []
        for node in node_ids:
            node_sequences.append(list(self.get_node_sequence(node)))

        linear_ref_nodes = list(self.linear_ref_nodes())
        from_nodes = []
        to_nodes = []
        for from_node in node_ids:
            for to_node in self.get_edges(from_node):
                from_nodes.append(from_node)
                to_nodes.append(to_node)

        return node_ids, node_sequences, node_sizes, from_nodes, to_nodes, linear_ref_nodes, self.chromosome_start_nodes

    @classmethod
    def from_dicts(cls, node_sequences, edges, linear_ref_nodes):
        assert linear_ref_nodes is not None
        logging.info("Making graph from dicts")
        nodes = list(node_sequences.keys())
        logging.info("Making sequences")
        node_sequences = [list(node_sequences[node]) for node in nodes]
        node_sizes = [len(seq) for seq in node_sequences]
        from_nodes = []
        to_nodes = []
        i = 0
        for from_node, to_nodes_edges in edges.items():
            if i % 100000 == 0:
                logging.info("%d nodes processed" % i)
            for to_node in to_nodes_edges:
                from_nodes.append(from_node)
                to_nodes.append(to_node)

            i += 1

        to_nodes_set = set(to_nodes)
        logging.info("Done preparing data from dicts")
        logging.info("There are %d node sequences" % len(node_sequences))
        return Graph.from_flat_nodes_and_edges(nodes, node_sequences, node_sizes, np.array(from_nodes), np.array(to_nodes), np.array(linear_ref_nodes), [node for node in nodes if node not in to_nodes_set])

    @classmethod
    def from_flat_nodes_and_edges(cls, node_ids, node_sequences, node_sizes, from_nodes, to_nodes, linear_ref_nodes, chromosome_start_nodes):
        logging.info("Chromosome start nodes are: %s" % chromosome_start_nodes)

        """
        logging.info("Asserting linear ref nodes are not empty")
        for i, node in enumerate(node_ids):
            size = node_sizes[i]
            if size == 0 and node in linear_ref_nodes:
                raise Exception("Placeholder nodes (for insertions) cannot be marked as linear ref nodes, since that breaks ref position to node lookup. Best solution is to not put these in linear ref nodes")
        """


        max_node = np.max(node_ids)
        nodes = np.zeros(max_node+1, dtype=np.uint32)
        nodes[node_ids] = node_sizes

        # Node sequences
        logging.info("Making node sequences")
        new_node_positions = np.cumsum(node_sizes)
        node_sequence_index = np.zeros(max_node+1, dtype=np.uint64)
        node_sequence_index[node_ids[0]] = 0
        node_sequence_index[node_ids[1:]] = new_node_positions[:-1]
        logging.info("Merging sequences into one string, chaining them first")
        node_sequences = list(itertools.chain.from_iterable(node_sequences))
        logging.info("Total length of graph sequences: %d" % len(node_sequences))
        logging.info("Done chaining sequences")
        node_sequences_array = np.array(node_sequences).astype("<U1")
        logging.info("Done merging sequences")


        #node_sequences_array = np.zeros(max_node+0, dtype=object)
        #logging.info("Allocating empty sequence array for node sequences")
        #node_sequences_array = np.array([np.zeros(size, dtype="<U1") for size in nodes])  # Must allocate space
        #node_sequences_array[node_ids] = node_sequences

        logging.info("Sorting nodes")
        sorting = np.argsort(from_nodes)
        from_nodes = from_nodes[sorting]
        to_nodes = to_nodes[sorting]

        logging.info("Making node index")
        diffs = np.ediff1d(from_nodes, to_begin=1)
        positions_of_unique_nodes = np.nonzero(diffs)[0]
        unique_nodes = from_nodes[positions_of_unique_nodes]

        logging.info("Making node to edge index")
        node_to_edge_index = np.zeros(max_node+1, dtype=np.uint32)
        node_to_n_edges = np.zeros(max_node+1, dtype=np.uint8)
        node_to_edge_index[unique_nodes] = positions_of_unique_nodes
        n_edges_numbers = np.ediff1d(positions_of_unique_nodes, to_end=len(from_nodes)-positions_of_unique_nodes[-1])
        node_to_n_edges[unique_nodes] = n_edges_numbers

        #logging.info("Finding ref offsets")
        node_to_ref_offset = np.zeros(max_node+1, np.uint64)
        try:
            ref_node_sizes = nodes[linear_ref_nodes]
        except IndexError:
            logging.error("Problem with linear ref nodes, which are %s" % linear_ref_nodes)
            print(type(linear_ref_nodes))
            raise
        #print("Ref node sizes: %s"  % ref_node_sizes)
        ref_offsets = np.cumsum(ref_node_sizes)
        #print("Ref offsets: %s"  % ref_offsets)
        node_to_ref_offset[linear_ref_nodes[1:]] = ref_offsets[:-1]

        # Find last node to add to linear ref size
        last_node_in_graph = np.argmax(node_to_ref_offset)
        last_node_size = nodes[last_node_in_graph]
        #logging.info("Last node in graph is %d with size %d" % (last_node_in_graph, last_node_size))

        ref_offset_to_node = np.zeros(int(np.max(node_to_ref_offset)) + last_node_size)
        index_positions = np.cumsum(ref_node_sizes)[:-1]
        #logging.info("INdex positions: %s" % index_positions)
        #logging.info("Linear ref nodes: %s" % linear_ref_nodes)
        #logging.info("Diff linear ref nodes: %s" % np.diff(linear_ref_nodes))
        ref_offset_to_node[index_positions] = np.diff(linear_ref_nodes)
        ref_offset_to_node[0] = linear_ref_nodes[0]
        #logging.info("Ref offset to node 1: %s" % ref_offset_to_node)
        ref_offset_to_node = np.cumsum(ref_offset_to_node, dtype=np.uint32)

        return cls(nodes, node_sequence_index, node_sequences_array, node_to_edge_index,
                   node_to_n_edges, to_nodes, node_to_ref_offset, ref_offset_to_node, chromosome_start_nodes)

    @classmethod
    def from_vg_json_files(cls, file_names):

        node_sequences = []
        node_ids = []
        node_sizes = []
        edges_from = []
        edges_to = []
        linear_ref_nodes = []
        chromosomes = []
        chromosome_start_nodes = []
        n_nodes_added = 0
        n_edges_added = 0
        chromosome_offset = 0
        for file_name in file_names:
            # Hackish (send in chromosome on comand line later):
            path_found = False
            file = open(file_name)
            for line in file:
                json_object = json.loads(line)
                if "node" in json_object:
                    for node in json_object["node"]:
                        id = int(node["id"])
                        node_ids.append(id)
                        node_sizes.append(len(node["sequence"]))
                        node_sequence = list(re.sub(r'[^acgtn]', "n", node["sequence"].lower()))
                        node_sequences.append(node_sequence)
                        n_nodes_added += 1
                        if n_nodes_added % 100000 == 0:
                            logging.info("%d nodes added" % n_nodes_added)

                if "edge" in json_object:
                    for edge in json_object["edge"]:
                        from_node = int(edge["from"])
                        to_node = int(edge["to"])

                        edges_from.append(from_node)
                        edges_to.append(to_node)
                        n_edges_added += 1
                        if n_edges_added % 100000 == 0:
                            logging.info("%d edges added" % n_edges_added)

                if "path" in json_object:
                    for path in json_object["path"]:
                        assert not path_found, "Found multiple paths, not sure which is the reference path. Path now:" % path["name"]
                        logging.info("Found path %s, assuming this is the reference path" % path["name"])
                        nodes_in_path = [mapping["position"]["node_id"] for mapping in path["mapping"]]
                        linear_ref_nodes.extend(nodes_in_path)
                        chromosome_start_nodes.append(nodes_in_path[0])
                        path_found = True

                logging.info("Chromosome start nodes are: %s" % chromosome_start_nodes)

        node_ids = np.array(node_ids)
        node_sizes = np.array(node_sizes)
        edges_from = np.array(edges_from)
        edges_to = np.array(edges_to)
        linear_ref_nodes = np.array(linear_ref_nodes, dtype=np.uint32)
        chromosome_start_nodes = np.array(chromosome_start_nodes, dtype=np.uint32)
        logging.info("Chromosome start nodes are: %s" % chromosome_start_nodes)

        return cls.from_flat_nodes_and_edges(node_ids, node_sequences, node_sizes, edges_from, edges_to, linear_ref_nodes, chromosome_start_nodes)


    def get_snp_nodes(self, ref_offset, variant_bases, chromosome=1):
        node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        node_offset = self.get_node_offset_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        assert node_offset == 0, "Node offset is 0 at ref_offset %d. Variant bases: %s. Chromosome %d" % (ref_offset, variant_bases, chromosome)
        prev_node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset - 1)

        _, possible_snp_nodes = self.find_nodes_from_node_that_matches_sequence(prev_node, variant_bases, [], [])
        assert len(possible_snp_nodes) > 0, "Did not find any possible snp nodes for variant at ref offset %s with variant sequence %s" % (ref_offset, variant_bases)

        for possible_snp_node in possible_snp_nodes:
            # Not true if deletion before snp
            #assert len(possible_snp_node) == 1, "SNP nodes should only be one node"
            assert len([n for n in possible_snp_node if self.get_node_size(n) > 0]) == 1, "There should only be one non-empty SNP node"
            potential_next = possible_snp_node[-1]  # Get last node, this will be the snp node if there are more nodes
            assert self.get_node_size(potential_next) > 0
            # Also require that this node shares next node with our ref snp node, if not this could be a false match against an indel node
            next_from_snp_node = self.get_edges(node)
            if len([n for n in next_from_snp_node if n in self.get_edges(potential_next)]) > 0:
                return node, potential_next

        logging.error("Could not parse substitution at offset %d with bases %s" % (ref_offset, variant_bases))
        logging.error("Next nodes from snp node: %d" % next_from_snp_node)
        raise Exception("Parseerrror")

    def get_deletion_nodes(self, ref_offset, deletion_length, deleted_sequence, chromosome=1):
        ref_offset += 1
        # Node is the reference node that is deleted (the first one if there are multiple)
        node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        node_offset = self.get_node_offset_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        #logging.info("Processing deletion at ref pos %d with size %d. Node inside deletion: %d" % (ref_offset, deletion_length, node))

        assert node_offset == 0, "Node offset is %d, not 0" %  (node_offset)

        prev_node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset - 1)
        #logging.info("Node before deletion: %d" % prev_node)

        # Find next reference node with offset corresponding to the number of deleted base pairs
        next_ref_pos = ref_offset + deletion_length
        next_ref_node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset + deletion_length)
        #logging.info("Node after deletion: %d" % next_ref_node)
        if self.get_node_offset_at_chromosome_and_chromosome_offset(chromosome, next_ref_pos) != 0:
            logging.info("Could not find deletion at ref offset %d in graph" % ref_offset)
            #logging.error("Offset %d is not at beginning of node" % next_ref_pos)
            #logging.error("Node at %d: %d" % (next_ref_pos, next_ref_node))
            #logging.error("Ref length in deletion: %s" % deletion_length)
            #logging.info("Ref pos beginning of deletion: %d" % ref_offset)

            raise VariantNotFoundException("Deletion not in graph")

        # Find an empty node between prev_node and next_ref_node
        # This solution should handle cases where the reference has multiple nodes inside deleted segment (e.g. a SNP inside a deleted sequence)
        # This is not true anymore, since a deltion can be followed by an insertion, leading to two or more dummy nodes before reaching the linear ref again
        #deletion_nodes = [node for node in self.get_edges(prev_node) if next_ref_node in self.get_edges(node) and self.get_node_size(node) == 0]
        deletion_nodes = [node for node in self.get_edges(prev_node) if self.get_node_size(node) == 0]
        #debug = [(node, self.get_edges(node), self.get_node_size(node)) for node in self.get_edges(prev_node)]
        #logging.info("%s" % debug)




        if len(deletion_nodes) != 1:

            # There can be more than one deletion node if there is an insertion right after a deletion
            if len(deletion_nodes) > 1:
                # If there are more than 1 deletion node, we find all paths that matches the deleted sequence, and try to find the deletion node that starts before and ends after this path
                one_path, possible_deleted_paths = self.find_nodes_from_node_that_matches_sequence(prev_node, deleted_sequence, [], [])
                assert len(one_path) >  0 and len(possible_deleted_paths) > 0, "Found no paths from node %d that matches the deleted sequence %s" % (prev_node, deleted_sequence)
                deletion_nodes_orig = deletion_nodes.copy()

                #print(possible_deleted_paths)
                deletion_nodes = [n for n in deletion_nodes if len([path for path in possible_deleted_paths if len([edge_from_deletion_node for edge_from_deletion_node in self.get_edges(n) if edge_from_deletion_node in self.get_edges(path[-1])]) > 0]) > 0]

                if len(deletion_nodes) == 0:
                    logging.warning("Error on deletion at pos %s. Node: %d. Deletion nodes: %s. Orig deletion nodes: %s. Edges going out: %s" % (ref_offset, node, deletion_nodes, deletion_nodes_orig, self.get_edges(node)))
                    logging.error("Could not find any nodes among possible deletion nodes %s that wrapped around any of the paths %s" % (deletion_nodes_orig, possible_deleted_paths))
                    raise VariantNotFoundException("Deletion not found")
                elif len(deletion_nodes) > 1:
                    logging.warning("Deletion at ref offset %s has multiple empty nodes going out from %d: %s. They all wrap around a deleted path in %s" % (ref_offset, prev_node, deletion_nodes, possible_deleted_paths))
                    if deletion_nodes[0] not in self.get_edges(node):
                        logging.warning("Error on deletion. Deletion nodes: %s. Edges going out: %s" % (deletion_nodes, self.get_edges(node)))
                        raise VariantNotFoundException("Deletion not found")
                    # Assert that there is an empty node between this deletion node and next ref node
                    # This assertion won't always hold, e.g. if there are two insertions after a deltion, so just skip it
                    #assert len([n for n in self.get_edges(deletion_nodes[0]) if next_ref_node in self.get_edges(n) and self.get_node_size(n) == 0]), "There is more than one deletion node, and the next node
            else:
                logging.warning("Could not find deletion at ref offset %d in graph" % ref_offset)
                logging.warning("There should be only one deletion node between %d and %d for variant at ref pos %d. There are %d. Edges out from %d are: %s. Edges out from %d are %s" % (prev_node, next_ref_node, ref_offset, len(deletion_nodes), prev_node, self.get_edges(prev_node), node, self.get_edges(node)))
                raise VariantNotFoundException("Deletion not in graph")

        deletion_node = deletion_nodes[0]
        # Not true anymore, e.g. if there is a snp right after a deletion
        #assert len(self.get_edges(deletion_node)) == 1, "Deletion node %d has not exactly 1 edge. Variant at ref offset %d. Edges are %s" % (deletion_node, ref_offset, self.get_edges(deletion_node))
        return (node, deletion_node)

    def get_insertion_nodes(self, variant, chromosome=1):
        ref_offset = variant.position
        #logging.info("Ref offset: %d" % ref_offset)
        # Find node right before insertion node
        node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset-1)
        #logging.info("Node at ref offset: %d" % node)
        node_offset = self.get_node_offset_at_chromosome_and_chromosome_offset(chromosome, ref_offset-1)
        node_size = self.get_node_size(node)
        insertion_length = len(variant.variant_sequence) - 1
        insertion_sequence = variant.variant_sequence[1:]  # First base is not included

        assert node_offset == node_size - 1, "Node offset %d is not at end of node %d which has size %d. Insertion not found in graph." % (node_offset, node, node_size)

        # Find out which next node matches the insertion
        # NB! There might be more than one node with the same sequence going out before the variant, since a SNP can have same sequence as the insertion
        next_ref_node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        variant_node = None
        for potential_next in self.get_edges(node):
            #if potential_next in self.linear_ref_nodes():
            if self.is_linear_ref_node(potential_next):
                continue  # Next should not be in linear ref

            #print("Processing insertion at ref offset %d with base %s" % (ref_offset, base))
            #print("  Node %d with offset %d. Node size: %d" % (node, node_offset, self.blocks[node].length()))

            if self.get_node_size(potential_next) == 0:
                continue  # This is an empty deletion node

            #next_bases = self.get_node_sequence(potential_next)[0:insertion_length].upper()
            next_bases = self.get_node_subsequence(potential_next, 0, insertion_length).upper()
            if next_bases == insertion_sequence.upper():
                if next_ref_node in self.get_edges(potential_next):
                    variant_node = potential_next
                    break

        assert variant_node is not None, "Could not find insertion node. No nodes going out from %d has sequence %s" % (node, variant.variant_sequence)

        # Find linear reference node (empty dummy node)
        assert next_ref_node in self.get_edges(variant_node), "Failed parsing insertion %s. Found %d as next ref node after variant node, but there is no edge from variantnode %d to %d" % (variant, next_ref_node, variant_node, next_ref_node)
        # If there are multiple insertions, one could find the correct dummy node by choosing the one that goes to next node with lowest ref pos
        dummy_nodes = [node for node in self.get_edges(node) if next_ref_node in self.get_edges(node) and self.get_node_size(node) == 0]
        assert len(dummy_nodes) == 1, "Error when parsing insertion %s. There are not exactly 1 insertion node between %d and %d. Nodes of length 0 between are %s" % (variant, node, next_ref_node, dummy_nodes)
        insertion_node = dummy_nodes[0]
        assert next_ref_node in self.get_edges(insertion_node), "Failed parsing insertion %s. Found %d as next ref node after dummy node, but there is no edge from dumy node %d to %d" % (variant, next_ref_node, insertion_node, next_ref_node)

        return (insertion_node, variant_node)

    def get_variant_nodes(self, variant):
        if variant.type == "SNP":
            return self.get_snp_nodes(variant.position-1, variant.variant_sequence, variant.chromosome)
        elif variant.type == "DELETION":
            return self.get_deletion_nodes(variant.position-1, len(variant.ref_sequence)-1, variant.get_deleted_sequence(), variant.chromosome)
        elif variant.type == "INSERTION":
            return self.get_insertion_nodes(variant, variant.chromosome)

        raise Exception("Invalid variant %s. Has no type set." % variant)

    def set_allele_frequencies_from_variants(self, variants, use_chromosome=None):
        frequencies = np.zeros(len(self.nodes), dtype=np.float16) + 1  # Set all to 1.0 initially

        for i, variant in enumerate(variants):
            if use_chromosome is not None:
                variant.chromosome = use_chromosome

            if i % 1000000 == 0:
                logging.info("%d variants processed" % i)
            try:
                reference_node, variant_node = self.get_variant_nodes(variant)
                reference_af = variant.get_reference_allele_frequency()
            except VariantNotFoundException:
                logging.info("Variant %s not found in graph" % variant)
                continue
            variant_af = variant.get_variant_allele_frequency()
            frequencies[reference_node] = reference_af
            frequencies[variant_node] = variant_af

        self.allele_frequencies = frequencies

    def get_variant_nodes_in_region(self, chromosome, linear_ref_start, linear_ref_end):
        nodes = set(self.get_node_at_chromosome_and_chromosome_offset(chromosome, pos) for pos in range(linear_ref_start, linear_ref_end))
        return [self.get_edges(node) for node in nodes if len(self.get_edges(node)) == 2]

    def get_first_node(self):
        for candidate in [0, 1]:
            if len(self.get_edges(candidate)) > 0:
                return candidate

        raise Exception("Couldn't find first graph node")

    def get_haplotype_node_paths_for_haplotypes(self, variants, limit_to_n_haplotypes=10):
        # First find all variant nodes that the haplotype has
        haplotypes = list(range(0, limit_to_n_haplotypes))
        variant_nodes_in_haplotype = defaultdict(set)
        for i, variant in enumerate(variants):
            if i % 1000 == 0:
                logging.info("%d variants processed" % i)
            reference_node, variant_node = self.get_variant_nodes(variant)

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
        nodes = np.zeros((len(haplotypes), len(self.nodes)), dtype=np.uint32)
        for haplotype in haplotypes:
            current_node = self.get_first_node()
            i = 0
            while True:
                #nodes[haplotype].append(current_node)
                nodes[haplotype, i] = current_node

                next_nodes = self.get_edges(current_node)
                if len(next_nodes) == 0:
                    break

                next_node = None
                if len(next_nodes) == 1:
                    next_node = next_nodes[0]
                else:
                    for potential_next in next_nodes:
                        if potential_next in variant_nodes_in_haplotype[haplotype]:
                            next_node = potential_next

                assert next_node is not None

                current_node = next_node
                i += 1

            nodes[haplotype, i] = current_node

        return nodes

    def get_reverse_edges_dict(self):
        reverse_edges = defaultdict(list)
        for node in self.get_all_nodes():
            for edge in self.get_edges(node):
                reverse_edges[edge].append(node)

        return reverse_edges

    def convert_chromosome_ref_offset_to_graph_ref_offset(self, chromosome_offset, chromosome):
        # Add the graph ref offset at this chromosome to the chromosome offset
        chromosome_position = chromosome - 1
        return int(self.node_to_ref_offset[self.chromosome_start_nodes[chromosome_position]] + chromosome_offset)


    def find_nodes_from_node_that_matches_sequence(self, from_node, sequence, nodes_found, all_paths_found):
        #logging.info("== Searching from node %d with sequence %s. Variant type %s. Nodes found: %s. All paths found: %s" % (from_node, sequence, variant_type, nodes_found, all_paths_found))
        if sequence == "":
            # All of the sequence is used successfully, return
            all_paths_found.append(nodes_found)
            #logging.info("   RETURNING. All paths found: %s" % all_paths_found)
            return nodes_found, all_paths_found

        next_nodes = self.get_edges(from_node)
        result = (nodes_found, all_paths_found)
        for possible_next in next_nodes:

            #print("Checking next node %d" % possible_next)

            node_size = self.get_node_size(possible_next)
            if node_size == 0 or self.get_node_sequence(possible_next).lower() == sequence[0:node_size].lower():
                # This node is a match, we can continue searching
                new_sequence = sequence[node_size:]
                new_nodes_found = nodes_found.copy()
                new_nodes_found.append(possible_next)

                #print("   Matching sequence. New sequence is now %s" % new_sequence)
                result = self.find_nodes_from_node_that_matches_sequence(possible_next, new_sequence, new_nodes_found, all_paths_found)
                if not result:
                    continue
            else:
                #print("   No sequence match")
                continue

        return result

