import json
import logging
import re
import numpy as np

class Graph:
    def __init__(self, nodes, node_to_sequence_index, node_sequences, node_to_edge_index, node_to_n_edges, edges,
                 node_to_ref_offset, ref_offset_to_node, chromosome_start_nodes):
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

    def get_node_size(self, node):
        return self.nodes[node]

    def get_node_sequence(self, node):
        index_position = self.node_to_sequence_index[node]
        return ''.join(self.node_sequences[index_position:index_position+self.nodes[node]])

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

    def linear_ref_length(self):
        return len(self.ref_offset_to_node)

    def get_node_at_ref_offset(self, ref_offset):
        return self.ref_offset_to_node[ref_offset]

    def get_ref_offset_at_node(self, node):
        return self.node_to_ref_offset[node]

    def get_node_offset_at_ref_offset(self, ref_offset):
        node = self.get_node_at_ref_offset(ref_offset)
        offset_at_node = self.get_ref_offset_at_node(node)
        return ref_offset - offset_at_node

    def to_file(self, file_name):
        np.savez(file_name,
                 nodes=self.nodes,
                 node_to_sequence_index=self.node_to_sequence_index,
                 node_sequences=self.node_sequences,
                 node_to_edge_index=self.node_to_edge_index,
                 node_to_n_edges=self.node_to_n_edges,
                 edges=self.edges,
                 node_to_ref_offset=self.node_to_ref_offset,
                 ref_offset_to_node=self.ref_offset_to_node,
                 chromosome_start_nodes=self.chromosome_start_nodes
                 )

        logging.info("Saved to file %s" % file_name)

    @classmethod
    def from_file(cls, file_name):
        logging.info("Reading file from %s" % file_name)
        data = np.load(file_name)
        return cls(data["nodes"],
                   data["node_to_sequence_index"],
                   data["node_sequences"],
                   data["node_to_edge_index"],
                   data["node_to_n_edges"],
                   data["edges"],
                   data["node_to_ref_offset"],
                   data["ref_offset_to_node"],
                   data["chromosome_start_nodes"])


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
        nodes = list(node_sequences.keys())
        node_sequences = [list(node_sequences[node]) for node in nodes]
        node_sizes = [len(seq) for seq in node_sequences]
        from_nodes = []
        to_nodes = []
        for from_node, to_nodes_edges in edges.items():
            for to_node in to_nodes_edges:
                from_nodes.append(from_node)
                to_nodes.append(to_node)

        return Graph.from_flat_nodes_and_edges(nodes, node_sequences, node_sizes, np.array(from_nodes), np.array(to_nodes), np.array(linear_ref_nodes), [node for node in nodes if node not in to_nodes])

    @classmethod
    def from_flat_nodes_and_edges(cls, node_ids, node_sequences, node_sizes, from_nodes, to_nodes, linear_ref_nodes, chromosome_start_nodes):

        logging.info("Asserting linear ref nodes are not empty")
        for i, node in enumerate(node_ids):
            size = node_sizes[i]
            if size == 0 and node in linear_ref_nodes:
                raise Exception("Placeholder nodes (for insertions) cannot be marked as linear ref nodes, since that breaks ref position to node lookup. Best solution is to not put these in linear ref nodes")


        max_node = np.max(node_ids)
        nodes = np.zeros(max_node+1, dtype=np.uint8)
        nodes[node_ids] = node_sizes

        # Node sequences
        new_node_positions = np.cumsum(node_sizes)
        node_sequence_index = np.zeros(max_node+1, dtype=np.uint64)
        node_sequence_index[node_ids[0]] = 0
        node_sequence_index[node_ids[1:]] = new_node_positions[:-1]
        node_sequences_array = np.concatenate(node_sequences).astype("<U1")


        #node_sequences_array = np.zeros(max_node+0, dtype=object)
        #logging.info("Allocating empty sequence array for node sequences")
        #node_sequences_array = np.array([np.zeros(size, dtype="<U1") for size in nodes])  # Must allocate space
        #node_sequences_array[node_ids] = node_sequences

        #logging.info("Sorting nodes")
        sorting = np.argsort(from_nodes)
        from_nodes = from_nodes[sorting]
        to_nodes = to_nodes[sorting]

        diffs = np.ediff1d(from_nodes, to_begin=1)
        positions_of_unique_nodes = np.nonzero(diffs)[0]
        unique_nodes = from_nodes[positions_of_unique_nodes]

        node_to_edge_index = np.zeros(max_node+1, dtype=np.uint32)
        node_to_n_edges = np.zeros(max_node+1, dtype=np.uint8)
        node_to_edge_index[unique_nodes] = positions_of_unique_nodes
        n_edges_numbers = np.ediff1d(positions_of_unique_nodes, to_end=len(from_nodes)-positions_of_unique_nodes[-1])
        node_to_n_edges[unique_nodes] = n_edges_numbers

        #logging.info("Finding ref offsets")
        node_to_ref_offset = np.zeros(max_node+1, np.uint64)
        ref_node_sizes = nodes[linear_ref_nodes]
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


        node_ids = np.array(node_ids)
        node_sizes = np.array(node_sizes)
        edges_from = np.array(edges_from)
        edges_to = np.array(edges_to)
        linear_ref_nodes = np.array(linear_ref_nodes, dtype=np.uint32)
        chromosome_start_nodes = np.array(chromosome_start_nodes, dtype=np.uint32)

        return cls.from_flat_nodes_and_edges(node_ids, node_sequences, node_sizes, edges_from, edges_to, linear_ref_nodes, chromosome_start_nodes)


    def get_snp_nodes(self, ref_offset, variant_bases, chromosome=1):
        node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        node_offset = self.get_node_offset_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        assert node_offset == 0
        prev_node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset - 1)

        # Try to find next node that matches read base
        for potential_next in self.get_edges(prev_node):
            if potential_next == node:
                continue

            node_seq = self.get_node_sequence(potential_next)[0]
            if node_seq.lower() == variant_bases.lower():
                return node, potential_next

        logging.error("Could not parse substitution at offset %d with bases %s" % (ref_offset, variant_bases))
        raise Exception("Parseerrror")

    def get_deletion_nodes(self, ref_offset, deletion_length, chromosome=1):
        ref_offset += 1
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
            logging.error("Offset %d is not at beginning of node" % next_ref_pos)
            logging.error("Node at %d: %d" % (next_ref_pos, next_ref_node))
            logging.error("Ref length in deletion: %s" % deletion_length)
            logging.info("Ref pos beginning of deletion: %d" % ref_offset)
            raise Exception("Deletion not in graph")

        # Find an empty node between prev_node and next_ref_node
        deletion_nodes = [node for node in self.get_edges(prev_node) if next_ref_node in self.get_edges(node) and self.get_node_size(node) == 0]
        #debug = [(node, self.get_edges(node), self.get_node_size(node)) for node in self.get_edges(prev_node)]
        #logging.info("%s" % debug)
        assert len(deletion_nodes) == 1, "There should be only one deletion node between %d and %d. There are %d. Edges out from %d are: %s. Edges out from %d are %s" % (prev_node, next_ref_node, len(deletion_nodes), prev_node, self.get_edges(prev_node), node, self.get_edges(node))

        deletion_node = deletion_nodes[0]
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
        variant_node = None
        for potential_next in self.get_edges(node):
            if potential_next in self.linear_ref_nodes():
                continue  # Next should not be in linear ref

            #print("Processing insertion at ref offset %d with base %s" % (ref_offset, base))
            #print("  Node %d with offset %d. Node size: %d" % (node, node_offset, self.blocks[node].length()))

            next_bases = self.get_node_sequence(potential_next)[0:insertion_length].upper()
            if next_bases == insertion_sequence.upper():
                variant_node = potential_next
                break

        assert variant_node is not None, "Could not find insertion node. No nodes going out from %d has sequence %s" % (node, variant.variant_sequence)

        # Find linear reference node
        next_ref_pos = ref_offset + insertion_length
        next_ref_node = self.get_node_at_chromosome_and_chromosome_offset(chromosome, ref_offset)
        assert next_ref_node in self.get_edges(variant_node), "Failed parsing insertion %s. Found %d as next ref node after variant node, but there is no edge from variantnode %d to %d" % (variant, next_ref_node, variant_node, next_ref_node)
        # If there are multiple insertions, one could find the correct dummy node by choosing the one that goes to next node with lowest ref pos
        dummy_nodes = [node for node in self.get_edges(node) if next_ref_node in self.get_edges(node) and self.get_node_size(node) == 0]
        assert len(dummy_nodes) == 1, "There are not exactly 1 insertion node between %d and %d. Nodes of length 0 between are %s" % (node, next_ref_node, dummy_nodes)
        insertion_node = dummy_nodes[0]
        assert next_ref_node in self.get_edges(insertion_node), "Failed parsing insertion %s. Found %d as next ref node after dummy node, but there is no edge from dumy node %d to %d" % (variant, next_ref_node, insertion_node, next_ref_node)

        return (insertion_node, variant_node)

    def get_variant_nodes(self, variant, chromosome=1):
        if variant.type == "SNP":
            return self.get_snp_nodes(variant.position-1, variant.variant_sequence, chromosome)
        elif variant.type == "DELETION":
            return self.get_deletion_nodes(variant.position-1, len(variant.ref_sequence)-1, chromosome)
        elif variant.type == "INSERTION":
            return self.get_insertion_nodes(variant, chromosome)

        raise Exception("Invalid variant %s. Has no type set." % variant)

