import json
import logging
import re
import numpy as np

class Graph:
    def __init__(self, nodes, node_to_sequence_index, node_sequences, node_to_edge_index, node_to_n_edges, edges, node_to_ref_offset, ref_offset_to_node):
        self.nodes = nodes
        self.node_to_sequence_index = node_to_sequence_index
        self.node_sequences = node_sequences
        self.node_to_edge_index = node_to_edge_index
        self.node_to_n_edges = node_to_n_edges
        self.edges = edges
        self.node_to_ref_offset = node_to_ref_offset
        self.ref_offset_to_node = ref_offset_to_node
        self._linear_ref_nodes_cache = None

    def get_node_size(self, node):
        return self.nodes[node]

    def get_node_sequence(self, node):
        index_position = self.node_to_sequence_index[node]
        return ''.join(self.node_sequences[index_position:index_position+self.nodes[node]])

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
                 ref_offset_to_node=self.ref_offset_to_node)

        logging.info("Saved to file %s" % file_name)

    @classmethod
    def from_file(cls, file_name):
        data = np.load(file_name)
        return cls(data["nodes"],
                   data["node_to_sequence_index"],
                   data["node_sequences"],
                   data["node_to_edge_index"],
                   data["node_to_n_edges"],
                   data["edges"],
                   data["node_to_ref_offset"],
                   data["ref_offset_to_node"])

    @classmethod
    def from_flat_nodes_and_edges(cls, node_ids, node_sequences, node_sizes, from_nodes, to_nodes, linear_ref_nodes):
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

        logging.info("Sorting nodes")
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

        logging.info("Finding ref offsets")
        node_to_ref_offset = np.zeros(max_node+1, np.uint64)
        ref_node_sizes = nodes[linear_ref_nodes]
        ref_offsets = np.cumsum(ref_node_sizes) - ref_node_sizes[0]
        node_to_ref_offset[linear_ref_nodes] = ref_offsets

        ref_offset_to_node = np.zeros(int(np.max(node_to_ref_offset)) + 32)
        index_positions = np.cumsum(ref_node_sizes)[:-1]
        ref_offset_to_node[index_positions] = np.diff(linear_ref_nodes)
        ref_offset_to_node[0] = linear_ref_nodes[0]
        ref_offset_to_node = np.cumsum(ref_offset_to_node, dtype=np.uint32)

        return cls(nodes, node_sequence_index, node_sequences_array, node_to_edge_index,
                   node_to_n_edges, to_nodes, node_to_ref_offset, ref_offset_to_node)

    @classmethod
    def from_vg_json_files(cls, file_names):

        node_sequences = []
        node_ids = []
        node_sizes = []
        edges_from = []
        edges_to = []
        linear_ref_nodes = []
        n_nodes_added = 0
        n_edges_added = 0
        for file_name in file_names:
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
                        path_found = True


        node_ids = np.array(node_ids)
        node_sizes = np.array(node_sizes)
        edges_from = np.array(edges_from)
        edges_to = np.array(edges_to)
        linear_ref_nodes = np.array(linear_ref_nodes, dtype=np.uint32)

        return cls.from_flat_nodes_and_edges(node_ids, node_sequences, node_sizes, edges_from, edges_to, linear_ref_nodes)



