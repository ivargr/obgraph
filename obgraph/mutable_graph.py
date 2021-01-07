from collections import defaultdict


class MutableGraph:
    def __init__(self, nodes, node_sequences, edges, linear_ref_nodes=None, node_to_ref_offset=None, ref_offset_to_node=None, chromosome_start_nodes=None, allele_frequencies=None):
        self.nodes = nodes
        self.node_sequences = node_sequences
        self.edges = edges
        self.linear_ref_nodes = linear_ref_nodes
        print("Linear ref nodes: %s" % self.linear_ref_nodes)
        self.node_to_ref_offset = node_to_ref_offset
        self.ref_offset_to_node = ref_offset_to_node
        self.chromosome_start_nodes = chromosome_start_nodes
        self.allele_frequencies = allele_frequencies
        self.reverse_edges = defaultdict(list)
        self.get_reverse_edges()

    def __str__(self):
        description = "Nodes: %s " % (self.node_sequences)
        description += "\nEdges: %s" % (self.edges)
        return description

    def get_reverse_edges(self):
        for node in self.get_all_nodes():
            for edge in self.get_edges(node):
                self.reverse_edges[edge].append(node)

    def get_all_nodes(self):
        return list(self.nodes.keys())

    def get_node_size(self, node):
        return self.nodes[node]

    def get_node_sequence(self, node):
        return self.node_sequences[node]

    def get_edges(self,node):
        if node in self.edges:
            return self.edges[node]
        else:
            return []

    def remove_edge(self, from_node, to_node):
        self.edges[from_node].remove(to_node)
        self.reverse_edges[to_node].remove(from_node)

    def add_edge(self, from_node, to_node):
        if from_node not in self.edges:
            self.edges[from_node] = []

        self.edges[from_node].append(to_node)
        print("   Adding edge from %d to %d" % (from_node, to_node))

        if to_node not in self.reverse_edges:
            self.reverse_edges[to_node] = []

        self.reverse_edges[to_node].append(from_node)

    def add_node(self, id, sequence=""):
        self.nodes[id] = len(sequence)
        self.node_sequences[id] = sequence

    def get_nodes_before(self, node):
        return self.reverse_edges[node]

    def find_nodes_from_node_that_matches_sequence(self, from_node, sequence, nodes_found=[]):

        if sequence == "":
            # All of the sequence is used successfully, return
            return nodes_found

        next_nodes = self.get_edges(from_node)

        for possible_next in next_nodes:
            #print("Checking next node %d" % possible_next)
            node_size = self.get_node_size(possible_next)
            if node_size == 0 or self.get_node_sequence(possible_next).lower() == sequence[0:node_size].lower():
                # This node is a match, we can continue searching
                new_sequence = sequence[node_size:]
                new_nodes_found = nodes_found.copy()
                new_nodes_found.append(possible_next)
                #print("   Matching sequence. New sequence is now %s" % new_sequence)
                result = self.find_nodes_from_node_that_matches_sequence(possible_next, new_sequence, new_nodes_found)
                if not result:
                    continue
                else:
                    return result
            else:
                #print("   No sequence match")
                continue

        # No match in this search path
        return False

