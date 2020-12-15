import numpy as np
from .graph import Graph


def add_indel_dummy_nodes(graph):
    node_ids, node_sequences, node_sizes, from_nodes, to_nodes, linear_ref_nodes, chromosome_start_nodes = graph.get_flat_graph()

    linear_ref_set = set(linear_ref_nodes)

    node_counter = node_ids[-1] + 1

    change_edges = {}

    # traverse graph, detect indels
    for node in node_ids:
        if node in linear_ref_set:
            edges = graph.get_edges(node)
            if len(edges) > 1:

                # Find all edges going to linear ref node where
                for next_node in edges:
                    if next_node in linear_ref_set and len([other_next for other_next in edges if other_next != next_node and next_node in graph.get_edges(other_next)]) > 0:
                        node_ids.append(node_counter)
                        node_sequences.append([""])
                        node_sizes.append(0)
                        from_nodes.append(node)
                        to_nodes.append(node_counter)
                        change_edges[(node, next_node)] = (node_counter, next_node)
                        node_counter += 1



    # Change edges
    for i in range(len(from_nodes)):
        from_node = from_nodes[i]
        to_node = to_nodes[i]
        if (from_node, to_node) in change_edges:
            from_nodes[i] = change_edges[(from_node, to_node)][0]
            to_nodes[i] = change_edges[(from_node, to_node)][1]

    return Graph.from_flat_nodes_and_edges(np.array(node_ids), node_sequences, np.array(node_sizes), np.array(from_nodes), np.array(to_nodes), linear_ref_nodes, chromosome_start_nodes)
