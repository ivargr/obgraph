import logging


def traverse_graph_by_following_nodes(graph, follow_nodes, only_add_nodes_in_follow_nodes=False):
    assert type(follow_nodes) == set

    logging.info("N nodes in follow set: %d" % len(follow_nodes))

    nodes_found = []
    linear_ref_nodes = graph.linear_ref_nodes()
    current_node = graph.get_first_node()
    i = 0
    while True:
        if current_node in follow_nodes and only_add_nodes_in_follow_nodes:
            nodes_found.append(current_node)
        else:
            nodes_found.append(current_node)

        next_nodes = graph.get_edges(current_node)
        if len(next_nodes) == 0:
            break

        if i % 100000 == 0:
            logging.info("%i nodes traversed (total is approx %d)" % (i, len(linear_ref_nodes)))

        next_node = None
        if len(next_nodes) == 1:
            next_node = next_nodes[0]
        else:
            for potential_next in next_nodes:
                if potential_next in follow_nodes or (next_node is None and potential_next in linear_ref_nodes):
                    next_node = potential_next

            # IF did not find any nodes, choose an empty node if there is one
            if next_node is None:
                next_node = [n for n in next_nodes if graph.get_node_size(n) == 0]
                assert len(next_node) == 1, "There are multiple empty nodes from node %d: %s" % (current_node, next_nodes)
                next_node = next_node[0]


        assert next_node is not None, "Could not find any new nodes from %d. Edges are %s" % (current_node, str(next_nodes))
        current_node = next_node
        i += 1

    if current_node in follow_nodes or not only_add_nodes_in_follow_nodes:
        nodes_found.append(current_node)

    return nodes_found
