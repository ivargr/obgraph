import numpy as np
from .graph import Graph
import logging


def add_indel_dummy_nodes(graph):
    node_ids, node_sequences, node_sizes, from_nodes, to_nodes, linear_ref_nodes, chromosome_start_nodes = graph.get_flat_graph()

    linear_ref_set = set(linear_ref_nodes)

    node_counter = node_ids[-1] + 1

    change_edges = {}
    reverse_edges = graph.get_reverse_edges_dict()
    # traverse graph, detect indels
    for node in node_ids:
        if node in linear_ref_set:
            edges = graph.get_edges(node)
            if len(edges) > 1:
                linear_ref_pos_end_of_node = graph.get_ref_offset_at_node(node) + graph.get_node_size(node)

                # Find all edges going to linear ref node where
                for next_node in edges:
                    # deletions and insertions

                    # New smarter detection that handles cases with overlapping variant
                    # if edge goes from linear ref node to another linear ref node where next linear ref node ref pos is not directly after
                    #   this must be a deletion, add a dummy node
                    # if edge goes directly to another linear ref node that follows on the linear ref and there are other edges that do not go to linear ref nodes and these go to the same linear ref node that is the next
                    #  there must be an insertion, this edge should be replace
                    #   this solution does not support insertions with more than one node in the insertion, but that does not happen??
                    next_node_linear_ref_pos = graph.get_ref_offset_at_node(next_node)
                    #logging.info("Checking edge from %d to %d. Linear ref end: %d. Next linear ref pos: %d" % (node, next_node, linear_ref_pos_end_of_node, next_node_linear_ref_pos))
                    """
                    if (next_node in linear_ref_set and next_node_linear_ref_pos > linear_ref_pos_end_of_node) \
                            or (next_node in linear_ref_set and next_node_linear_ref_pos == linear_ref_pos_end_of_node and \
                                len([other_next for other_next in edges if other_next != next_node and other_next not in linear_ref_set and next_node in graph.get_edges(other_next)]) > 0):
                        #logging.info("Adding edge between node %d and %d" % (node, next_node))
                        node_ids.append(node_counter)
                        node_sequences.append([""])
                        node_sizes.append(0)
                        from_nodes.append(node)
                        to_nodes.append(node_counter)
                        change_edges[(node, next_node)] = (node_counter, next_node)
                        node_counter += 1
                    """

                    # Alternative method that also deals with cases where an insertion node has multiple edges in
                    # insertions
                    if next_node in linear_ref_set and next_node_linear_ref_pos == linear_ref_pos_end_of_node:
                        #print("Checking insertion between %d and %d" % (node, next_node))
                        insertion_nodes = [n for n in edges if n != next_node and n not in linear_ref_set and next_node in graph.get_edges(n)]
                        #print("  Insertion nodes: %s" % insertion_nodes)
                        for insertion_node in insertion_nodes:
                            # There might be multiple insertions we assume
                            # For every node having an edge into the insertion node, we want to add a dummy node between that node and the next_node (next node on lienar ref)

                            nodes_in = reverse_edges[insertion_node]
                            #print("  Nodes going in to %d: %s" % (insertion_node, nodes_in))
                            if len(nodes_in) > 0:
                                node_ids.append(node_counter)
                                node_sequences.append([""])
                                node_sizes.append(0)
                                # Add edge from dummy node to next
                                from_nodes.append(node_counter)
                                to_nodes.append(next_node)
                                # Change all old edges from previous node to next_node to go to dummy node
                                for node_in in nodes_in:
                                    #print("   Changing node from %d-%d to %d-%d" % (node_in, next_node, node_in, node_counter))
                                    change_edges[(node_in, next_node)] = (node_in, node_counter)
                                node_counter += 1

                    # Deletion
                    if next_node in linear_ref_set and next_node_linear_ref_pos > linear_ref_pos_end_of_node:
                        #print("CHecking deletion from %d to %d" % (node, next_node))
                        # Insertion nodes can only be linear reference nodes that are before the next node
                        insertion_nodes = [n for n in edges if n != next_node and n in linear_ref_set and graph.get_ref_offset_at_node(n) < graph.get_ref_offset_at_node(next_node)]
                        #print("  Insertion nodes on linear ref: %s" % insertion_nodes)
                        for insertion_node in insertion_nodes:
                            # There cannot be multiple deletions to the same other ref node? Still looping
                            # For every node having an edge into the insertion node, we want to add a dummy node between that node and the next_node (next node on lienar ref)

                            nodes_in = reverse_edges[insertion_node]
                            if len(nodes_in) > 0:
                                node_ids.append(node_counter)
                                node_sequences.append([""])
                                node_sizes.append(0)
                                # Add edge from dummy node to next
                                from_nodes.append(node_counter)
                                to_nodes.append(next_node)
                                #print("   Adding new dummy node %d and edge from %d to %d" % (node_counter, node_counter, next_node))
                                # Change all old edges from previous node to next_node to go to dummy node
                                for node_in in nodes_in:
                                    #print("   Changing edte %d-%d to %d-%d" % (node_in, next_node, node_in, node_counter))
                                    change_edges[(node_in, next_node)] = (node_in, node_counter)
                                node_counter += 1


                    """
                    if next_node in linear_ref_set and len([other_next for other_next in edges if other_next != next_node and next_node in graph.get_edges(other_next)]) > 0:
                        node_ids.append(node_counter)
                        node_sequences.append([""])
                        node_sizes.append(0)
                        from_nodes.append(node)
                        to_nodes.append(node_counter)
                        change_edges[(node, next_node)] = (node_counter, next_node)
                        node_counter += 1
                    """

    # Change edges
    for i in range(len(from_nodes)):
        from_node = from_nodes[i]
        to_node = to_nodes[i]
        if (from_node, to_node) in change_edges:
            from_nodes[i] = change_edges[(from_node, to_node)][0]
            to_nodes[i] = change_edges[(from_node, to_node)][1]

    return Graph.from_flat_nodes_and_edges(np.array(node_ids), node_sequences, np.array(node_sizes), np.array(from_nodes), np.array(to_nodes), linear_ref_nodes, chromosome_start_nodes)
