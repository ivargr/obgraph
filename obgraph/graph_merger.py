import logging
import numpy as np
from .graph import Graph


def merge_graphs(graphs):


    # Makes list of all data, then concatenates those in the end into numpy arrays that can be used to create a new graph
    new_nodes = []
    new_nodes_to_sequence_index = []
    new_node_sequences = []
    new_node_to_edge_index = []
    new_node_to_n_edges = []
    new_edges = []
    new_node_to_ref_offset = []
    new_ref_offset_to_node = []
    new_chromosome_start_nodes = []
    new_allele_frequencies = []

    node_offset = 0
    node_sequence_offset = 0
    edge_index_offset = 0
    ref_offset = 0

    for i, graph in enumerate(graphs):
        logging.info("Processing graph %i" % i)
        logging.info("Node offset is now %d" % node_offset)

        # Add nodes, these are simply just appended
        new_nodes.append(graph.nodes)

        # All node ids should be increased by node_offset

        # Node to sequence index can just be appended and index positions should increase by the length
        # of the index up to now (new indexes needs to point this offset further in the sequence array)
        new_nodes_to_sequence_index.append(graph.node_to_sequence_index+node_sequence_offset)

        new_node_sequences.append(graph.node_sequences)
        new_node_to_edge_index.append(graph.node_to_edge_index+edge_index_offset)
        new_edges.append(graph.edges+node_offset)
        new_node_to_n_edges.append(graph.node_to_n_edges)

        new_node_to_ref_offset.append(graph.node_to_ref_offset+ref_offset)
        new_ref_offset_to_node.append(graph.ref_offset_to_node+node_offset)

        new_chromosome_start_nodes.append(graph.get_first_node() + node_offset)

        if graph.allele_frequencies is not None:
            new_allele_frequencies.append(graph.allele_frequencies)

        # Increase offsets
        node_offset += len(graph.nodes)
        node_sequence_offset += len(graph.node_sequences)
        edge_index_offset += len(graph.edges)

        # increase the length of the linear reference genome
        ref_offset += len(graph.ref_offset_to_node)

        logging.info("Ref offset is now %d" % ref_offset)

    logging.info("Concatenating all data")
    new_nodes = np.concatenate(new_nodes)
    new_nodes_to_sequence_index = np.concatenate(new_nodes_to_sequence_index)
    new_node_sequences = np.concatenate(new_node_sequences)
    new_node_to_edge_index = np.concatenate(new_node_to_edge_index)
    new_node_to_n_edges = np.concatenate(new_node_to_n_edges)
    new_edges = np.concatenate(new_edges)
    new_node_to_ref_offset = np.concatenate(new_node_to_ref_offset)
    new_ref_offset_to_node = np.concatenate(new_ref_offset_to_node)
    new_chromosome_start_nodes = np.array(new_chromosome_start_nodes)
    if len(new_allele_frequencies) > 0:
        new_allele_frequencies = np.concatenate(new_allele_frequencies)
    else:
        new_allele_frequencies = None


    return Graph(new_nodes, new_nodes_to_sequence_index, new_node_sequences,
                 new_node_to_edge_index, new_node_to_n_edges, new_edges,
                 new_node_to_ref_offset, new_ref_offset_to_node, new_chromosome_start_nodes,
                 new_allele_frequencies)
