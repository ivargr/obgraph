import logging
from .mutable_graph import MutableGraph
from .graph import Graph
import numpy as np


def create_graph_from_gfa_file(file_name):

    graph = MutableGraph()

    path_lines = []

    with open(file_name) as f:
        for line in f:
            l = line.split()
            if line.startswith("S"):
                id = int(l[1])
                sequence = l[2]
                graph.add_node(id, sequence)
            elif line.startswith("L"):
                from_node = int(l[1])
                to_node = int(l[3])

                assert l[2] == "+" and l[4] == "+", "Only links from positive side to positive side are supported"
                assert l[5] == "*", "Overlaps between segments are not supported"

                graph.add_edge(from_node, to_node)
            elif line.startswith("P"):
                if not l[1].startswith("_alt"):
                    path_lines.append(l)

    #assert len(path_lines) == 1, "Supports only one path, which should be reference path. Found %d" % len(path_lines)

    # find reference node from path
    linear_ref_nodes = []
    chromosome_start_nodes = []
    for path_line in path_lines:
        path_name = path_line[1]
        logging.info("Processing path %s" % path_name)
        nodes = path_line[2].split(",")
        nodes = [int(node.replace("+", "")) for node in nodes]
        chromosome_start_nodes.append(nodes[0])
        linear_ref_nodes.extend(nodes)
        logging.info("There are %d linear ref nodes in path %s"  % (len(nodes), path_name))

    graph.set_linear_ref_nodes(linear_ref_nodes)
    graph.chromosome_start_nodes = np.array(chromosome_start_nodes)

    logging.info("Chromosome start nodes: %s" % graph.chromosome_start_nodes)
    return Graph.from_mutable_graph(graph)


