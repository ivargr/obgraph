from obgraph.util import add_indel_dummy_nodes
from obgraph import Graph
import numpy as np


def test_add_indel_dummy_nodes():
    graph = Graph.from_flat_nodes_and_edges(
        np.array([1, 2, 3]),
        [list("ACT"), list("ACT"), list("GGG")],
        np.array([3, 3, 3]),
        np.array([1, 1, 2]),
        np.array([2, 3, 3]),
        np.array([1, 3]),
        None
    )


    new_graph = add_indel_dummy_nodes(graph)

    assert list(new_graph.get_edges(1)) == [2, 4]
    assert list(new_graph.get_edges(4)) == [3]

    assert new_graph.get_node_size(4) == 0
    assert new_graph.get_node_sequence(4) == ""

def test_add_indel_dummy_nodes2():
    graph = Graph.from_flat_nodes_and_edges(
        np.array([1, 2, 3, 4, 5, 6, 7]),
        [list("ACT"), list("ACT"), list("GGG"), list("A"), list("A"),  list("A"),  list("A") ],
        np.array([3, 3, 3, 1, 1, 1, 1]),
        np.array([1, 1, 2, 3, 4, 4, 5, 6]),
        np.array([2, 3, 3, 4, 5, 6, 7, 7]),
        np.array([1, 3, 4, 5, 7]),
        None
    )


    new_graph = add_indel_dummy_nodes(graph)

    assert list(new_graph.get_edges(1)) == [2, 8], "Has edges %s" % new_graph.get_edges(1)
    assert list(new_graph.get_edges(8)) == [3]

    assert new_graph.get_node_size(8) == 0
    assert new_graph.get_node_sequence(8) == ""

    assert list(new_graph.get_edges(4)) == [5, 6], new_graph.get_edges(4)

test_add_indel_dummy_nodes2()