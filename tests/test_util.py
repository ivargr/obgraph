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


def test_insertion_with_snp_before():
    graph = Graph.from_dicts(
        {1: "A", 2: "A", 3: "A", 4: "A", 5: "A"},
        {1: [2, 3], 2: [4, 5], 3: [4, 5], 5: [4]},
        [1, 2, 4]
    )
    new_graph = add_indel_dummy_nodes(graph)
    assert list(new_graph.get_edges(2)) == [6, 5], new_graph.get_edges(2)
    assert list(new_graph.get_edges(3)) == [6, 5], new_graph.get_edges(3)
    assert list(new_graph.get_edges(6)) == [4], new_graph.get_edges(6)

def test_deletion_with_snp_before_and_snp_inside():
    graph = Graph.from_dicts(
        {1: "A", 2: "A", 3: "A", 4: "A", 5: "A", 6: "A", 6: "A", 7: "A"},
        {1: [2, 3], 2: [4, 7], 3: [4, 7], 5: [5, 6], 5: [7], 6: [7]},
        [1, 3, 4, 6, 7]
    )
    new_graph = add_indel_dummy_nodes(graph)
    assert list(new_graph.get_edges(3)) == [4, 8], new_graph.get_edges(3)
    assert list(new_graph.get_edges(2)) == [4, 8], new_graph.get_edges(3)
    assert list(new_graph.get_edges(8)) == [7], new_graph.get_edges(8)


def test_complex_deletions():
    graph = Graph.from_dicts(
        {1: "T", 2: "GTGA", 3: "GT", 4: "AAA"},
        {1: [2, 3, 4], 2: [3, 4], 3: [4]},
        [1, 2, 3, 4]
    )
    new_graph = add_indel_dummy_nodes(graph)
    assert 3 in list(new_graph.get_edges(2))
    #print(list(new_graph.get_edges(2)))

test_insertion_with_snp_before()
test_deletion_with_snp_before_and_snp_inside()
test_add_indel_dummy_nodes2()
test_complex_deletions()