import logging
logging.basicConfig(level=logging.INFO)
from obgraph import Graph
import numpy as np


def test_from_dicts():
    g = Graph.from_dicts(
        {1: "ACTG", 2: "A", 3: "G", 4: "AAA"},
        {1: [2, 3],
         2: [4],
         3: [4]},
        [1, 2, 4]
    )

    assert g.get_node_size(1) == 4
    assert g.get_node_size(2) == 1
    assert g.get_node_size(3) == 1
    assert g.get_node_size(4) == 3

    assert g.get_node_sequence(2) == "A"

    assert list(g.get_edges(1)) == [2, 3]


    assert list(g.get_numeric_node_sequence(2)) == [0]

    print(g.get_numeric_node_sequence(np.array([1, 2, 3, 4])))
    assert list(g.get_numeric_node_sequences(np.array([1, 2, 3, 4]))) == [0, 1, 2, 3, 0, 3, 0, 0, 0]

    assert g.get_node_at_ref_offset(0) == 1
    assert g.get_node_at_ref_offset(3) == 1
    assert g.get_node_at_ref_offset(4) == 2
    assert g.get_node_at_ref_offset(5) == 4

    assert g.get_ref_offset_at_node(4) == 5
    assert g.get_ref_offset_at_node(2) == 4


def test_sparse_graph():
    g = Graph.from_dicts(
        {1: "AGGG", 4: "CACCT"},
        {1: [4]},
        [1, 4]
    )
    assert list(g.get_edges(1)) == [4]
    assert g.get_node_at_ref_offset(4) == 4
    assert g.get_nodes_sequence([1, 4]) == "AGGGCACCT"


def test2():
    g = Graph.from_dicts(
        {0: "ACT", 1: "A", 2: "", 3: "GGG"},
        {0: [1, 2], 1: [3], 2: [3]},
        [0, 3]
    )

    assert set(g.get_edges(0)) == set([1, 2])

test_from_dicts()
test_sparse_graph()
test2()
