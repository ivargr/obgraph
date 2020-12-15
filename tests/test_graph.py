from obgraph import Graph


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

    assert list(g.get_edges(1)) == [2, 3]

    assert g.get_node_sequence(2) == "A"


test_from_dicts()
