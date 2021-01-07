from obgraph.mutable_graph import MutableGraph

def test_create():
    graph = MutableGraph({1: 4, 2: 3, 3: 1, 4: 1}, {1: "ACTG", 2: "A", 3: "C", 4: "AAAA"}, {1: [2, 3], 3: [4], 2: [4]})

    assert graph.get_node_size(1) == 4
    assert graph.get_edges(1) == [2, 3]
    assert graph.get_nodes_before(2) == [1]
    assert 2 in graph.get_nodes_before(4)
    assert 3 in graph.get_nodes_before(4)
    assert graph.get_node_sequence(2) == "A"


def test_get_nodes_matching_sequence_single_node():
    graph = MutableGraph({1: 4, 2: 3, 3: 1, 4: 1}, {1: "ACTG", 2: "A", 3: "C", 4: "AAAA"}, {1: [2, 3], 3: [4], 2: [4]})

    nodes = graph.find_nodes_from_node_that_matches_sequence(1, "A")
    assert nodes == [2]

def test_get_nodes_matching_sequence_double_deletion_and_snp():
    graph = MutableGraph(
        {1: 4, 2: 1, 3: 1, 4: 1, 5: 3, 6: 1},
        {1: "ACTG", 2: "A", 3: "C", 4: "T", 5: "AAA", 6: "G"},
        {
            1: [2, 5, 6],
            2: [3, 4],
            3: [5, 6],
            4: [5, 6],
            5: [6]
        }
    )

    nodes = graph.find_nodes_from_node_that_matches_sequence(1, "AT")
    assert nodes == [2, 4]

    nodes = graph.find_nodes_from_node_that_matches_sequence(1, "AAA")
    assert nodes == [5]


test_create()
test_get_nodes_matching_sequence_single_node()
test_get_nodes_matching_sequence_double_deletion_and_snp()
