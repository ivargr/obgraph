from obgraph import Graph, DummyNodeAdder
from alignment_free_graph_genotyper.variants import GenotypeCalls, VariantGenotype

def test_simple_insertion():
    graph = Graph.from_dicts(
        {1: "ACTG", 2: "C", 3: "AAAA"},
        {1: [2, 3], 2: [3]},
        [1, 3]
    )

    variants = GenotypeCalls([VariantGenotype(4, "G", "GC", type="INSERTION")])
    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()

    assert new_graph.node_has_edges(5, [3])
    assert new_graph.node_has_edges(1, [2, 5])
    assert new_graph.node_has_edges(2, [3])
    #print(new_graph)


def test_double_deletion_with_snp_inside_first_deletion():

    graph = Graph.from_dicts(
        {1: "ACTG", 2: "A", 3: "C", 4: "T", 5: "AAA", 6: "G"},
        {
            1: [2, 5, 6],
            2: [3, 4],
            3: [5, 6],
            4: [5, 6],
            5: [6]
        },
        [1, 2, 4, 6]
    )

    variants = GenotypeCalls([VariantGenotype(1, 4, "GAT", "G", type="DELETION"), VariantGenotype(1, 6, "TAAA", "T", type="DELETION")])
    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()
    print(new_graph)

def test_double_deletion_with_snp_inside_first_deletion2():

    repeated_sequence = "AGGTCCCAGGTCCATCT"
    graph = Graph.from_dicts(
        {1: "TTTT", 2: "AGGTCC", 3: "C", 4: "A", 5: repeated_sequence, 6: repeated_sequence},
        {
            1: [2, 5, 6],
            2: [3, 4],
            3: [5, 6],
            4: [5, 6],
            5: [6]
        },
        [1, 2, 3, 5, 6]
    )

    variants = GenotypeCalls([VariantGenotype(4, "TAGGTCCC", "T", type="DELETION"), VariantGenotype(11, "CAGGTCCCAGGTCCATCT", "C", type="DELETION")])
    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()
    print(new_graph)

def test_insertion_with_multiple_paths():

    graph = Graph.from_dicts(
        {1: "AAAG", 2: "GAGT", 3: "GA", 4: "C", 5: "G", 6: "T"},
        {
            1: [2, 3],
            2: [3],
            3: [4, 5],
            4: [4, 6],
            5: [6]
        },
        [1, 3, 5, 6]
    )

    variants = GenotypeCalls([VariantGenotype(4, "G", "GGAGT", type="INSERTION")])
    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()
    print(new_graph)


def test_tricky_case_nested_deletions():
    graph = Graph.from_dicts(
        {1: "TATAT", 2: "AT", 3: "A", 4: "T", 5: "A", 6: "A", 7: "T", 8: "A", 9: "GG" },
        {
            1: [2, 6],
            2: [3, 6],
            3: [4, 5],
            4: [6],
            5: [6],
            6: [7, 8],
            7: [9],
            8: [9]
        },
        [1, 2, 3, 5, 6, 8, 9]
    )

    variants = GenotypeCalls(
        [
            VariantGenotype(1, 5, "TATAA", "T", type="DELETION"),
            VariantGenotype(1, 7, "TAA", "T", type="DELETION"),
            VariantGenotype(1, 5, "A", "T", type="SNP"),
        ]
    )

    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()
    print(new_graph)

    assert list(new_graph.get_edges(1)) == [2, 11]
    assert list(new_graph.get_edges(2)) == [3, 12]
    assert list(new_graph.get_edges(11)) == [6]
    assert list(new_graph.get_edges(12)) == [6]

def test_complex_case_multiple_equal_paths():
    pass

test_tricky_case_nested_deletions()
#test_simple_insertion()
#test_double_deletion_with_snp_inside_first_deletion()
#test_insertion_with_multiple_paths()