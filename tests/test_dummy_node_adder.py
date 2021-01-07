from obgraph import Graph, DummyNodeAdder
from alignment_free_graph_genotyper.variants import GenotypeCalls, VariantGenotype

def test_simple_insertion():
    graph = Graph.from_dicts(
        {1: "ACTG", 2: "C", 3: "AAAA"},
        {1: [2, 3], 2: [3]},
        [1, 3]
    )

    variants = GenotypeCalls([VariantGenotype(3, "G", "GC", type="INSERTION")])
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

    variants = GenotypeCalls([VariantGenotype(3, "GAT", "G", type="DELETION"), VariantGenotype(5, "TAAA", "T", type="DELETION")])
    dummy_node_adder = DummyNodeAdder(graph, variants)
    new_graph = dummy_node_adder.create_new_graph_with_dummy_nodes()
    print(new_graph)

#test_simple_insertion()
test_double_deletion_with_snp_inside_first_deletion()