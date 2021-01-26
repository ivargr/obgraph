from obgraph.graph_construction import GraphConstructor as GraphConstructor
from alignment_free_graph_genotyper.variants import GenotypeCalls, VariantGenotype


def test_single_snp():
    reference = "ACTGGG"

    variants = GenotypeCalls(
        [VariantGenotype(1, 2, "C", "T", type="SNP")]
    )

    constructor = GraphConstructor(reference, variants)
    graph = constructor.get_graph()

    assert list(graph.get_edges(1)) == [2, 3], graph.get_edges(1)
    assert list(graph.get_edges(2)) == [4], graph.get_edges(2)
    assert list(graph.get_edges(3)) == [4], graph.get_edges(3)
    assert list(graph.get_edges(4)) == [] , graph.get_edges(4)

    assert graph.get_node_size(1) == 1
    assert graph.get_node_size(2) == 1
    assert graph.get_node_size(3) == 1
    assert graph.get_node_size(4) == 4

    assert graph.get_node_sequence(1) == "A"
    assert graph.get_node_sequence(2) == "T"
    assert graph.get_node_sequence(3) == "C"
    assert graph.get_node_sequence(4) == "TGGG"

    assert graph.linear_ref_nodes() == set([1, 3, 4])


def test_single_deletion():
    reference = "AATTGG"

    variants = GenotypeCalls(
        [VariantGenotype(1, 2, "ATT", "A", type="DELETION")]
    )

    constructor = GraphConstructor(reference, variants)
    graph = constructor.get_graph()
    print(graph)

    """
    assert list(graph.get_edges(1)) == [2, 3]
    assert list(graph.get_edges(2)) == [4]
    assert list(graph.get_edges(3)) == [4]

    assert graph.get_node_sequence(2) == ""
    assert graph.get_node_sequence(3) == "TT"
    assert graph.get_node_sequence(1) == "AA"
    assert graph.get_node_sequence(4) == "GG"

    assert graph.linear_ref_nodes() == set([1, 3, 4])
    """


def test_single_insertion():
    reference = "AATTGG"

    variants = GenotypeCalls(
        [VariantGenotype(1, 2, "A", "AAA", type="INSERTION")]
    )

    constructor = GraphConstructor(reference, variants)
    graph = constructor.get_graph()

    print(graph)

def test_double_deletion_with_snp_inside_first_deletion():

    reference = "ACTGATAAAG"
    variants = GenotypeCalls([VariantGenotype(1, 4, "GAT", "G", type="DELETION"), VariantGenotype(1, 6, "T", "C", type="SNP"), VariantGenotype(1, 6, "TAAA", "T", type="DELETION")])
    constructor = GraphConstructor(reference, variants)
    graph = constructor.get_graph()
    print(graph.get_node_at_ref_offset(0))
    print(graph.get_node_at_ref_offset(1))
    graph_with_dummy_nodes = constructor.get_graph_with_dummy_nodes()
    print(graph_with_dummy_nodes)

def test_insertion_with_snp_right_before():
    reference = "AAAAAA"
    variants = GenotypeCalls([VariantGenotype(1, 2, "A", "T", type="SNP"), VariantGenotype(1, 2, "A", "AC", type="INSERTION")])
    constructor = GraphConstructor(reference, variants)
    graph = constructor.get_graph()

    assert list(graph.get_edges(1)) == [2, 3]
    assert list(graph.get_edges(2)) == [4, 5]
    assert list(graph.get_edges(3)) == [4, 5]
    assert list(graph.get_edges(4)) == [6]
    assert list(graph.get_edges(5)) == [6]
    assert graph.get_node_sequence(5) == ""
    assert graph.get_node_sequence(4) == "C"
    assert graph.linear_ref_nodes() == set([1, 3, 6])


def test_insertion_with_snp_right_before_and_right_after():
    reference = "AAAAAA"
    variants = GenotypeCalls([VariantGenotype(1, 2, "A", "T", type="SNP"), VariantGenotype(1, 2, "A", "AC", type="INSERTION"), VariantGenotype(1, 3, "A", "C", type="SNP")])
    constructor = GraphConstructor(reference, variants)
    graph = constructor.get_graph()
    print(graph)

def test_deletion_with_snp_right_before_and_right_after():
    reference = "AAAAAA"
    variants = GenotypeCalls([VariantGenotype(1, 2, "A", "T", type="SNP"), VariantGenotype(1, 2, "AA", "A", type="DELETION"), VariantGenotype(1, 3, "A", "C", type="SNP")])
    constructor = GraphConstructor(reference, variants)
    graph = constructor.get_graph()
    graph_with_dummy_nodes = constructor.get_graph_with_dummy_nodes()
    print(graph)
    print(graph_with_dummy_nodes)


def test_messy_graph():
    reference = "GCATATTTT"
    variants = GenotypeCalls(
        [
            VariantGenotype(1, 2, "CAT", "C", type="DELETION"),
            VariantGenotype(1, 3, "A", "G", type="SNP"),
            VariantGenotype(1, 4, "TA", "T", type="DELETION"),
            VariantGenotype(1, 5, "A", "AT", type="INSERTION"),
        ]
    )

    constructor = GraphConstructor(reference, variants)
    graph = constructor.get_graph()
    print(graph)
    print(constructor.get_graph_with_dummy_nodes())

    ref, var = graph.get_variant_nodes(variants[2])
    assert ref == 5
    assert var == 11


#test_single_snp()
#test_single_deletion()
#test_single_insertion()
#test_double_deletion_with_snp_inside_first_deletion()
#test_insertion_with_snp_right_before_and_right_after()
#test_deletion_with_snp_right_before_and_right_after()
test_messy_graph()

