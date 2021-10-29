from obgraph import Graph
from obgraph.variants import VcfVariants, VcfVariant

def test_find_insertion_nodes():
    g = Graph.from_dicts(
        {
            1: "CTACCA",
            2: "AA",
            3: "TAAATAA",
            4: ""
        },
        {
            1: [2, 4],
            2: [3],
            4: [3]

        },
        [1, 3]
    )
    variant = VcfVariant(1, 6, "A", "AAA", "", "INSERTION")


    ref_node, variant_node = g.get_variant_nodes(variant)
    assert ref_node == 4
    assert variant_node == 2

test_find_insertion_nodes()