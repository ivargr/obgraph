from obgraph.util import create_coordinate_map
from obgraph import Graph
import numpy as np


def test_create_coordinate_map():

    g = Graph.from_dicts(
        {1: "AAAA", 2: "G", 3: "GGGGGGG", 4: "ACA", 5: "C"},
        {1: [2, 3], 2: [4], 3: [4], 4: [5]},
        [1, 2, 4, 5]
    )

    path = np.array([1, 3, 4, 5])

    m = create_coordinate_map(path, g, chromosome_index=0)

    assert m[0] == 0
    assert m[11] == 5
    assert m[14] == 8

test_create_coordinate_map()

