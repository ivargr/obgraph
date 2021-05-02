

class NeighbourHaplotypeNodes:
    def __init__(self, index):
        self._index = index

    @classmethod
    def from_variants(cls, variants, variant_to_nodes):
        for variant_index in range(len(variants)):
            # Iterate last 4 variants,


