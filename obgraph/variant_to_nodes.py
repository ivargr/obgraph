import numpy as np
from .graph import VariantNotFoundException
import logging


class VariantToNodes:
    properties = {"ref_nodes", "var_nodes"}
    def __init__(self, ref_nodes=None, var_nodes=None):
        self.ref_nodes = ref_nodes
        self.var_nodes = var_nodes

    @classmethod
    def from_file(cls, file_name):
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npz")

        return cls(data["ref_nodes"], data["var_nodes"])

    def to_file(self, file_name):
        np.savez(file_name, ref_nodes=self.ref_nodes, var_nodes=self.var_nodes)

    @classmethod
    def from_graph_and_variants(cls, graph, variants):
        n_variants = len(variants)
        var_nodes = np.zeros(n_variants, dtype=np.uint32)
        ref_nodes = np.zeros(n_variants, dtype=np.uint32)

        for i, variant in enumerate(variants):
            if i % 100000 == 0:
                logging.info("%d variants processed" % i)
            try:
                ref_node, var_node = graph.get_variant_nodes(variant)
            except VariantNotFoundException:
                continue

            var_nodes[i] = var_node
            ref_nodes[i] = ref_node

        return cls(ref_nodes, var_nodes)
