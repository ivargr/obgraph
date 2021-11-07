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

    def slice(self, from_variant, to_variant):
        return VariantToNodes(self.ref_nodes[from_variant:to_variant], self.var_nodes[from_variant:to_variant])

    @classmethod
    def from_graph_and_variants(cls, graph, variants):
        n_variants = len(variants)
        var_nodes = np.zeros(n_variants, dtype=np.uint32)
        ref_nodes = np.zeros(n_variants, dtype=np.uint32)

        max_graph_node = graph.max_node_id()
        for i, variant in enumerate(variants):
            if i % 100000 == 0:
                logging.info("%d variants processed" % i)
            try:
                ref_node, var_node = graph.get_variant_nodes(variant)
            except VariantNotFoundException:
                logging.erro("Could not find variant, aborting")
                raise

            var_nodes[i] = var_node
            ref_nodes[i] = ref_node

            assert var_node <= max_graph_node
            assert ref_node <= max_graph_node, "Found ref node %d which is not <= max graph node %d. Variant %s" % (ref_node, max_graph_node, variant)



        return cls(ref_nodes, var_nodes)


class NodeToVariants:
    properties = {"index"}
    def __init__(self, index):
        self.index = index

    @classmethod
    def from_file(cls, file_name):
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npz")

        return cls(data["index"])

    def to_file(self, file_name):
        np.savez(file_name, index=self.index)

    def get_variant_at_node(self, node):
        if self.index[node] == -1:
            return None
        return self.index[node]

    @classmethod
    def from_variant_to_nodes(cls, variant_to_nodes):
        n_nodes = max(np.max(variant_to_nodes.ref_nodes), np.max(variant_to_nodes.var_nodes))
        index = np.zeros(n_nodes+1, dtype=np.int32) - 1

        for variant in range(len(variant_to_nodes.ref_nodes+1)):
            if variant % 100000 == 0:
                logging.info("%d variants processed" % variant)

            ref_node = variant_to_nodes.ref_nodes[variant]
            var_node = variant_to_nodes.var_nodes[variant]

            index[ref_node] = variant
            index[var_node] = variant

        return cls(index)
