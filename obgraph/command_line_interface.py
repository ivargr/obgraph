import logging
logging.basicConfig(level=logging.INFO, format='%(module)s %(asctime)s %(levelname)s: %(message)s')
import sys
import argparse
from . import Graph
from .util import add_indel_dummy_nodes
from alignment_free_graph_genotyper.variants import GenotypeCalls
from .haplotype_nodes import HaplotypeToNodes, NodeToHaplotypes
from .dummy_node_adder import DummyNodeAdder
from .haplotype_nodes import NodeToHaplotypes
from .genotype_matrix import GenotypeMatrix, GenotypeMatrixAnalyser, GenotypeFrequencies

from . import cython_traversing

#from .cython_traversing import traverse_graph_by_following_nodes

import time

def make(args):
    logging.info("Will create from files %s" % args.vg_json_files)
    graph = Graph.from_vg_json_files(args.vg_json_files)
    graph.to_file(args.out_file_name)

def add_indel_nodes(args):
    graph = Graph.from_file(args.graph_file_name)
    new_graph = add_indel_dummy_nodes(graph)
    new_graph.to_file(args.out_file_name)

def add_indel_nodes2(args):
    variants = GenotypeCalls.from_vcf(args.vcf_file_name)
    graph = Graph.from_file(args.graph_file_name)
    adder = DummyNodeAdder(graph, variants)
    new_graph = adder.create_new_graph_with_dummy_nodes()
    new_graph.to_file(args.out_file_name)

def add_allele_frequencies(args):
    logging.info("Reading graph")
    graph = Graph.from_file(args.graph_file_name)
    graph.set_allele_frequencies_from_vcf(args.vcf_file_name)
    graph.to_file(args.graph_file_name)
    logging.info("Wrote modified graph to the same file %s" % args.graph_file_name)


def make_haplotype_to_nodes(args):
    graph = Graph.from_file(args.graph_file_name)
    variants = GenotypeCalls.from_vcf(args.vcf_file_name)
    haplotype_to_nodes = HaplotypeToNodes.from_graph_and_variants(graph, variants, args.n_haplotypes)
    logging.info("Saving to file")
    haplotype_to_nodes.to_file(args.out_file_name)
    logging.info("Wrote to file %s" % args.out_file_name)


def main():
    run_argument_parser(sys.argv[1:])


def run_argument_parser(args):
    parser = argparse.ArgumentParser(
        description='Obgrapph.',
        prog='obgraph',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=100))

    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser("make")

    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-j", "--vg-json-files", nargs='+', required=True)
    subparser.set_defaults(func=make)

    subparser = subparsers.add_parser("add_indel_nodes2")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-g", "--graph-file-name", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.set_defaults(func=add_indel_nodes2)

    subparser = subparsers.add_parser("add_allele_frequencies")
    subparser.add_argument("-g", "--graph-file-name", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.set_defaults(func=add_allele_frequencies)

    subparser = subparsers.add_parser("make_haplotype_to_nodes")
    subparser.add_argument("-g", "--graph-file-name", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.add_argument("-n", "--n-haplotypes", type=int, required=True)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.set_defaults(func=make_haplotype_to_nodes)

    def make_node_to_haplotypes_lookup(args):
        haplotype_nodes = HaplotypeNodes.from_file(args.haplotype_nodes)
        n = NodeToHaplotypes.from_haplotype_nodes(haplotype_nodes)
        n.to_file(args.out_file_name)
        logging.info("Saved to %s" % args.out_file_name)

    subparser = subparsers.add_parser("make_node_to_haplotypes_lookup")
    subparser.add_argument("-H", "--haplotype_nodes", required=True)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.set_defaults(func=make_node_to_haplotypes_lookup)

    def make_genotype_matrix(args):
        from .genotype_matrix import GenotypeMatrix
        graph = Graph.from_file(args.graph)
        variants = GenotypeCalls.from_vcf(args.vcf_file_name, skip_index=True, limit_to_n_lines=None)

        if args.node_to_haplotypes is not None:
            nodes_to_haplotypes = NodeToHaplotypes.from_file(args.node_to_haplotypes)
            matrix = GenotypeMatrix.from_nodes_to_haplotypes_and_variants(nodes_to_haplotypes, variants, graph, args.n_individuals)
        else:
            logging.info("Making genotype matrix directly from vcf")
            matrix = GenotypeMatrix.from_variants(variants, args.n_individuals)

        matrix.to_file(args.out_file_name)

    subparser = subparsers.add_parser("make_genotype_matrix")
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.add_argument("-n", "--n-individuals", type=int, required=True)
    subparser.add_argument("-N", "--node_to_haplotypes", required=False)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.set_defaults(func=make_genotype_matrix)

    def analyse_genotype_matrix(args):

        matrix = GenotypeMatrix.from_file(args.genotype_matrix)
        analyser = GenotypeMatrixAnalyser(matrix)
        lookup = analyser.analyse()
        lookup.to_file(args.out_file_name)
        logging.info("Wrote lookup of most similar genotype to file %s" % args.out_file_name)


    subparser = subparsers.add_parser("analyse_genotype_matrix")
    subparser.add_argument("-G", "--genotype-matrix", required=True)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.set_defaults(func=analyse_genotype_matrix)

    def traverse(args):
        g = Graph.from_file(args.graph)
        haplotype_to_nodes = HaplotypeToNodes.from_file(args.haplotype_nodes)
        #from .traversing import traverse_graph_by_following_nodes

        for haplotype in range(0, args.n_haplotypes):
            nodes_to_follow = set(haplotype_to_nodes.get_nodes(haplotype))
            start_time = time.time()
            new_nodes = cython_traversing.traverse_graph_by_following_nodes(g, nodes_to_follow)
            logging.info("Got %d nodes" % len(new_nodes))
            end_time = time.time()
            logging.info("Time spent on haplotype %d: %.5f" % (haplotype, end_time - start_time))

    subparser = subparsers.add_parser("traverse")
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-T", "--type", required=False, default="correct_haplotype_nodes")
    subparser.add_argument("-H", "--haplotype_nodes", required=False)
    subparser.add_argument("-n", "--n_haplotypes", type=int, default=1, required=False)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.set_defaults(func=traverse)

    def get_genotype_frequencies(args):
        matrix = GenotypeMatrix.from_file(args.genotype_matrix)
        frequencies = GenotypeFrequencies.from_genotype_matrix(matrix)
        frequencies.to_file(args.out_file_name)
        logging.info("Wrote frequencies to file %s" % args.out_file_name)

    subparser = subparsers.add_parser("get_genotype_frequencies")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-g", "--genotype-matrix", required=True)
    subparser.set_defaults(func=get_genotype_frequencies)

    def make_random_haplotypes(args):
        graph = Graph.from_file(args.graph)
        variants = GenotypeCalls.from_vcf(args.vcf_file_name, skip_index=True)
        haplotype_nodes = HaplotypeToNodes.make_from_n_random_haplotypes(graph, variants, n_haplotypes=args.n_haplotypes)
        logging.info("Making new haplotypenodes by traversing full graph for each haplotype")
        new = haplotype_nodes.get_new_by_traversing_graph(graph, args.n_haplotypes)
        new.to_file(args.out_file_name)
        logging.info("Wrote haplotypenodes to %s" % args.out_file_name)

    subparser = subparsers.add_parser("make_random_haplotypes")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.add_argument("-n", "--n-haplotypes", type=int, required=False, default=10)
    subparser.set_defaults(func=make_random_haplotypes)


    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)
    args.func(args)

