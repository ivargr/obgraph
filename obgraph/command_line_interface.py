import logging
logging.basicConfig(level=logging.INFO, format='%(module)s %(asctime)s %(levelname)s: %(message)s')
import pyximport; pyximport.install()
import sys
import argparse
from . import Graph
from .util import add_indel_dummy_nodes
from .variants import VcfVariants
from .haplotype_nodes import HaplotypeToNodes, NodeToHaplotypes
from .dummy_node_adder import DummyNodeAdder
from .haplotype_nodes import NodeToHaplotypes
from .genotype_matrix import GenotypeMatrix, GenotypeMatrixAnalyser, GenotypeFrequencies
from pyfaidx import Fasta
from .graph_construction import GraphConstructor
from .graph_merger import merge_graphs
import numpy as np
from graph_kmer_index.shared_mem import from_shared_memory, to_shared_memory, SingleSharedArray, remove_shared_memory_in_session
from multiprocessing import Pool
import time
from itertools import repeat
from graph_kmer_index.flat_kmers import letter_sequence_to_numeric


def np_letter_sequence_to_numeric(letter_sequence):
    return letter_sequence_to_numeric(letter_sequence.astype("|S1").view(np.int8)).astype(np.uint8)

def get_numeric_node_sequence_single_thread(interval):
    from_pos, to_pos = interval
    start_time = time.time()
    graph = from_shared_memory(Graph, "graph_shared")
    numeric_node_sequences = from_shared_memory(SingleSharedArray, "numeric_node_sequences")
    result = np_letter_sequence_to_numeric(graph.node_sequences[from_pos:to_pos])
    numeric_node_sequences.array[from_pos:to_pos] = result
    logging.info("Spent %.3f s on interval" % (time.time()-start_time))
    return from_pos, to_pos


def merge_graphs_command(args):
    graphs = [Graph.from_file(graph) for graph in args.graphs]
    logging.info("Done reading graphs")

    merged_graph = merge_graphs(graphs)
    merged_graph.to_file(args.out_file_name)


def make(args):
    if args.vcf is not None:
        logging.info("Will create from vcf file")
        reference = Fasta(args.reference_fasta_file)

        chromosome = args.chromosome
        numeric_chromosome = chromosome
        if chromosome == "X":
            numeric_chromosome = "23"
        elif chromosome == "Y":
            numeric_chromosome = "24"

        ref_sequence = str(reference[args.chromosome])
        logging.info("Extracted sequence for chromosome %s. Length is: %d" % (chromosome, len(ref_sequence)))
        variants = VcfVariants.from_vcf(args.vcf, limit_to_chromosome=numeric_chromosome)
        logging.info("There are %d variants in chromosome" % len(variants))

        constructor = GraphConstructor(ref_sequence, variants)
        graph = constructor.get_graph_with_dummy_nodes()
        graph.to_file(args.out_file_name)
    else:
        logging.info("Will create from files %s" % args.vg_json_files)
        graph = Graph.from_vg_json_files(args.vg_json_files)
        graph.to_file(args.out_file_name)

def add_indel_nodes(args):
    graph = Graph.from_file(args.graph_file_name)
    new_graph = add_indel_dummy_nodes(graph)
    new_graph.to_file(args.out_file_name)

def add_indel_nodes2(args):
    variants = VcfVariants.from_vcf(args.vcf_file_name)
    graph = Graph.from_file(args.graph_file_name)
    adder = DummyNodeAdder(graph, variants)
    new_graph = adder.create_new_graph_with_dummy_nodes()
    new_graph.to_file(args.out_file_name)

def add_allele_frequencies(args):
    logging.info("Reading graph")
    graph = Graph.from_file(args.graph_file_name)
    variants = VcfVariants.from_vcf(args.vcf_file_name, limit_to_chromosome=args.chromosome, skip_index=True)
    graph.set_allele_frequencies_from_variants(variants, use_chromosome=1)  # Use chromosome 1 because we always assume this is a single-chromosome graph
    graph.to_file(args.graph_file_name)
    logging.info("Wrote modified graph to the same file %s" % args.graph_file_name)


def make_haplotype_to_nodes(args):
    graph = Graph.from_file(args.graph_file_name)
    variants = VcfVariants.from_vcf(args.vcf_file_name)
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
    subparser.add_argument("-j", "--vg-json-files", nargs='+', required=False)
    subparser.add_argument("-v", "--vcf", required=False)
    subparser.add_argument("-r", "--reference_fasta_file", required=False)
    subparser.add_argument("-c", "--chromosome", required=False)
    subparser.set_defaults(func=make)

    subparser = subparsers.add_parser("add_indel_nodes2")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-g", "--graph-file-name", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.set_defaults(func=add_indel_nodes2)

    subparser = subparsers.add_parser("add_allele_frequencies")
    subparser.add_argument("-g", "--graph-file-name", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.add_argument("-c", "--chromosome", default="1", required=False, help="If vcf contains multiple chromsomes, use this to limit to the chromosome that the graph is made from")
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

        if args.n_variants is not None:
            n_variants = args.n_variants
            assert args.n_individuals is not None, "n_individuals must be specified when n_variants is"
            n_individuals = args.n_individuals
        else:
            logging.info("N variants and individuals will be counted from the vcf file. Can take some time.")
            from .util import get_number_of_variants_and_individuals_from_vcf
            n_variants, n_individuals = get_number_of_variants_and_individuals_from_vcf(args.vcf_file_name)

        variants = VcfVariants.from_vcf(args.vcf_file_name, skip_index=True, limit_to_n_lines=None, make_generator=True)

        if args.node_to_haplotypes is not None:
            graph = Graph.from_file(args.graph)
            nodes_to_haplotypes = NodeToHaplotypes.from_file(args.node_to_haplotypes)
            matrix = GenotypeMatrix.from_nodes_to_haplotypes_and_variants(nodes_to_haplotypes, variants, graph, args.n_individuals)
        else:
            logging.info("Making genotype matrix directly from vcf")
            matrix = GenotypeMatrix.from_variants(variants, n_individuals, n_variants, n_threads=args.n_threads, chunk_size=args.chunk_size)

        matrix.to_file(args.out_file_name)

    subparser = subparsers.add_parser("make_genotype_matrix")
    subparser.add_argument("-g", "--graph", required=False)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.add_argument("-n", "--n-individuals", type=int, required=False)
    subparser.add_argument("-N", "--node-to-haplotypes", required=False)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.add_argument("-m", "--n-variants", required=False, type=int)
    subparser.add_argument("-t", "--n-threads", required=False, type=int, default=32, help="Number of threads used to fill matrix")
    subparser.add_argument("-c", "--chunk-size", required=False, type=int, default=10000, help="Number of variants to process in each job")
    subparser.set_defaults(func=make_genotype_matrix)

    def make_haplotype_matrix(args):
        from .haplotype_matrix import HaplotypeMatrix
        variants = VcfVariants.from_vcf(args.vcf_file_name, skip_index=True, limit_to_n_lines=None,
                                        make_generator=True)
        matrix = HaplotypeMatrix.from_variants(variants, args.n_individuals, args.n_variants,
                                                  n_threads=args.n_threads, chunk_size=args.chunk_size)

        matrix.to_file(args.out_file_name)

    subparser = subparsers.add_parser("make_haplotype_matrix")
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.add_argument("-n", "--n-individuals", type=int, required=True)
    subparser.add_argument("-N", "--node-to-haplotypes", required=False)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.add_argument("-m", "--n-variants", required=True, type=int)
    subparser.add_argument("-t", "--n-threads", required=False, type=int, default=6,
                           help="Number of threads used to fill matrix")
    subparser.add_argument("-c", "--chunk-size", required=False, type=int, default=10000,
                           help="Number of variants to process in each job")
    subparser.set_defaults(func=make_haplotype_matrix)

    def analyse_genotype_matrix(args):

        whitelist_array = None
        if args.whitelist_array is not None:
            whitelist_array = np.load(args.whitelist_array)

        matrix = GenotypeMatrix.from_file(args.genotype_matrix)
        analyser = GenotypeMatrixAnalyser(matrix, whitelist_array=whitelist_array)
        lookup = analyser.analyse(args.n_threads)
        lookup.to_file(args.out_file_name)
        logging.info("Wrote lookup of most similar genotype to file %s" % args.out_file_name)


    subparser = subparsers.add_parser("analyse_genotype_matrix")
    subparser.add_argument("-G", "--genotype-matrix", required=True)
    subparser.add_argument("-w", "--whitelist-array", required=False, help="Array of whitelist variants")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-t", "--n-threads", required=False, type=int, help="Number of threads to use", default=8)
    subparser.set_defaults(func=analyse_genotype_matrix)


    def make_transition_probabilities(args):
        from .genotype_matrix import GenotypeTransitionProbabilities, MostSimilarVariantLookup
        probs = GenotypeTransitionProbabilities.from_most_similar_variants_and_matrix(
            MostSimilarVariantLookup.from_file(args.most_similar_variants),
            GenotypeMatrix.from_file(args.genotype_matrix),
            n_threads=args.n_threads
        )
        probs.to_file(args.out_file_name)
        logging.info("Wrote to file %s" % args.out_file_name)

    subparser = subparsers.add_parser("make_transition_probabilities")
    subparser.add_argument("-G", "--genotype-matrix", required=True)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-m", "--most-similar-variants", required=True)
    subparser.add_argument("-t", "--n-threads", required=False, type=int, help="Number of threads to use", default=8)
    subparser.set_defaults(func=make_transition_probabilities)


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
        if args.vcf_file is not None:
            logging.warning("Creating naively from af field of vcf file")
            variants = VcfVariants.from_vcf(args.vcf_file)
            frequencies = GenotypeFrequencies.create_naive_from_vcf_af_field(variants)
        else:
            matrix = GenotypeMatrix.from_file(args.genotype_matrix)
            frequencies = GenotypeFrequencies.from_genotype_matrix(matrix, args.n_threads)

        frequencies.to_file(args.out_file_name)
        logging.info("Wrote frequencies to file %s" % args.out_file_name)

    subparser = subparsers.add_parser("get_genotype_frequencies")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-g", "--genotype-matrix", required=False)
    subparser.add_argument("-v", "--vcf-file", required=False, help="If specified, a naive approach will be used, computing from the AF field")
    subparser.add_argument("-t", "--n-threads", default=10, type=int, required=False)
    subparser.set_defaults(func=get_genotype_frequencies)

    def make_random_haplotypes(args):
        graph = Graph.from_file(args.graph)
        variants = VcfVariants.from_vcf(args.vcf_file_name, skip_index=True)
        haplotype_nodes = HaplotypeToNodes.make_from_n_random_haplotypes(graph, variants, n_haplotypes=args.n_haplotypes, weight_by_allele_frequency=not args.no_frequency_weighting)
        logging.info("Making new haplotypenodes by traversing full graph for each haplotype")
        new = haplotype_nodes.get_new_by_traversing_graph(graph, args.n_haplotypes)
        new.to_file(args.out_file_name)
        logging.info("Wrote haplotypenodes to %s" % args.out_file_name)

    subparser = subparsers.add_parser("make_random_haplotypes")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.add_argument("-n", "--n-haplotypes", type=int, required=False, default=10)
    subparser.add_argument("-e", "--no-frequency-weighting", type=bool, required=False, default=False, help="Set to True to not weight haplotypes by allele frequency")
    subparser.set_defaults(func=make_random_haplotypes)


    def validate_graph(args):
        variants = VcfVariants.from_vcf(args.vcf)
        graph = Graph.from_file(args.graph)

        for i, variant in enumerate(variants):
            if i % 1000 == 0:
                logging.info("%d variants processed" % i)

            ref_node, var_node = graph.get_variant_nodes(variant)

    subparser = subparsers.add_parser("validate_graph")
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-v", "--vcf", required=True)
    subparser.set_defaults(func=validate_graph)


    subparser = subparsers.add_parser("merge_graphs")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-g", "--graphs", nargs="+", required=True)
    subparser.set_defaults(func=merge_graphs_command)

    def make_variant_to_nodes(args):
        from .variant_to_nodes import VariantToNodes
        graph = Graph.from_file(args.graph)
        variants = VcfVariants.from_vcf(args.vcf)
        variant_to_nodes = VariantToNodes.from_graph_and_variants(graph, variants)
        variant_to_nodes.to_file(args.out_file_name)
        logging.info("Wrote to file %s" % args.out_file_name)

    subparser = subparsers.add_parser("make_variant_to_nodes")
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-v", "--vcf", required=True)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.set_defaults(func=make_variant_to_nodes)

    def make_node_to_variants(args):
        from .variant_to_nodes import NodeToVariants, VariantToNodes
        variant_to_nodes = VariantToNodes.from_file(args.variant_to_nodes)
        node_to_variants = NodeToVariants.from_variant_to_nodes(variant_to_nodes)
        node_to_variants.to_file(args.out_file_name)

    subparser = subparsers.add_parser("make_node_to_variants")
    subparser.add_argument("-v", "--variant_to_nodes", required=True)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.set_defaults(func=make_node_to_variants)

    def set_numeric_node_sequences(args):
        graph = Graph.from_file(args.graph)
        to_shared_memory(graph, "graph_shared")
        pool = Pool(args.n_threads)

        numeric_node_sequences = SingleSharedArray(np.zeros(len(graph.node_sequences), dtype=np.uint8))
        to_shared_memory(numeric_node_sequences, "numeric_node_sequences")


        intervals = list([int(i) for i in np.linspace(0, len(graph.node_sequences), args.n_threads+1)])
        intervals = [(from_pos, to_pos) for from_pos, to_pos in zip(intervals[0:-1], intervals[1:])]
        logging.info("Intervals: %s" % intervals)

        for from_pos, to_pos in pool.imap(get_numeric_node_sequence_single_thread, intervals):
            logging.info("Done processing interval %d-%d. Inserting into full array" % (from_pos, to_pos))

        logging.info("Done with all intervals. Saving new graph")
        numeric_node_sequences = from_shared_memory(SingleSharedArray, "numeric_node_sequences")
        graph.numeric_node_sequences = numeric_node_sequences.array
        graph.to_file(args.graph)
        logging.info("Saved to the same file %s" % args.graph)

    subparser = subparsers.add_parser("set_numeric_node_sequences")
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-t", "--n-threads", required=False, default=10, type=int)
    subparser.set_defaults(func=set_numeric_node_sequences)

    def create_coordinate_converter(args):
        from .coordinate_converter import CoordinateConverter
        converter = CoordinateConverter.from_graph(Graph.from_file(args.graph))
        converter.to_file(args.out_file_name)
        logging.info("Wrote to file %s" % args.out_file_name)

    subparser = subparsers.add_parser("create_coordinate_converter")
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.set_defaults(func=create_coordinate_converter)

    def intersect_vcfs(args):
        variants1 = VcfVariants.from_vcf(args.vcf1)
        variants2 = VcfVariants.from_vcf(args.vcf2)

        new = variants1.intersect(variants2)
        new.to_vcf_file(args.out_vcf, sample_name_output="")

    subparser = subparsers.add_parser("intersect_vcfs")
    subparser.add_argument("-a", "--vcf1", help="Vcf to intersect with other. Lines from this vcf will be kept")
    subparser.add_argument("-b", "--vcf2")
    subparser.add_argument("-o", "--out_vcf", help="Write resulting vcf to this file")
    subparser.set_defaults(func=intersect_vcfs)


    def make_numpy_variants(args):
        from .numpy_variants import NumpyVariants
        n = NumpyVariants.from_vcf(args.vcf)
        n.to_file(args.out_file_name)

    subparser = subparsers.add_parser("make_numpy_variants")
    subparser.add_argument("-v", "--vcf", required=True)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.set_defaults(func=make_numpy_variants)


    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)
    args.func(args)
    remove_shared_memory_in_session()

