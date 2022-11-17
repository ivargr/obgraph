import logging
import pickle

logging.basicConfig(level=logging.INFO, format='%(module)s %(asctime)s %(levelname)s: %(message)s')
from obgraph.cython_traversing import traverse_graph_by_following_nodes
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
from shared_memory_wrapper.shared_memory import from_shared_memory, to_shared_memory, SingleSharedArray, remove_shared_memory_in_session, to_file, from_file, get_shared_pool, close_shared_pool
from multiprocessing import Pool
import time
from itertools import repeat
from .util import create_coordinate_map
from .util import fill_zeros_with_last
from .variant_to_nodes import VariantToNodes
from npstructures import RaggedArray




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
    variants = VcfVariants.from_vcf(args.vcf_file_name)
    graph = Graph.from_file(args.graph_file_name)
    adder = DummyNodeAdder(graph, variants)
    new_graph = adder.create_new_graph_with_dummy_nodes()
    edge_mapping = adder.get_edge_mapping()
    new_graph.to_file(args.out_file_name)
    logging.info("Chromosome start nodes in new graph: %s" % new_graph.chromosome_start_nodes)
    logging.info("Wrote new graph to file %s" % args.out_file_name)
    with open(args.out_file_name + ".edge_mapping", "wb") as f:
        pickle.dump(edge_mapping, f)
        logging.info("Wrote edge mapping from old edges to new dummy nodes to file %s.edge_mapping" % args.out_file_name)

def add_allele_frequencies(args):
    logging.info("Reading graph")
    graph = Graph.from_file(args.graph_file_name)
    variants = VcfVariants.from_vcf(args.vcf_file_name, limit_to_chromosome=args.chromosome, skip_index=True)
    graph.set_allele_frequencies_from_variants(variants, use_chromosome=args.chromosome)  # Use chromosome 1 because we always assume this is a single-chromosome graph
    graph.to_file(args.graph_file_name)
    logging.info("Wrote modified graph to the same file %s" % args.graph_file_name)


def make_haplotype_to_nodes(args):
    get_shared_pool(args.n_threads)

    graph = Graph.from_file(args.graph_file_name)
    variants = VcfVariants.from_vcf(args.vcf_file_name, make_generator=True, skip_index=True)
    haplotype_to_nodes = HaplotypeToNodes.from_graph_and_variants(graph, variants, args.n_haplotypes, n_threads=args.n_threads)
    # todo: Option to add dummy haplotypes for "all variant nodes" and "all reference nodes" (two for each, to mimic an individual)

    #haplotype_to_nodes = haplotype_to_nodes.get_new_by_traversing_graph(graph, args.n_haplotypes)
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

    subparser = subparsers.add_parser("add_indel_nodes")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-g", "--graph-file-name", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.set_defaults(func=add_indel_nodes)

    subparser = subparsers.add_parser("add_allele_frequencies")
    subparser.add_argument("-g", "--graph-file-name", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.add_argument("-c", "--chromosome", required=False, help="If vcf contains multiple chromsomes, use this to limit to the chromosome that the graph is made from")
    subparser.set_defaults(func=add_allele_frequencies)

    subparser = subparsers.add_parser("make_haplotype_to_nodes")
    subparser.add_argument("-g", "--graph-file-name", required=True)
    subparser.add_argument("-v", "--vcf-file-name", required=True)
    subparser.add_argument("-n", "--n-haplotypes", type=int, required=True)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-t", "--n-threads", type=int, default=8, required=False)
    subparser.set_defaults(func=make_haplotype_to_nodes)


    def make_haplotype_to_nodes_bnp(args):
        from .haplotype_nodes import make_ragged_haplotype_to_nodes
        phased_genotype_matrix = from_file(args.phased_genotype_matrix).matrix
        variant_to_nodes = VariantToNodes.from_file(args.variant_to_nodes)

        if args.make_disc_backed:
            from .haplotype_nodes import DiscBackedHaplotypeToNodes
            logging.info("Making disc backed")
            result = DiscBackedHaplotypeToNodes.from_phased_genotype_matrix(phased_genotype_matrix, variant_to_nodes, args.out_file_name)
            result.to_file(args.out_file_name)
        else:
            result = make_ragged_haplotype_to_nodes(variant_to_nodes, phased_genotype_matrix, args.n_threads)
            to_file(result, args.out_file_name)

    subparser = subparsers.add_parser("make_haplotype_to_nodes_bnp")
    subparser.add_argument("-g", "--variant-to-nodes", required=True)
    subparser.add_argument("-v", "--phased-genotype-matrix", required=True)
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-t", "--n-threads", type=int, default=8, required=False)
    subparser.add_argument("-n", "--n-haplotypes", type=int, required=False)
    subparser.add_argument("-d", "--make-disc-backed", type=bool, required=False, default=False, help="Uses less memory")
    subparser.set_defaults(func=make_haplotype_to_nodes_bnp)


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
        from .genotype_matrix import GenotypeMatrix, PhasedGenotypeMatrix

        if args.make_phased_matrix:
            logging.info("Making phased matrix")
            if args.vcf_file_name.endswith(".txt.gz"):
                logging.info("Input is txt")
                matrix = PhasedGenotypeMatrix.from_txt(args.vcf_file_name, args.n_variants, args.n_individuals)
            else:
                matrix = PhasedGenotypeMatrix.from_vcf(args.vcf_file_name, args.n_variants, args.n_individuals)
            to_file(matrix, args.out_file_name)
            return


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
    subparser.add_argument("-p", "--make-phased-matrix", required=False, type=bool, default=False)
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

            nodes_to_follow = np.zeros(len(g.nodes), dtype=np.uint8)
            nodes_to_follow[haplotype_to_nodes.get_nodes(haplotype)] = 1
            start_time = time.time()
            new_nodes = traverse_graph_by_following_nodes(g, nodes_to_follow)
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
            if i % 10000 == 0:
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
        variants = VcfVariants.from_vcf(args.vcf, skip_index=True)
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


    def make_position_id(args):
        from .position_id import PositionId
        graph = Graph.from_file(args.graph)
        position_id = PositionId.from_graph(graph)
        to_file(position_id, args.out_file_name)


    subparser = subparsers.add_parser("make_position_id")
    subparser.add_argument("-g", "--graph", required=True)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.set_defaults(func=make_position_id)


    def get_haplotype_sequence(args):
        # will traverse graph and follow these nodes
        # gets sequence for each chromosome in graph
        # writes fasta with sequences
        variant_nodes = np.load(args.nodes)
        nodes_to_follow = np.zeros(len(args.graph.nodes), dtype=np.uint8)
        nodes_to_follow[variant_nodes] = 1
        path_nodes, chromosome_indexes = traverse_graph_by_following_nodes(args.graph, nodes_to_follow, True)

        chromosome_indexes.append(len(path_nodes))

        chromosome_chunks = list(zip(chromosome_indexes[0:-1], chromosome_indexes[1:]))
        logging.info("Chromosome chunks: %s" % chromosome_chunks)

        chromosome_id = 1  # assume chromosomes are sorted
        fasta_lines = []
        coordinate_maps = {}
        refpos_to_node_maps = {}  # mapping from a ref pos in the haplotype coordinate space to node

        chromosome_index = 0
        for start_index, end_index in chromosome_chunks:
            nodes = path_nodes[start_index:end_index]
            fasta_lines.append(">" + str(chromosome_id) + "\n")
            haplotype_sequence = args.graph.get_nodes_sequence(nodes)
            fasta_lines.append(haplotype_sequence)
            fasta_lines.append("\n")

            # create a coordinate-map, a lookup from path pos to approx linear ref pos in graph
            coordinate_maps[str(chromosome_id)] = create_coordinate_map(nodes, args.graph, chromosome_index)

            refpos_to_node = np.zeros(len(haplotype_sequence), np.uint32)
            offsets = np.cumsum(args.graph.nodes[nodes])
            refpos_to_node[0] = nodes[0]
            refpos_to_node[offsets[:-1]] = nodes[1:]
            refpos_to_node = fill_zeros_with_last(refpos_to_node)
            refpos_to_node_maps[str(chromosome_id)] = refpos_to_node

            chromosome_index += 1
            chromosome_id += 1

        # write all sequences to fasta
        with open(args.out_file_name + ".fa", "w") as f:
            f.writelines(fasta_lines)
            logging.info("Wrote sequences to %s" % (args.out_file_name + ".fa"))

        # also write nodes in path to file
        np.save(args.out_file_name + ".nodes", path_nodes)
        logging.info("Wrote nodes in haplotype path to %s" % args.out_file_name + ".nodes.npy")

        to_file(coordinate_maps, args.out_file_name + ".coordinate_maps")
        logging.info("Wrote coordinate maps to %s" % args.out_file_name + ".coordinate_maps")

        to_file(refpos_to_node_maps, args.out_file_name + ".refpos_to_node")
        logging.info("Wrote refpos to node map to %s" % args.out_file_name + ".refpos_to_node")


    subparser = subparsers.add_parser("get_haplotype_sequence")
    subparser.add_argument("-n", "--nodes", required=True)
    subparser.add_argument("-g", "--graph", required=True, type=Graph.from_file)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.set_defaults(func=get_haplotype_sequence)

    def from_gfa(args):
        from .gfa import create_graph_from_gfa_file
        graph = create_graph_from_gfa_file(args.gfa)
        graph.to_file(args.out_file_name)

    subparser = subparsers.add_parser("from_gfa")
    subparser.add_argument("-g", "--gfa", required=True)
    subparser.add_argument("-o", "--out-file-name", required=True)
    subparser.set_defaults(func=from_gfa)


    def convert_gfa_ids_to_numeric_command(args):
        from .gfa import convert_gfa_ids_to_numeric
        convert_gfa_ids_to_numeric(args.gfa, args.out_base_name)


    subparser = subparsers.add_parser("convert_gfa_ids_to_numeric")
    subparser.add_argument("-g", "--gfa", required=True)
    subparser.add_argument("-o", "--out-base-name")
    subparser.set_defaults(func=convert_gfa_ids_to_numeric_command)

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)
    args.func(args)
    remove_shared_memory_in_session()

