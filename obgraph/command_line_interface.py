import logging
logging.basicConfig(level=logging.INFO, format='%(module)s %(asctime)s %(levelname)s: %(message)s')
import sys
import argparse
from . import Graph
from .util import add_indel_dummy_nodes

def make(args):
    logging.info("Will create from files %s" % args.vg_json_files)
    graph = Graph.from_vg_json_files(args.vg_json_files)
    graph.to_file(args.out_file_name)

def add_indel_nodes(args):
    graph = Graph.from_file(args.graph_file_name)
    new_graph = add_indel_dummy_nodes(graph)
    new_graph.to_file(args.out_file_name)

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

    subparser = subparsers.add_parser("add_indel_nodes")
    subparser.add_argument("-o", "--out_file_name", required=True)
    subparser.add_argument("-g", "--graph-file-name", required=True)
    subparser.set_defaults(func=add_indel_nodes)


    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)
    args.func(args)

