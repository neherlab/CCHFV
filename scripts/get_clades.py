import json
from Bio import Phylo
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Get clades", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--node-data", type=str, required=True, help="output clade node data"
    )
    parser.add_argument("--clade-name", type=str, required=True, help="clade name")
    parser.add_argument("--tree", type=str, required=True, help="input newick tree")
    args = parser.parse_args()

    nodes = {}
    T = Phylo.read(args.tree, "newick")
    for node in T.find_clades():
        nodes[node.name] = {"clade_membership": args.clade_name}

    with open(args.node_data, "w") as fh:
        json.dump({"nodes": nodes}, fh)
