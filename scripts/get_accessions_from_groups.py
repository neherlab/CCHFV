import pandas as pd
import argparse


def create_accession_list(metadata_name, groups_name, output_name):
    groups = pd.read_csv(groups_name, sep="\t")
    metadata = pd.read_csv(metadata_name, sep="\t")
    selected = metadata.loc[metadata["group_id"].isin(groups["group_id"].values)]
    accessions = list(selected["accession"])
    with open(output_name, "w") as file:
        for element in accessions:
            file.write(element + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="create accession_list",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--metadata",
        type=str,
        required=True,
        help="name of metadata of segment to get accessions from",
    )
    parser.add_argument(
        "--groups",
        type=str,
        required=True,
        help="name of group_ids file for which accessions should be selected",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="name of output text file with accession numbers",
    )
    args = parser.parse_args()

    create_accession_list(args.metadata, args.groups, args.output)
