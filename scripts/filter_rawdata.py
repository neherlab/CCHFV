import pandas as pd
from Bio import SeqIO
import argparse
import re
import numpy as np
import pathlib

segments = ["M", "L", "S"]
min_length_dic = {"S": 1000, "M": 4000, "L": 10000}


def get_segment_data(sequences, metadata_df):
    records = SeqIO.parse(sequences, "fasta")
    segment_dic = {}
    for record in records:
        found = False
        for segment in segments:
            re_input = re.compile(".*segment {0}.*".format(segment), re.IGNORECASE)
            x = re_input.search(record.description)
            if x:
                segment_dic[record.id] = segment
                found = True
                break
        if not found:
            segment_dic[record.id] = np.nan
    segment_list = []
    for accession in metadata_df["Accession"]:
        segment_list.append(segment_dic[accession])

    metadata_df.insert(2, "Segment", segment_list)
    return metadata_df


def group_name(isolate, date, location):
    country = str(location).split(":")[0]
    return (
        str(isolate).replace(" ", "-")
        + "/"
        + str(date)
        + "/"
        + country.replace(" ", "-")
    )


def group_metadata(df):
    # Delete rows without date or segment
    df_filtered = df.dropna(subset=["Isolate Collection date"])
    df_filtered = df.dropna(subset=["Segment"])

    # Only keep sequences with appropriate lengths
    df_filtered = df_filtered.loc[df_filtered["Length"] > 250]
    for segment in segments:
        min_length = min_length_dic[segment]
        df_filtered = df_filtered.drop(
            df_filtered[
                (df_filtered["Segment"] == segment)
                & (df_filtered["Length"] < min_length)
            ].index
        )

    # Group sequences according to isolate and collection date
    df_grouped = df_filtered
    grouped = df_grouped.groupby(
        ["Isolate Lineage", "Isolate Collection date", "Geographic Location"]
    )
    groups = grouped.groups.keys()

    group_id = "None"
    df_grouped["group_id"] = group_id
    number_of_groups = 0
    for g in groups:
        isolate, date, location = g
        df_grouped["group_id"] = df_grouped.apply(
            lambda row: group_name(isolate, date, location)
            if row["Isolate Lineage"] == isolate
            and row["Isolate Collection date"] == date
            and row["Geographic Location"] == location
            else row["group_id"],
            axis=1,
        )
        number_of_groups += 1

    print("Number of groups: ", number_of_groups)

    # Remove groups with group_id = "None" -> no isolate given
    df_grouped = df_grouped.loc[df_grouped["group_id"] != "None"]

    # Add tag if all segments are present
    all_segments = 0
    for g in groups:
        isolate, date, location = g
        group_id = group_name(isolate, date, location)
        group_g = df_grouped.loc[df_grouped["group_id"] == group_id]
        if (
            "S" in group_g["Segment"].values
            and "M" in group_g["Segment"].values
            and "L" in group_g["Segment"].values
            and len(group_g) == 3
        ):
            df_grouped.loc[df_grouped["group_id"] == group_id, "nr_segments"] = "all"
            all_segments += 1
    print("Number of groups with all segments: ", all_segments)

    df_grouped = df_grouped.loc[df_grouped["nr_segments"] == "all"]

    # Rename columns
    df_grouped["country"] = df_grouped["Geographic Location"].apply(
        lambda x: x.split(":")[0] if isinstance(x, str) else None
    )
    df_grouped = df_grouped.rename(
        columns={
            "Virus Name": "virus",
            "Accession": "accession",
            "Isolate Collection date": "date",
            "Geographic Region": "region",
            "Submitter Names": "author"
        }
    )

    # Write to directory: metadata only containing sequences where all segments are present
    path = pathlib.Path("data")
    path.mkdir(parents=True, exist_ok=True)
    df_grouped.to_csv("data/all_sequences_grouped.tsv", sep="\t", index=False)
    return df_grouped


def write_segment_metadata(df_all):
    for segment in segments:
        df_segment = df_all.loc[df_all["Segment"] == segment]
        df_segment.to_csv(
            "data/metadata_{0}.tsv".format(segment), sep="\t", index=False
        )


def write_fasta(all_sequences_path, segment, df_all):
    metadata = df_all.loc[df_all["Segment"] == segment]
    groups = metadata["group_id"].tolist()
    sequences_segment_path = "data/sequences_{0}.fasta".format(segment)
    with open(sequences_segment_path, "w") as seq:
        records = SeqIO.parse(all_sequences_path, "fasta")
        for record in records:
            group = metadata.loc[metadata["accession"] == record.id]
            group_ids = list(group["group_id"])
            if len(group_ids) > 0 and group_ids[0] in groups:
                name = str(group_ids[0])
                record.description = record.id
                record.id = name
                SeqIO.write(record, seq, "fasta")


def write_segment_fasta(raw_sequences, df_all):
    all_sequences_path = "data/all_sequences_renamed.fasta"

    with open(all_sequences_path, "w") as renamed:
        records = SeqIO.parse(raw_sequences, "fasta")
        for record in records:
            record.description = record.id
            SeqIO.write(record, renamed, "fasta")
    for segment in segments:
        write_fasta(all_sequences_path, segment, df_all)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="filter raw data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--metadata", type=str, required=True, help="csv file containing all sequences"
    )
    parser.add_argument(
        "--sequences",
        type=str,
        required=True,
        help="fasta file containing all sequences",
    )
    args = parser.parse_args()
    df_metadata = pd.read_csv(args.metadata, sep="\t", on_bad_lines="warn")
    df_metadata = get_segment_data(args.sequences, df_metadata)

    df_grouped = group_metadata(df_metadata)
    write_segment_metadata(df_grouped)
    write_segment_fasta(args.sequences, df_grouped)
