import pandas as pd
from Bio import SeqIO
import argparse
import re
import numpy as np
import pathlib

segments = ['M', 'L', 'S']

def get_segment_data(sequences, metadata_df):
    records = SeqIO.parse(sequences, 'fasta')  
    segment_dic = {}
    for record in records: 
        found = False
        for segment in segments:
            re_input = re.compile('.*segment {0}.*'.format(segment), re.IGNORECASE)
            x = re_input.search(record.description)
            if x:
                segment_dic[record.id] = segment
                found = True
                break
        if not found:
            segment_dic[record.id] = np.nan
    segment_list = []
    for accession in metadata_df['Accession']:
        segment_list.append(segment_dic[accession])

    metadata_df.insert(2, "Segment", segment_list)


def group_metadata(df):
    min_length_s = 1000
    min_length_m = 4000
    min_length_l = 10000

    # Delete rows without date
    df_filtered = df.dropna(subset=['Isolate Collection date'])
    # Only rows with length > 200
    df_filtered = df_filtered.loc[df_filtered['Length'] > 250]

    # Delete rows without segments
    df_filtered = df.dropna(subset=['Segment'])

    # Only keep sequences with appropriate lengths
    df_filtered = df_filtered.drop(df_filtered[(df_filtered['Segment'] == 'S') & (df_filtered['Length'] < min_length_s)].index)
    df_filtered = df_filtered.drop(df_filtered[(df_filtered['Segment'] == 'M') & (df_filtered['Length'] < min_length_m)].index)
    df_filtered = df_filtered.drop(df_filtered[(df_filtered['Segment'] == 'L') & (df_filtered['Length'] < min_length_l)].index)

    # Create group_id
    df_grouped = df_filtered
    df_grouped['group_id']= 0

    # Group sequences according to isolate and collection date
    grouped = df_grouped.groupby(["Isolate Lineage", "Isolate Collection date"])
    groups = grouped.groups.keys()
    group_id = 1
    for g in groups: 
        isolate, date = g
        df_grouped['group_id']=df_grouped.apply(lambda row:row['group_id']+group_id if row["Isolate Lineage"] == isolate and row["Isolate Collection date"] == date else row['group_id'], axis=1)
        group_id += 1

    number_of_groups = group_id
    # Remove groups with group_id = 0 -> no isolate given
    df_grouped = df_grouped.loc[df_grouped['group_id'] != 0]

    # Add tag if all segments are present
    all_Seg = ['S','M','L']
    for g in range(1, number_of_groups+1): 
        group_g = df_grouped.loc[df_grouped['group_id'] == g]
        if 'S' in group_g['Segment'].values and 'M' in group_g['Segment'].values and 'L' in group_g['Segment'].values:
            df_grouped.loc[df_grouped['group_id'] == g, 'nr_segments'] = 'all'
            
    # Rename columns
    df_grouped = df_grouped.rename(columns= {"Virus Name":"virus", "Accession": "accession", "Isolate Collection date": "date",
                                            "Geographic Region": "region", "Geographic Location": "country"} )

    # Write to directory: metadata only containing sequences where all segments are present
    path = pathlib.Path('data')
    path.mkdir(parents=True, exist_ok=True)
    df_grouped.to_csv('data/all_sequences_grouped.tsv', sep="\t")

    
def create_segment_metadata():
    df_all = pd.read_csv('data/all_sequences_grouped.tsv', sep="\t")
    #Segment S
    metadata_S = df_all.loc[df_all['Segment'] == 'S']
    metadata_S.to_csv('data/metadata_S.tsv', sep="\t")
    #Segment M
    metadata_M = df_all.loc[df_all['Segment'] == 'M']
    metadata_M.to_csv('data/metadata_M.tsv', sep="\t")
    #Segment L
    metadata_L = df_all.loc[df_all['Segment'] == 'L']
    metadata_L.to_csv('data/metadata_L.tsv', sep="\t")

def create_segment_fasta(raw_sequences):
    all_sequences = "data/all_sequences_renamed.fasta"
    metadata_S = pd.read_csv("data/metadata_S.tsv", sep='\t')       
    groups_S = metadata_S["group_id"].tolist()
    sequences_S = "data/sequences_S.fasta"
    metadata_M = pd.read_csv("data/metadata_M.tsv", sep='\t')       
    groups_M = metadata_M["group_id"].tolist()
    sequences_M = "data/sequences_M.fasta"
    metadata_L = pd.read_csv("data/metadata_L.tsv", sep='\t')       
    groups_L = metadata_L["group_id"].tolist()
    sequences_L = "data/sequences_L.fasta"

    with open(raw_sequences) as raw, open(all_sequences, 'w') as renamed:
        records = SeqIO.parse(raw_sequences, 'fasta')
        for record in records:
            record.description = record.id
            SeqIO.write(record, renamed, 'fasta')

    with open(all_sequences) as all_Seq, open(sequences_S, 'w') as seq_S:
        records = SeqIO.parse(all_sequences, 'fasta')    
        for record in records: 
            seq = metadata_S.loc[metadata_S['accession'] == record.id]
            group = list(seq['group_id'])
            if len(group) > 0 and group[0] in groups_S:
                name = str(group[0])
                record.id = name
                record.description = record.id
                SeqIO.write(record, seq_S, 'fasta')

    with open(all_sequences) as all_Seq, open(sequences_M, 'w') as seq_M:
        records = SeqIO.parse(all_sequences, 'fasta')    
        for record in records: 
            seq = metadata_M.loc[metadata_M['accession'] == record.id]
            group = list(seq['group_id'])
            if len(group) > 0 and group[0] in groups_M:
                name = str(group[0])
                record.id = name
                record.description = record.id
                SeqIO.write(record, seq_M, 'fasta')
    
    with open(all_sequences) as all_Seq, open(sequences_L, 'w') as seq_L:
        records = SeqIO.parse(all_sequences, 'fasta')    
        for record in records: 
            seq = metadata_L.loc[metadata_L['accession'] == record.id]
            group = list(seq['group_id'])
            if len(group) > 0  and group[0] in groups_L:
                name = str(group[0])
                record.id = name
                record.description = record.id
                SeqIO.write(record, seq_L, 'fasta')

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="filter raw data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--metadata', type=str, required=True, help="csv file containing all sequences")
    parser.add_argument('--sequences', type=str, required=True, help="fasta file containing all sequences")
    args = parser.parse_args()
    df_metadata = pd.read_csv(args.metadata, sep="\t", on_bad_lines='warn')
    get_segment_data(args.sequences, df_metadata)

    group_metadata(df_metadata)
    create_segment_metadata()
    create_segment_fasta(args.sequences)