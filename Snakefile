rule all:
    input:
        #tree = expand("results/tree_{Segment}.nwk", Segment = glob_wildcards("./data/metadata_{Segment}.tsv").Segment)
        auspice_json = expand("auspice/CCHF_{Segment}_treeknit.json", Segment = glob_wildcards("./data/metadata_{Segment}.tsv").Segment)


input_fasta = "data/sequences_{Segment}.fasta",
input_metadata = "data/metadata_{Segment}.tsv",
input_fasta_index = "data/sequences_M.fasta",
input_metadata_index = "data/metadata_M.tsv",
dropped_strains = "config/dropped_strains.txt",
reference = "config/outgroup_{Segment}.gb",
colors = "config/colors.tsv",
lat_longs = "config/lat_longs.tsv",
auspice_config = "config/auspice_config.json"

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences = input_fasta_index
    output:
        sequence_index = "results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences = input_fasta_index,
        sequence_index = "results/sequence_index.tsv",
        metadata = input_metadata_index,
        exclude = dropped_strains
    output:
        filtered_metadata = "results/initial_filtered_metadata.tsv"
    params:
        group_by = "country year",
        sequences_per_group = 20,
        id_column = "group_id",
        min_date = 1900,
        query = "nr_segments == 'all'"
    
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output-metadata {output.filtered_metadata} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --metadata-id-columns {params.id_column} \
            --min-date {params.min_date} \
            --query {params.query:q}
        
        """

rule extract_selected_groups:
    message: 
        """
        Select the group id of filtered sequences
        """
    input:
        metadata = rules.filter.output.filtered_metadata
    output:
        filtered_groups_list = "results/filtered_groups.tsv",
        filtered_groups = "results/filtered_groups.txt"
    shell:
        """
        tsv-select -H -f group_id {input.metadata} > {output.filtered_groups_list}
        cp {output.filtered_groups_list} {output.filtered_groups}
        """

rule select_sequences:
    message: 
        """
        Select all sequences according to accessions list
        """
    input: 
        groups = rules.extract_selected_groups.output.filtered_groups,
        sequences = input_fasta
    output:
        selected_sequences = "results/selected_sequences_{Segment}.fasta"
    shell:
        """
        seqkit grep -f {input.groups} < {input.sequences} > {output.selected_sequences}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.select_sequences.output.selected_sequences,
        reference = reference
    output:
        alignment = "results/aligned_{Segment}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --remove-reference\
            --output {output.alignment} \
            --fill-gaps
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw_{Segment}.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --tree-builder-args="-czb" \
            --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = input_metadata
    output:
        tree = "results/tree_{Segment}.nwk",
        node_data = "results/branch_lengths_{Segment}.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        id_column = "group_id",
        root = "mid_point"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
            --metadata-id-columns {params.id_column}
        """

rule use_treeknit:
    message:
        """
        Using treeknit to detangle trees
        """
    input: 
        tree_S = "results/tree_S.nwk",
        tree_M = "results/tree_M.nwk",
        tree_L = "results/tree_L.nwk"
    output: 
        tree_S_res = "results/tree_S_resolved.nwk",
        tree_M_res = "results/tree_M_resolved.nwk",
        tree_L_res = "results/tree_L_resolved.nwk"
    shell:
        """
        julia scripts/treeknit_cli.jl \
            {input.tree_S} \
            {input.tree_M} \
            {input.tree_L} \
            --output_S {output.tree_S_res} \
            --output_M {output.tree_M_res} \
            --output_L {output.tree_L_res}
        """

    

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = "results/tree_{Segment}_resolved.nwk",
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts_{Segment}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = "results/tree_{Segment}_resolved.nwk",
        node_data = rules.ancestral.output.node_data,
        reference = reference
    output:
        node_data = "results/aa_muts_{Segment}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = "results/tree_{Segment}_resolved.nwk",
        metadata = input_metadata
    output:
        node_data = "results/traits_{Segment}.json",
    params:
        columns = "region country",
        id_column = "group_id"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --metadata-id-columns {params.id_column}
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = "results/tree_{Segment}_resolved.nwk",
        metadata = input_metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = colors,
        lat_longs = lat_longs,
        auspice_config = auspice_config
    output:
        auspice_json = "auspice/CCHF_{Segment}_treeknit.json",
    params: 
        id_column = "accession"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence \
            --output {output.auspice_json} \
            --metadata-id-columns {params.id_column}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
