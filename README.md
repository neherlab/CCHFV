# Nextstrain build for Crimean-Congo hemorrhagic fever virus

This repository provides scripts for the phylogenetic analysis of the CCHF virus.

## Installation
To run the workflow using TreeKnit you need to set up an environment containing Nextstrain and Julia with TreeKnit and TreeTools. The workflow can then be build from within this environment. 

## Get the data
Download the raw data of CCHFV sequences e.g. from the NCBI Virus database.

## Filter and group the data
For filtering the data the python script `scripts/filter_rawdata.py` can be used. This will create separate fasta files and metadata files for each segment of the CCHF virus in the `data`  folder, which are needed as input for the workflow. 

## Run the Workflow
There are three different snakefile versions provided within this repository. 

The standard snakefile called `Snakefile` contains the whole workflow including a resolving step using TreeKnit. To use TreeKnit it is necessary to call a Julia script and thus Julia has to be installed in the environment. The workflow could be run e.g. within a micromamba environment containing Nextstrain, Julia and TreeKnit with the command: 
```
nextstrain build --ambient .
```
If julia is not available or TreeKnit should be excluded from the workflow the snakefile `Snakefile_without_treeknit` can be used alternatively, by using the file option of snakemake: 
```
-s FILE
```

The third snakefile `Snakefile_without_linking` creates the phylogenetic trees of all segments separately without ensuring that the sequences for each segment are chosen from the same set of samples, making it harder to create tanglegrams of the trees. 
