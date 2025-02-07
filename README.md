# Nextstrain build for Crimean-Congo hemorrhagic fever virus

This repository provides scripts for the phylogenetic analysis of the CCHF virus.

## Installation

To run the workflow using TreeKnit you need to set up an environment containing [Nextstrain](https://docs.nextstrain.org/en/latest/install.html) and Julia with [TreeKnit](https://github.com/PierreBarrat/TreeKnit.jl) and TreeTools. The workflow can then be build from within this environment.

### WARNING
TreeKnit currently only runs on julia 1.7 - this does not exist micromamba for arm, 
```
micromamba create -n julia17_trial -c conda-forge julia=1.7 --platform osx-64
```
creates an environment with julia 1.7, however this isn't compatible with the other dependencies. I installed julia 1.7.3 via dmg and then compiled the treeknit package to run these scripts.

## Get the data

Download the raw data of CCHFV sequences e.g. from the NCBI Virus database.

## Filter and group the data

For filtering the data the python script `scripts/filter_rawdata.py` can be used. This will create separate fasta files and metadata files for each segment of the CCHF virus in the `data` folder, which are needed as input for the workflow.

## Run the Workflow

There are three different snakefile versions provided within this repository.

The standard snakefile called `Snakefile` contains the whole workflow including a resolving step using TreeKnit. To use TreeKnit it is necessary to call a Julia script and thus Julia has to be installed in the environment. The workflow could be run e.g. within a micromamba environment containing Nextstrain, Julia and TreeKnit with the command:

```
nextstrain build .
```

If julia is not available or TreeKnit should be excluded from the workflow the snakefile `Snakefile_without_treeknit` can be used alternatively, by using the file option of snakemake:

```
-s FILE
```

The third snakefile `Snakefile_without_linking` creates the phylogenetic trees of all segments separately without ensuring that the sequences for each segment are chosen from the same set of samples, making it harder to create tanglegrams of the trees.

## Local Development

For local development activate the micromamba environment using:

```
    micromamba create -f environment.yml --platform=linux/amd64
    micromamba activate cchfv
```

You can format your snakemake files using

```
snakefmt ${file}
```

You can visualize the workflow using:

```
snakemake -s ${file} --dag | dot -Tpng > dag.png
```
