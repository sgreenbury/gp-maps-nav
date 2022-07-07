#!/bin/bash

RUN_TYPE="${1:-test}"
if [[ $RUN_TYPE == "test" ]]
then
    sources=1
    targets=1
    threshold=10000
else
    sources=20
    targets=50
    threshold=1000000
fi

# Print help
./bazel-bin/src/landscape/nav_disk --help

# RNA12
./bazel-bin/src/landscape/nav_disk \
    --base 4 \
    --dimension 12 \
    --length 12 \
    --file_name rna12 \
    --fitness_assignment target \
    --nd_ref 0 \
    --neutral_mutations \
    --gp_map rna \
    --geno_fname gp_maps/RNA_12/geno_list0.txt \
    --run 0 \
    --seed 1 \
    --swaps 100000 \
    --sources $sources \
    --targets $targets \
    --threshold $threshold \
    --outpath scratch/
