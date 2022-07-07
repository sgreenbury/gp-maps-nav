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
./bazel-bin/src/landscape/nav_alg --help

# RNA, L=15 example
./bazel-bin/src/landscape/nav_alg \
    --base 4 \
    --dimension 15 \
    --length 15 \
    --file_name rna15 \
    --fitness_assignment target \
    --nd_ref 0 \
    --neutral_mutations \
    --gp_map rna \
    --pheno_fname gp_maps/RNA_15/pheno_list0.txt \
    --seed 1 \
    --sources $sources \
    --targets $targets \
    --threshold $threshold \
    --outpath scratch/

# Polyomino, S_{3,8} example
./bazel-bin/src/landscape/nav_alg \
    --base 8 \
    --dimension 12 \
    --length 12 \
    --file_name s38 \
    --fitness_assignment target \
    --nd_ref 1 \
    --neutral_mutations \
    --gp_map polyomino \
    --pheno_fname gp_maps/s_3_8/phenotypes/ \
    --assembly_tests 20 \
    --seed 1 \
    --sources $sources \
    --targets $targets \
    --threshold $threshold \
    --outpath scratch/
