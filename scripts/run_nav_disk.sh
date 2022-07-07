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
    --sources $sources \
    --targets $targets \
    --threshold $threshold \
    --outpath scratch/

# Polyomino S_{2,8}
./bazel-bin/src/landscape/nav_disk \
    --base 4 \
    --dimension 8 \
    --length 8 \
    --file_name s28 \
    --fitness_assignment target \
    --nd_ref 0 \
    --neutral_mutations \
    --gp_map rna \
    --geno_fname gp_maps/s_2_8/geno_list0.txt \
    --run 0 \
    --seed 1 \
    --sources $sources \
    --targets $targets \
    --threshold $threshold \
    --outpath scratch/


# HP compact: HP5x5
./bazel-bin/src/landscape/nav_alg \
    --base 2 \
    --dimension 25 \
    --length 25 \
    --file_name hp5x5 \
    --fitness_assignment target \
    --nd_ref 0 \
    --neutral_mutations \
    --geno_fname gp_maps/HP5x5s/geno_list0.txt \
    --run 0 \
    --seed 1 \
    --sources $sources \
    --targets $targets \
    --threshold $threshold \
    --outpath scratch/

# HP compact: HP3x3x3
./bazel-bin/src/landscape/nav_alg \
    --base 2 \
    --dimension 27 \
    --length 27 \
    --file_name hp3x3x3 \
    --fitness_assignment target \
    --nd_ref 0 \
    --neutral_mutations \
    --geno_fname gp_maps/HP3x3x3s/geno_list0.txt \
    --run 0 \
    --seed 1 \
    --sources $sources \
    --targets $targets \
    --threshold $threshold \
    --outpath scratch/

# HP non-compact: HP20
./bazel-bin/src/landscape/nav_alg \
    --base 2 \
    --dimension 20 \
    --length 20 \
    --file_name hp20 \
    --fitness_assignment target \
    --nd_ref 0 \
    --neutral_mutations \
    --geno_fname gp_maps/HP_20/geno_list0.txt \
    --run 0 \
    --seed 1 \
    --sources $sources \
    --targets $targets \
    --threshold $threshold \
    --outpath scratch/

# HP non-compact: HP25
./bazel-bin/src/landscape/nav_alg \
    --base 2 \
    --dimension 25 \
    --length 25 \
    --file_name hp25 \
    --fitness_assignment target \
    --nd_ref 0 \
    --neutral_mutations \
    --geno_fname gp_maps/HP25/geno_list0.txt \
    --run 0 \
    --seed 1 \
    --sources $sources \
    --targets $targets \
    --threshold $threshold \
    --outpath scratch/
