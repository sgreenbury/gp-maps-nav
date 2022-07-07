#!/bin/bash

threshold=$1
fitness_measure=$2
outpath=$3

# fRNA neutrals
nohup python scripts/run_evolution_mono.py \
      --executable "./bazel-bin/src/landscape/nav_rna_dfs" \
      --sources 20 \
      --targets 1 \
      --n_runs 50 \
      --lengths 20 25 30 35 40 \
      --population_sizes -1 \
      --generations $threshold \
      --fitness_measures $2 \
      --shape_level 0 \
      --parallel_level 0\
      --outpath_root \
      "data/raw_output/landscape/frna/${threshold}/${outpath}/" \
      --norecord_diffs_only \
      --neutral_mutations \
      --frna
wait

# fRNA, no neutral mutations
nohup python scripts/run_evolution_mono.py \
      --executable "./bazel-bin/src/landscape/nav_rna_dfs" \
      --sources 20 \
      --targets 1 \
      --n_runs 50 \
      --lengths 20 25 30 35 40 \
      --population_sizes -1 \
      --generations $threshold \
      --fitness_measures $2 \
      --shape_level 0 \
      --parallel_level 0\
      --outpath_root \
      "data/raw_output/landscape/frna/${threshold}/${outpath}/" \
      --norecord_diffs_only \
      --noneutral_mutations \
      --frna

# EOF
