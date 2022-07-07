#!/bin/bash

# Shape level 0 (dot-bracket, with neutral mutations)
# L20, 30, 40; Pop size = 10, 1000, 100000
nohup python scripts/run_evolution_mono.py \
      --sources 20 \
      --targets 1 \
      --n_runs 50 \
      --lengths 20 30 40 \
      --population_sizes 10 1000 100000 \
      --generations 50000 \
      --fitness_measures random hamming \
      --shape_level 0 \
      --parallel_level 1\
      --outpath_root \
      "data/raw_output/rna_mono/" \
      --record_diffs_only \
      --neutral_mutations \
      --frna
wait

# Shape level 0 (dot-bracket, no neutral mutations)
# L20, 30, 40; Pop size = 10, 1000, 100000
nohup python scripts/run_evolution_mono.py \
      --sources 20 \
      --targets 1 \
      --n_runs 50 \
      --lengths 20 30 40 \
      --population_sizes 10 1000 100000\
      --generations 50000 \
      --fitness_measures random hamming \
      --shape_level 0 \
      --parallel_level 1\
      --outpath_root \
      "data/raw_output/rna_mono_no-neut-mut" \
      --record_diffs_only \
      --noneutral_mutations \
      --frna
wait

# Shape level 1
# L60, 100; Pop size = 10, 1000, 100000
nohup python scripts/run_evolution_mono.py \
      --sources 20 \
      --targets 1 \
      --n_runs 50 \
      --lengths 60 100 \
      --population_sizes 10 1000 100000 \
      --generations 1000 \
      --fitness_measures random lev \
      --shape_level 1 \
      --parallel_level 1\
      --outpath_root \
      "data/raw_output/rna_mono_shape_l1" \
      --record_diffs_only \
      --neutral_mutations \
      --frna
wait

# Shape level 3
# L60, 100, 140; Pop size = 10, 1000, 100000
nohup python scripts/run_evolution_mono.py \
      --sources 20 \
      --targets 1 \
      --n_runs 50 \
      --lengths 60 100 \
      --population_sizes 10 1000 100000 \
      --generations 1000 \
      --fitness_measures random lev \
      --shape_level 3 \
      --parallel_level 1\
      --outpath_root \
      "data/raw_output/rna_mono_shape_l3" \
      --record_diffs_only \
      --neutral_mutations \
      --frna

# EOF
