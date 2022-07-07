#!/bin/bash

# Polymorphic evolutionary simulations: N=100; random and Hamming fitness
# assignments; L = 20, 30, 40; muL = 0.01, 0.1, 1
# mul = 1
nohup python scripts/run_evolution_polymorphic.py \
      --executable "./bazel-bin/src/evolution/nav_rna" \
      --sources 20 \
      --targets 10 \
      --targets_per_run 3 \
      --lengths 20 30 40 \
      --population_size 100 \
      --mul 1 \
      --generations 20000 \
      --fitness_measures random hamming \
      --parallel_level 1 \
      --source_fitness "random" \
      --outpath_root \
      "data/raw_output/rna_poly_mul_1" \
      --neutral_mutations \
      > nohup_mul_1.out &


# mul = 0.1
nohup python scripts/run_evolution_polymorphic.py \
      --executable "./bazel-bin/src/evolution/nav_rna" \
      --sources 20 \
      --targets 10 \
      --targets_per_run 3 \
      --lengths 20 30 40 \
      --population_size 100 \
      --mul 0.1 \
      --generations 20000 \
      --fitness_measures random hamming \
      --parallel_level 1 \
      --source_fitness "random" \
      --outpath_root \
      "data/raw_output/rna_poly_mul_01" \
      --neutral_mutations \
      > nohup_mul_01.out &

# mul = 0.01
nohup python scripts/run_evolution_polymorphic.py \
      --executable "./bazel-bin/src/evolution/nav_rna" \
      --sources 20 \
      --targets 10 \
      --targets_per_run 3 \
      --lengths 20 30 40 \
      --population_size 100 \
      --mul 0.01 \
      --generations 20000 \
      --fitness_measures random hamming \
      --parallel_level 1 \
      --source_fitness "random" \
      --outpath_root \
      "data/raw_output/rna_poly_mul_001" \
      --neutral_mutations \
      > nohup_mul_001.out &
wait

# EOF
