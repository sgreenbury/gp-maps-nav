#!/bin/bash

# NB. Script set to performm dry-run. Replace `--dry-run` with `--nodry-run` to run all

# Swaps and Dim for each GP map (Fig. 3A,3B,3C)
echo "> Running navigability for correlations and dimensions experiments experiments..."
for i in \
    "RNA_12","gp_maps/RNA_12/geno_list0.txt",4,12,0 \
    "s_2_8","gp_maps/s_2_8/geno_list0.txt",8,8,2 \
    "HP5x5s","gp_maps/HP5x5s/geno_list0_nd1.txt",2,25,0 \
    "HP3x3x3s","gp_maps/HP3x3x3s/geno_list0_nd1_reindex.txt",2,27,0;
do
    # Set parameters
    IFS=","
    set -- $i
    echo $1, $2, $3, $4, $5
    gp_map=$1
    geno_fname=$2
    base=$3
    length=$4
    nd_ref=$5

    # Run swaps script
    echo "> ${gp_map} correlations..."
    nohup \
	python scripts/run_nav_swaps_and_dim.py \
	--dry-run \
	--gp_map $gp_map \
	--experiment "correlations" \
	--seed_type "swaps" \
	--min_swaps -1 \
	--max_swaps -1 \
	--num_swaps 100 \
	--length $length \
	--base $base \
	--geno_fname  $geno_fname \
	--nd_ref $nd_ref \
	--min_dim $length \
	--max_dim $length \
	--num_dim 1 \
	--sources 20 \
	--targets 50 \
	--parallel_level 0 \
	> "nohup_${gp_map}_swaps.out" &
    # Wait for swaps to finish before running dimensions
    wait

    # Run dimensions script
    echo "> ${gp_map} dimensions..."
    nohup \
	python scripts/run_nav_swaps_and_dim.py \
	--dry-run \
	--gp_map $gp_map \
	--experiment "dimensionality" \
	--seed_type "dim" \
	--min_swaps 0 \
	--max_swaps 0 \
	--num_swaps 1 \
	--length $length \
	--base $base \
	--geno_fname  $geno_fname \
	--nd_ref $nd_ref \
	--min_dim 1 \
	--max_dim $length \
	--num_dim $length \
	--sources 20 \
	--targets 50 \
	--nostop_at_fittest \
	--parallel_level 0 \
	> "nohup_${gp_map}_dim.out" &
    wait

done

# Correlations properties for swaps cases: Swaps and Dim for each GP map
# (Fig. 3A, x-axis)
echo "> Running GP map stats for correlations experiments..."
for i in \
    "RNA_12","gp_maps/RNA_12/geno_list0.txt",4,12,0 \
    "s_2_8","gp_maps/s_2_8/geno_list0.txt",8,8,2 \
    "HP5x5s","gp_maps/HP5x5s/geno_list0_nd1.txt",2,25,0 \
    "HP3x3x3s","gp_maps/HP3x3x3s/geno_list0_nd1_reindex.txt",2,27,0;
do
    # Set parameters
    IFS=","
    set -- $i
    echo $1, $2, $3, $4, $5
    gp_map=$1
    geno_fname=$2
    base=$3
    length=$4
    nd_ref=$5

    echo "${gp_map} stats"

    # Run swaps script
    nohup \
	python scripts/run_nav_swaps_and_dim.py \
	--dry-run \
	--gp_map $gp_map \
	--executable_mode "gp_map_stats" \
	--seed_type "swaps" \
	--experiment "correlations" \
	--min_swaps -1 \
	--max_swaps -1 \
	--num_swaps 100 \
	--length $length \
	--base $base \
	--geno_fname  $geno_fname \
	--nd_ref $nd_ref \
	--parallel_level 0 \
	> "nohup_${gp_map}_swaps_pheno_stats.out" &
    # Wait for swaps to finish before running dimensions
    wait
done


# Mixed swaps and dimensions for RNA_12 (Fig. 3D,3E)
# Three dimensions (D=2,5,12) and three swaps (swaps=0, 4000000, 2**31-1) are run.
echo "${gp_map} mixed correlations and dimensions..."
for i in \
    "RNA_12","gp_maps/RNA_12/geno_list0.txt",4,12,0;
do
    # Set parameters
    IFS=","
    set -- $i
    echo $1, $2, $3, $4, $5
    gp_map=$1
    geno_fname=$2
    base=$3
    length=$4
    nd_ref=$5

    # Run swaps script
    nohup \
	python scripts/run_nav_swaps_and_dim.py \
	--dry-run \
	--gp_map $gp_map \
	--seed_type "random" \
	--swaps_spec 0 4000000 2147483647 \
	--dim_spec 2 5 12 \
	--length $length \
	--base $base \
	--geno_fname  $geno_fname \
	--nd_ref $nd_ref \
	--sources 20 \
	--targets 50 \
	--parallel_level 0 \
	> "nohup_RNA_12_swaps_and_dim.out" &
    # Wait for swaps to finish before running dimensions
    wait
done

# EOF
