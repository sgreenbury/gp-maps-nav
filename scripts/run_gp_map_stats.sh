#!/bin/bash

# Print help
./bazel-bin/src/properties/gp_map_stats --help

# Run for all GP maps
# GP map properties for maps available on disk in Table 1/Fig. 2
echo "> Running GP map stats properties calculations (no component stats by default)..."
for i in \
    "RNA_12","gp_maps/RNA_12/geno_list0.txt",4,12,0 \
    "s_2_8","gp_maps/s_2_8/geno_list0.txt",8,8,2 \
    "HP5x5s","gp_maps/HP5x5s/geno_list0_nd1.txt",2,25,0 \
    "HP3x3x3s","gp_maps/HP3x3x3s/geno_list0_nd1_reindex.txt",2,27,0 \
    "HP_20","gp_maps/HP_20/geno_list0.txt",2,20,0 \
    "HP25","gp_maps/HP25/geno_list0.txt",2,25,0;
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
    outpath="data/raw_output/properties/${gp_map}/"
    mkdir -p $outpath

    nohup \
	./bazel-bin/src/properties/gp_map_stats \
	--base $base \
	--length $length \
	--file_name 0 \
	--nd_ref $nd_ref \
	--nd_include \
	--geno_fname $geno_fname \
	--seed 0 \
	--swaps 0 \
	--outpath $outpath \
	--estimate_component_stats \
	> "nohup_${gp_map}_overall_gp_map_stats.out" &
    wait
done
