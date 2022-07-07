#!/bin/bash

# Run analysis script on all evolutionary outputs
python scripts/evo_analysis/run_analysis_script_batch.py\
       --paths\
       data/raw_output/rna_poly_mul_1 \
       data/raw_output/rna_poly_mul_01 \
       data/raw_output/rna_poly_mul_001 \
       data/raw_output/rna_mono \
       data/raw_output/rna_mono_no-neut-mut \
       data/raw_output/rna_mono_shape_l1 \
       data/raw_output/rna_mono_shape_l3 \
       --executable scripts/evo_analysis/run_analysis_script.py \
       --by all
