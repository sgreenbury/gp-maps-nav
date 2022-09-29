# Navigability of genotype-phenotype (GP) maps
This code accompanies our Nature Ecology & Evolution
[paper](https://doi.org/10.1038/s41559-022-01867-z) and provides the
computational tools for investigating properties of genotype-phenotype (GP)
maps, specifically:
- Estimating GP map properties such as robustness and evolvability
- Estimate navigability of GP maps
- Simulate evolutionary dynamics

## Install
To install and compile the required tools, packages and binaries, follow the
instructions in [INSTALL.md](INSTALL.md).

## GP maps
### GP maps on disk
Utilising GP map models anlaysed and published in
[10.1371/journal.pcbi.1004773](https://doi.org/10.1371/journal.pcbi.1004773),
we generate the following small GP maps and store on disk to facilitate
faster computation:
  - `S_{2,8}`: Polyomino GP map with two tiles and eight colours
  - `RNA12`: RNA GP map with sequence length L=12
  - `HP5x5`: Compact HP GP map on a 5x5 grid (L=25)
  - `HP3x3x3`: Compact HP GP map on a 3x3x3 grid (L=27)
  - `HP20`: Non-compact HP GP map of length L=20
  - `HP25`: Non-compact HP GP map of length L=25

The GP maps are available after extracting:
```bash
tar -zxf gp_maps/gp_maps.tar.gz -C gp_maps
cd gp_maps
```

The GP maps have two files:
- `geno_list*txt`: line `n` is the phenotype of the the genotype that
   integer `n` is when represented in base `K` of the system.
- `pheno_list*txt`:  the integer `n` in `geno_list*.txt` corresponds to
  the phenotype that is written on `n`th line OR for polyomino the file
  `pheno*n*.txt` which contains a diagrammatic representation of the
  polyomino.

The following numerical maps are used to convert character sequences to base `K`
integer sequences:
```text
- S_{2,8} : {'0': 0, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7}
- RNA12   : {'A': 0, 'U': 1, 'C': 2, 'G': 3}
- HP5x5   : {'H': 0, 'P': 1}
- HP3x3x3 : {'H': 0, 'P': 1}
- HP20    : {'H': 1, 'P': 0}
- HP25    : {'H': 1, 'P': 0}
- fRNA    : {'C': 0, 'G': 1, 'A': 2, 'U': 3}
```

For example, the phenotype on line `12557964` of `gp_maps/RNA_12/geno_list0.txt`
is phenptype index `56`. This corresponds to phenotype on line `57` of
`gp_maps/RNA_12/pheno_list0.txt` which is `((.((...))))`. The mapping can be
expressed as:
```text
line of geno file =      12557964   ->  genotype index    =      12557963
genotype index    =      12557963   ->  base K=4 sequence = '320223123332'
base K=4 sequence = '320223123332'  ->  RNA sequence      = 'GCACCGUCGGGC'
RNA sequence      = 'GCACCGUCGGGC'  ->  RNA dot-bracket   = '((.((...))))'
```
The above mapping can be applied to all GP maps on disk in the same fashion
using the corresponding `char` to `int` maps described above.

### GP maps constructed algorithmically
The `S_{3,8}`, `RNA15` and `fRNA` GP maps are larger and the phenotype of a
given genotype is algorithmically estimated as it is not computationally
feasible to work with these larger GP maps on disk.

### fRNA extraction
The functional RNA database (fRNAdb) is available to
[download](https://doi.org/10.18908/lsdba.nbdc00452-001).

We extracted all genotype sequences of length `L` with `20 <= L <= 140` from the
file `frnadb_summary.csv`, outputting each sequence of a given `<LENGTH>` to
`fRNA/data/cs_l<LENGTH>.txt`. Each sequence in `fRNA/data/cs_l<LENGTH>.txt` was
then folded using ViennaRNA 1.8.5, with the dot-bracket structures saved in the
same line by line order to `fRNA/data/ps_l<LENGTH>.txt`.

## Experiments
### GP maps stats
The properties of the GP maps on disk can be measured with the binary:
```bash
./bazel-bin/src/properties/gp_map_stats --help
```

For example, with RNA12 and component level estimates:
```bash
./bazel-bin/src/properties/gp_map_stats \
  --geno_fname gp_maps/RNA_12/geno_list0.txt \
  --outpath scratch/ \
  --base 4 \
  --length 12 \
  --nd_ref 0 \
  --nd_include 0 \
  --nd_threshold 1 \
  --estimate_component_stats \
  --verbose
```

To measure properties given random genotype swaps (e.g. 100,000):
```
./bazel-bin/src/properties/gp_map_stats \
  --geno_fname gp_maps/RNA_12/geno_list0.txt \
  --outpath scratch/ \
  --base 4 \
  --length 12 \
  --nd_ref 0 \
  --swaps 100000 \
  --seed 1 \
  --verbose
```

[scripts/run_gp_map_stats.sh](scripts/run_gp_map_stats.sh) script to run the
experiments for Table 1 and Fig. 2.

### Navigability for GP maps on disk
Navigability can be estimated with any GP map in `gp_maps/` with:
```bash
./bazel-bin/src/landscape/nav_disk --help
```
For example, with RNA12:
```bash
./bazel-bin/src/landscape/nav_disk \
  --geno_fname gp_maps/RNA_12/geno_list0.txt \
  --outpath scratch/ \
  --base 4 \
  --length 12 \
  --targets 50 \
  --sources 20 \
  --nd_ref 0 \
  --seed 1 \
  --dimension 4 \
  --swaps 100000 \
  --verbose
```

### Navigability for GP maps not on disk
We estimate the navigability for larger GP maps that are not stored on disk by
algorithmically folding/assembling. These are available for `S_{3,8}` and
`RNA15` here:

```bash
./bazel-bin/src/landscape/nav_alg --help
```

E.g. for `S_{3,8}`:
```bash
./bazel-bin/src/landscape/nav_alg \
    --base 8 \
    --length 12 \
    --file_name s38 \
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
```

### Navigability for restricted neutral correlations and dimensionality
[scripts/run_nav_swaps_and_dim.sh](scripts/run_nav_swaps_and_dim.sh) provides a
bash script to perform the experiments exploring how navigability changes under
removing neutral correlations and restricting dimensionality is provided.

### Navigability in fRNA
A binary for performing the navigability estimation in fRNA is available running:
```bash
./bazel-bin/src/landscape/nav_rna_dfs --help
```

Complete bash scripts for performing all fRNA experiments are available at
[scripts/run_nav_rna_dfs.sh](scripts/run_nav_rna_dfs.sh). For example,
```bash
./scripts/run_nav_rna_dfs.sh 20000 hamming outs
```
will perform navigability estimation on fRNA at all lengths with and without
neutral mutations with a threshold of `20000` total genotypes considered with
Hamming fitness function applied.

### Evolutionary navigability
#### Monormorphic evolutionary dynamics
A binary to run monormorphic evolutionary dynamics on RNA is available with:
```bash
./bazel-bin/src/evolution/nav_rna_mono --help
```
An example run is given below:
```bash
./bazel-bin/src/evolution/nav_rna_mono \
	--length 40 \
	--outpath scratch/mono/ \
	--fitness_measure hamming \
	--threshold 10000 \
	--sources 1 \
	--targets 1 \
	--pheno_fname "frna/data/ps_l40.txt" \
	--population_size 10000 \
	--norecord_diffs_only \
	--seed 1 \
	--neutral_mutations \
	--shape_level 0
```
A complete script for running the experiments is available at
[scripts/run_all_mono.sh](scripts/run_all_mono.sh) which utilises a python
subscript [scripts/run_evolution_mono.py](scripts/run_evolution_mono.py).


#### Polymorphic evolutionary dynamics
A binary to run monormorphic evolutionary dynamics on RNA is available with:
```bash
./bazel-bin/src/evolution/nav_rna --help
```
An example run is given below:
```bash
./bazel-bin/src/evolution/nav_rna \
	--length 40 \
	--outpath scratch/poly/ \
	--generations 10000 \
	--sources 1 \
	--targets 1 \
	--pheno_fname "frna/data/ps_l40.txt" \
	--population_size 100 \
	--seed 1
```
A complete script for running the experiments is available at
[scripts/run_all_poly.sh](scripts/run_all_poly.sh) which utilises a python
subscript
[scripts/run_evolution_polymorphic.py](scripts/run_evolution_polymorphic.py).


#### Analysis
The outputs of the evolutionary dynamics are in the form of sequential mutants
(monomorphic case) or population summary statistics (polymorphic case). To derive
estimates of the navigability additional analysis must be run on the output files
to measure navigability in each source-target pair.

Scripts to produce these outputs across batches of paths containing evolutionary
ouputs can be run with, for example:
```bash
python scripts/evo_analysis/run_analysis_script_batch.py\
       --paths \
	   data/raw_output/mono \
	   data/raw_output/poly \
       --executable scripts/evo_analysis/run_analysis_script.py \
       --by all
```
calling the subscripts
[scripts/evo_analysis/run_analysis_script.py](scripts/evo_analysis/run_analysis_script.py)
and
[scripts/evo_analysis/analysis_script.py](scripts/evo_analysis/analysis_script.py).
with the outputs written within each of `--paths` to subpath `analysis/`.


### Plotting
The plotting notebooks and scripts used to generate the figures and tables are
available in [plotting](plotting/).
