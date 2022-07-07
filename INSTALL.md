# Install

## GP maps navigability library
Ensure the following are available:
  - C++ compiler (`g++`)
  - Python 3 (`python`)
  - [Conda](https://docs.conda.io/en/latest/miniconda.html) (not required but
    used below)
  - [Bazel](https://docs.bazel.build/versions/main/install-os-x.html) build tool
	  Simplest options by OS:
  	- macOS: `brew install bazel`
  	- linux: `conda install -c conda-forge bazel`

If conda is available, make a new environment:
```bash
conda create -y --name gp_maps_nav python=3.8
conda activate gp_maps_nav
```
And then install the requirements and build:
```bash
tar -zxf gp_maps/gp_maps.tar.gz -C gp_maps
export GP_MAPS_NAV=$(pwd)
mkdir scratch
python -m pip install -r requirements.txt
```

## ViennaRNA library
### Compile and install ViennaRNA 1.8.5
Source and compilation instructions are for [ViennaRNA](https://www.tbi.univie.ac.at/RNA/)
(`ViennaRNA 1.8.5`). Additional instructions for `ViennaRNA 2.5.0` are
[below](#viennarna-250) but note the current fold
[src/landscape/rna_fold.cpp](src/landscape/rna_fold.cpp) needs modification
to function with `ViennaRNA >= 2.0`.

To be able to make this on Ubuntu 16.04 successfully and macOS, the following
need to be run.

```bash
./configure \
    CC=gcc \
    --prefix=<MY_VIENNA_RNA_PATH> \
    CFLAGS="-std=gnu89 -g -O2 -fcommon" \
    --without-perl
make
make install
```
Specifically:
  - `--gnu89` needed otherwise error with extern c definitions
  - `--without-perl` required otherwise did not compile with recent `perl`
  - Needs the GNU compiler, not mac clang compiler. Therefore on macOS,
    set `CC=gcc-6` or something similar installed from `brew`

### Symlink
To link the ViennaRNA library, the following symlink is required in the project
root:

```bash
cd $GP_MAPS_NAV
ln -s <MY_VIENNA_RNA_PATH> ViennaRNA
```

## RNAshapes binary
`RNAshapes` currently only correctly interfaces when running on linux. This
means that only `--shape-level 0` can be passed to binaries with the flag
`--shape-level`.

### Linux
To use RNA shapes as phenotypes, RNA shapes binary is most easily installed with
conda in a seperate environment:
```bash
conda install -y -c bioconda rnashapes
```
with the binary then available from the command line:
```
RNAshapes
```

### macOS
TODO

## Build
To build all binaries:
```bash
bazel build //...
```

## Test
The following should exit ok if all parts of the installation have been
successful:
```bash
bazel test //...
```

## Optional: compile and install ViennaRNA >= 2.0
### ViennaRNA 2.5.0
After downloading the [source](https://www.tbi.univie.ac.at/RNA/), the following can be run to install on:
- macOS Catalina:
    ```bash
    ./configure \
        --prefix=<MY_VIENNA_RNA_PATH> \
        CC=gcc-6 \
        CFLAGS="-std=gnu99 -g -O2" \
        --without-perl
    make
    make install
    ```
    `-std=gnu99` is required to compile for this version.

- Ubuntu 16.04
    ```bash
    ./configure \
        --prefix=<MY_VIENNA_RNA_PATH> \
        CC=gcc-5 \
        CXX=g++-5 \
        --without-perl
    make
    make install
    ```
    compiles ok but without specifying `gcc-5` Link Time Optiization (LTO)
    breaks the build on version being used.

- Ubuntu 20.04
    ```bash
    ./configure \
        --prefix=<MY_VIENNA_RNA_PATH>
    make
    make install
    ```
