#!/usr/bin/env python

import os
import subprocess
import numpy as np
import argparse

from utils import get_unique_id, proc_com, print_args


# Set-up argparser
parser = argparse.ArgumentParser()
parser.add_argument("-f")
parser.add_argument(
    "--targets", type=int, default=20, help="Number of targets to run"
)
parser.add_argument(
    "--sources", type=int, default=20, help="Number of targets to run"
)
parser.add_argument("--population_size", type=int, default=100)
parser.add_argument("--generations", type=int, default=20000)
parser.add_argument("--mu", type=float, default=None)
parser.add_argument("--mul", type=float, default=1.0)
parser.add_argument("--targets_per_run", type=int, default=1)
parser.add_argument("--niceness", type=str, default="-5")
parser.add_argument("--lengths", type=int, nargs="+", default=20)
parser.add_argument(
    "--source_fitness", type=str, default="random", help="'random' or 'zero'"
)
parser.add_argument(
    "--fitness_measures", type=str, nargs="+", default="random"
)
parser.add_argument("--outpath_root", type=str, default="outs_rna_polymorphic")
parser.add_argument("--track_mrca", dest="track_mrca", action="store_true")
parser.add_argument("--notrack_mrca", dest="track_mrca", action="store_false")
parser.set_defaults(track_mrca=True)
parser.add_argument(
    "--neutral_mutations", dest="neutral_mutations", action="store_true"
)
parser.add_argument(
    "--noneutral_mutations", dest="neutral_mutations", action="store_false"
)
parser.set_defaults(neutral_mutations=True)
parser.add_argument("--parallel_level", type=int, default=0)
parser.add_argument(
    "--executable", type=str, default="./bazel-bin/evolution/nav_rna"
)
args = parser.parse_args()

# Print args
print_args(args)

# Change path to project root where the script should be called from
os.chdir("./")

# Set-up from argparse
executable = args.executable
generations = args.generations
population_size = args.population_size
sources = args.sources
targets = args.targets
targets_per_run = args.targets_per_run
niceness = args.niceness
outpath_root = args.outpath_root
fitness_measures = args.fitness_measures
lengths = args.lengths
track_mrca = "--track_mrca" if args.track_mrca else "--notrack_mrca"
neutral_mutations = (
    "--neutral_mutations"
    if args.neutral_mutations
    else "--noneutral_mutations"
)
for i_length, length in enumerate(lengths):
    procs = []
    outfiles = []

    # Set mu from argparse if available
    if args.mu is None:
        # Default is to set: μL = 1; with pop size N = 100, μLN = 100 >> 1
        mu = args.mul / length
    else:
        mu = args.mu

    pheno_fname = (
        f"frna/data/ps_l{length}.txt"
        if length >= 20
        else f"gp_maps/RNA_{length}/pheno_list0.txt"
    )
    for i_fitness_measure, fitness_measure in enumerate(fitness_measures):
        outpath = os.path.join(
            outpath_root, f"{fitness_measure}", f"l{length}/"
        )
        os.makedirs(outpath, exist_ok=True)
        for i_target, pheno_target in enumerate(np.arange(targets)):
            file_name = get_unique_id()
            seed = get_unique_id(integer=True)

            outfiles.append(
                open(
                    os.path.join(outpath, f"stdout_{file_name}.txt"),
                    "w",
                )
            )
            command_str = [
                "nice",
                niceness,
                executable,
                "--verbose",
                "--pheno_fname",
                str(pheno_fname),
                "--length",
                str(length),
                "--seed",
                str(seed),
                "--outpath",
                str(outpath),
                "--file_name",
                str(file_name),
                "--generations",
                str(generations),
                "--fitness_measure",
                str(fitness_measure),
                "--sources",
                str(sources),
                "--targets",
                str(targets_per_run),
                "--population_size",
                str(population_size),
                "--source_fitness",
                str(args.source_fitness),
                "--mu",
                str(mu),
                str(track_mrca),
                str(neutral_mutations),
            ]
            print(command_str)
            procs.append(
                subprocess.Popen(
                    command_str,
                    stdout=outfiles[-1],
                )
            )

        # Parallel over targets
        if args.parallel_level == 0:
            proc_com(procs)

    # Parallel over fitness_measures
    if args.parallel_level == 1:
        proc_com(procs)

# Parallel over lengths
if args.parallel_level == 2:
    proc_com(procs)
