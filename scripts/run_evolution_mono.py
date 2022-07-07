#!/usr/bin/env python

import argparse
import os
import subprocess
import uuid
import numpy as np


def get_unique_id(integer=False):
    """
    Use uuid1 or uuid4 from lib to generate unique ID
    https://docs.python.org/3/library/uuid.html
    """
    if not integer:
        return str(uuid.uuid4()).replace("-", "")
    # Return a seed that can be used as an unsigned long
    return str(uuid.uuid4().int)[:9]


def proc_com(procs):
    """
    Communicate with open processes in procs list, holding loop until
    process completes.
    """
    for proc in procs:
        proc.communicate()


def print_args(args):
    """Print args."""
    for arg, val in vars(args).items():
        try:
            print(f"{arg:>39} : {val:>39}")
        except Exception:
            pass


# Set-up argparser
parser = argparse.ArgumentParser()
parser.add_argument(
    "--targets", type=int, default=1, help="Number of targets to run"
)
parser.add_argument(
    "--sources", type=int, default=20, help="Number of targets to run"
)
parser.add_argument("--n_runs", type=int, default=20, help="Number of runs")
parser.add_argument("--generations", type=int, default=50000)
parser.add_argument("--base", type=int, default=4)
parser.add_argument("--niceness", type=str, default="-5")
parser.add_argument("--population_sizes", type=int, nargs="+", default=1000)
parser.add_argument("--lengths", type=int, nargs="+", default=20)
parser.add_argument("--shape_level", type=int, default=0)
parser.add_argument("--frna", dest="frna", action="store_true")
parser.add_argument("--nofrna", dest="frna", action="store_false")
parser.set_defaults(frna=True)
parser.add_argument(
    "--record_diffs_only", dest="record_diffs_only", action="store_true"
)
parser.add_argument(
    "--norecord_diffs_only", dest="record_diffs_only", action="store_false"
)
parser.set_defaults(record_diffs_only=True)

parser.add_argument(
    "--source_fitness", type=str, default="random", help="'random' or 'zero'"
)
parser.add_argument(
    "--fitness_measures", type=str, nargs="+", default="random"
)
parser.add_argument("--outpath_root", type=str, default="outs_rna_monomorphic")

parser.add_argument(
    "--neutral_mutations", dest="neutral_mutations", action="store_true"
)
parser.add_argument(
    "--noneutral_mutations", dest="neutral_mutations", action="store_false"
)
parser.set_defaults(neutral_mutations=True)

parser.add_argument(
    "--executable", type=str, default="./bazel-bin/evolution/nav_rna_mono"
)
parser.add_argument("--parallel_level", type=int, default=1)
args, unknown = parser.parse_known_args()

# Print args
print_args(args)

# Change path to location run from path calling
os.chdir("./")

# If fRNA, choose different geno and pheno list paths
if args.frna:
    gp_maps = [
        (f"frna/data/ps_l{length}.txt", f"frna/data/cs_l{length}.txt", length)
        for length in args.lengths
    ]
else:

    gp_maps = [
        ("gp_maps/RNA_12/pheno_list0.txt", "", 12),
        ("gp_maps/RNA_15/pheno_list0.txt", "", 15),
    ]

# Loop over length
for i_length, (pheno_fname, geno_fname, length) in enumerate(gp_maps):
    procs = []
    outfiles = []

    # Loop over fitnesses
    for i_fitness_measure, fitness_measure in enumerate(args.fitness_measures):

        # Loop over population sizes
        for i_population_size, population_size in enumerate(
            args.population_sizes
        ):

            # Make path for outfiles
            outpath = os.path.join(
                args.outpath_root,
                f"{fitness_measure}",
                f"l{length}",
                f"n{population_size}/",
            )
            # Make paths
            os.makedirs(outpath, exist_ok=True)

            for i_run, run in enumerate(np.arange(1, args.n_runs + 1)):

                # Get unique file ID
                file_name = get_unique_id()

                # Get seed
                seed = get_unique_id(integer=True)

                # Get the outfiles for stdout
                outfiles.append(
                    open(os.path.join(outpath, f"stdout_{file_name}.txt"), "w")
                )

                # Load genotypes from files as well as phenotypes for starting points
                load_genos = (
                    "--noload_genos" if geno_fname == "" else "--load_genos"
                )

                # Print to stdout
                command_str = [
                    "nice",
                    str(args.niceness),
                    str(args.executable),
                    "--verbose",
                    "--pheno_fname",
                    str(pheno_fname),
                    "--geno_fname",
                    str(geno_fname),
                    str(load_genos),
                    "--base",
                    str(args.base),
                    "--length",
                    str(length),
                    "--seed",
                    str(seed),
                    "--outpath",
                    str(outpath),
                    "--file_name",
                    str(file_name),
                    "--threshold",
                    str(args.generations),
                    "--fitness_measure",
                    str(fitness_measure),
                    "--sources",
                    str(args.sources),
                    "--targets",
                    str(args.targets),
                    "--shape_level",
                    str(args.shape_level),
                    (
                        "--neutral_mutations"
                        if args.neutral_mutations
                        else "--noneutral_mutations"
                    ),
                ]
                if args.record_diffs_only:
                    command_str += ["--record_diffs_only"]
                if int(population_size) >= 0:
                    command_str += [
                        "--population_size",
                        str(population_size),
                    ]

                print(" ".join(command_str))
                procs.append(
                    subprocess.Popen(
                        command_str,
                        stdout=outfiles[-1],
                    )
                )
            if args.parallel_level == 0:
                proc_com(procs)
        if args.parallel_level == 1:
            proc_com(procs)
    if args.parallel_level == 2:
        proc_com(procs)
