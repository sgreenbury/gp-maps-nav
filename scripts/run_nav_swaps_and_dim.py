#!/usr/bin/env python
"""Script for running all correlations, dimensionality and ruggedness experiments."""

import os
import subprocess
import argparse
import numpy as np
from utils import proc_com, get_unique_id, print_args


def get_spec(args, swaps=False, dim=False):
    """Get an np.array of items to loop over given args."""
    assert not (swaps and dim)
    if swaps:
        if args.swaps_spec is None:
            return np.linspace(
                args.min_swaps, args.max_swaps, args.num_swaps, dtype="int"
            )
        return np.array(args.swaps_spec)
    if dim:
        if args.dim_spec is None:
            return np.linspace(
                args.min_dim, args.max_dim, args.num_dim, dtype="int"
            )
        return np.array(args.dim_spec)
    else:
        raise ValueError("swaps or dim must be True")


# Change path to project root where the script should be called from
os.chdir("./")
niceness = "-5"
python = "python"

parser = argparse.ArgumentParser(
    description="Navigability for varying correlations and dimensions."
)
parser.add_argument("--gp_map", default="RNA_12", type=str)
parser.add_argument("--executable_mode", default="nav", type=str)
parser.add_argument("--seed_type", default="swaps", type=str)
parser.add_argument("--base", default=4, type=int)
parser.add_argument("--length", default=12, type=int)
parser.add_argument("--sources", default=20, type=int)
parser.add_argument("--targets", default=50, type=int)
parser.add_argument("--threshold", default=1000000, type=int)
parser.add_argument(
    "--swaps_spec",
    default=None,
    type=int,
    nargs="+",
    help="List of swaps to run.",
)
parser.add_argument("--num_swaps", default=1, type=int)
parser.add_argument("--min_swaps", default=0, type=int)
parser.add_argument("--max_swaps", default=0, type=int)
parser.add_argument(
    "--dim_spec",
    default=None,
    type=int,
    nargs="+",
    help="List of dimensions to run.",
)
parser.add_argument(
    "--num_dim",
    default=1,
    type=int,
    help="Number of dimension restrictions to test.",
)
parser.add_argument("--min_dim", default=12, type=int)
parser.add_argument("--max_dim", default=12, type=int)
parser.add_argument("--nd_ref", default=0, type=int)
parser.add_argument(
    "--experiment",
    default="correlations",
    type=str,
    help="'correlations' or 'dimensions'",
)
parser.add_argument("--geno_fname", default=None, type=str)
parser.add_argument("--seed", default=1, type=int)
parser.add_argument(
    "--stop_at_fittest", dest="stop_at_fittest", action="store_true"
)
parser.add_argument(
    "--nostop_at_fittest", dest="stop_at_fittest", action="store_false"
)
parser.set_defaults(stop_at_fittest=True)
parser.add_argument("--dry-run", dest="dry_run", action="store_true")
parser.add_argument("--nodry-run", dest="dry_run", action="store_false")
parser.set_defaults(dry_run=False)
parser.add_argument("--parallel_level", type=int, default=0)
args, unknown_args = parser.parse_known_args()

# Get executable
if args.executable_mode == "nav":
    executable = "./bazel-bin/src/landscape/nav_disk"
elif args.executable_mode == "gp_map_stats":
    executable = "./bazel-bin/src/properties/gp_map_stats"
else:
    raise ValueError(f"{args.executable_mode} not implemented.")

# Set swaps and dim if -1 used
args.min_swaps = 0 if args.min_swaps == -1 else args.min_swaps
args.max_swaps = (
    args.base ** args.length if args.max_swaps == -1 else args.max_swaps
)
args.min_dim = 0 if args.min_dim == -1 else args.min_dim
args.max_dim = args.length if args.max_dim == -1 else args.max_dim

# Print args
print_args(args)

# Outpath
if args.executable_mode == "nav":
    outpath_root = (
        f"data/raw_output/landscape/{args.experiment}/{args.gp_map}/"
    )
elif args.executable_mode == "gp_map_stats":
    outpath_root = (
        f"data/raw_output/landscape/{args.experiment}/{args.gp_map}/NC/"
    )
else:
    raise ValueError(f"{args.executable_mode} not implemented.")

# Make outpath
os.makedirs(outpath_root, exist_ok=True)

# List for stdouts and processes
procs, stdouts = [], []

# Run scripts
for i, swaps in enumerate(get_spec(args, swaps=True), start=0):
    for j, dim in enumerate(get_spec(args, dim=True), start=0):
        # Generate seed and file_name
        # Use the enumerate index for swaps/dim so GP map stats and navigability
        # experiments can be seeded with the same value and have the same GP map after
        # swaps
        if args.seed_type == "swaps":
            seed = i + 1
            file_name = str(i)
        elif args.seed_type == "dim":
            seed = j + 1
            file_name = str(j)
        elif args.seed_type == "random":
            seed = get_unique_id(integer=True)
            file_name = get_unique_id(integer=False)
        else:
            raise ValueError(f"{args.seed_type} not implemented.")

        # Get commmand for execution dependent on mode
        if args.executable_mode == "nav":
            # Write command string for nav case
            cmd_str = [
                "nice",
                niceness,
                executable,
                "--base",
                str(args.base),
                "--length",
                str(args.length),
                "--file_name",
                # str(i),
                str(file_name),
                "--fitness_assignment",
                "target",
                "--nd_ref",
                str(args.nd_ref),
                "--neutral_mutations",
                "--geno_fname",
                str(args.geno_fname),
                # "--run",
                # str(run),
                "--seed",
                str(seed),
                "--dimension",
                str(dim),
                "--swaps",
                str(swaps),
                "--sources",
                str(args.sources),
                "--targets",
                str(args.targets),
                "--threshold",
                str(args.threshold),
                "--outpath",
                outpath_root,
            ]
            # If measuring ruggedness, do not stop when reach target, otherwise biased
            if not args.stop_at_fittest:
                cmd_str.append("--nostop_at_fittest")

        elif args.executable_mode == "gp_map_stats":
            # Write command string for nav case
            cmd_str = [
                "nice",
                niceness,
                executable,
                "--base",
                str(args.base),
                "--length",
                str(args.length),
                "--file_name",
                # str(i),
                str(file_name),
                "--nd_ref",
                str(args.nd_ref),
                "--geno_fname",
                str(args.geno_fname),
                # "--run",
                # str(run),
                "--seed",
                str(seed),
                "--swaps",
                str(swaps),
                # "--threshold",
                #  str(args.threshold),
                "--outpath",
                outpath_root,
            ]
            # If measuring ruggedness, do not stop when reach target, otherwise biased
            if not args.stop_at_fittest:
                cmd_str.append("--nostop_at_fittest")
        else:
            raise ValueError(f"{args.executable_mode} not implemented.")

        # Show command
        print(" ".join(cmd_str))

        # Run the command
        if not args.dry_run:
            # Get std out file object
            stdouts.append(
                open(
                    os.path.join(outpath_root, f"out_{file_name}.txt"),
                    "w",
                )
            )

            procs.append(
                subprocess.Popen(
                    cmd_str,
                    stdout=stdouts[-1],
                )
            )

        # Parallel over dim if parallel_level > 0
        if args.parallel_level == 0:
            proc_com(procs)

    # Parallel over swaps if parallel_level > 0
    if args.parallel_level == 0:
        proc_com(procs)
