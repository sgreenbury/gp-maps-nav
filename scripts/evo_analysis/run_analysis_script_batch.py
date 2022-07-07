#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
from datetime import datetime
from pathlib import Path
import argparse


def get_list_of_args(list_of_paths):
    list_of_args = []
    for path in list_of_paths:
        for el in Path(path).glob("**/*"):
            if not el.is_file() and len(list(el.glob("*yaml"))) > 0:
                path = os.path.join(el)
                outpath_root = os.path.join(el, "analysis")
                list_of_args.append([path, outpath_root])
    return list_of_args


# Parser for list of paths
parser = argparse.ArgumentParser()
parser.add_argument("--executable", type=str, default="run_analysis_script.py")
parser.add_argument(
    "--paths",
    type=str,
    default="../data/raw_output/20220220_outs_rna_mul_1/",
    nargs="+",
)
parser.add_argument("--by", type=str, default="all")
parser.add_argument("--dry-run", dest="dry_run", action="store_true")
parser.add_argument("--nodry-run", dest="dry_run", action="store_false")
parser.set_defaults(dry_run=False)
args, unknown = parser.parse_known_args()


# Call script from project root
os.chdir("./")

# Get list of paths from argparse
list_of_paths = args.paths

# Get list of args
list_of_args = get_list_of_args(list_of_paths)


# Add to process list
for path, outpath_root in list_of_args:
    mono = "--mono" if "mono" in outpath_root else "--nomono"
    print("---")
    print(datetime.now().strftime("%Y-%m-%d | %H:%M"))
    print(f"path: {path}")
    print(f"outpath_root: {outpath_root}")
    command_str = [
        "python",
        args.executable,
        "--outpath_root",
        outpath_root,
        "--path",
        path,
        "--by",
        args.by,
        mono,
    ]
    print(" ".join(command_str))
    if not args.dry_run:
        print("...")
        subprocess.run(command_str, capture_output=True)
