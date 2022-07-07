import os
import subprocess
import argparse
import uuid


def get_unique_id(integer=False):
    """
    Use uuid1 or uuid4 from lib to generate unique ID
    https://docs.python.org/3/library/uuid.html
    """
    if not integer:
        return str(uuid.uuid4()).replace("-", "")
    # Return a seed that can be used as an unsigned long
    return str(uuid.uuid4().int)[:9]


def get_bys(args):
    by = args.by
    if by == "all":
        bys = ["source-target", "target"]
    else:
        bys = [by]
    return bys


# Make argument parser
parser = argparse.ArgumentParser(description="Analysis script.")
parser.add_argument("--by", default="target", type=str)
parser.add_argument(
    "--path", default="20211223_outs_rna_polymorphic/random/l20/", type=str
)
parser.add_argument("--mono_path", default="symlink_outs_mono", type=str)
parser.add_argument("--mono_pop_sizes", default="10,1000", type=str)
parser.add_argument(
    "--poly_path", default="symlinks/20220208_symlink_outs_poly", type=str
)
parser.add_argument("--outpath_root", default="analysis/scratch", type=str)
parser.add_argument("--threshold_stride", default=0.025, type=int)
parser.add_argument("--lengths", default="20,30,40", type=str)
parser.add_argument("--mono", action="store_true")
parser.add_argument("--nomono", action="store_false")
parser.set_defaults(mono=False)
args, unknown_args = parser.parse_known_args()


# Change path to project root where the script should be called from
os.chdir("./")
niceness = "-5"
python = "python"
executable = "scripts/evo_analysis/analysis_script.py"
path = args.path


# Set from argparse
bys = get_bys(args)
mono_path = args.mono_path
poly_path = args.poly_path
threshold_stride = args.threshold_stride
outpath_root = args.outpath_root

# Set fitness types to search over based on whether mono
if not args.mono:
    fitness_types = ["max_fit", "most_common", "mrca"]
    diff_thresholds = [0.0, 0.5]
    completed_thresholds = [0.0, 0.5]
else:
    fitness_types = ["mono"]
    diff_thresholds = [0.0]
    completed_thresholds = [0.0]

outfiles = []
for by in bys:
    for i_diff_threshold, diff_threshold in enumerate(diff_thresholds):
        for i_comp_threshold, completed_threshold in enumerate(
            completed_thresholds
        ):
            procs = []
            for i_comp_only, completed_only in enumerate([1, 0]):
                for i_fitness_type, fitness_type in enumerate(fitness_types):
                    file_name = get_unique_id()
                    seed = get_unique_id(integer=True)
                    outpath = os.path.join(
                        outpath_root,
                        f"fitness_type-{fitness_type}",
                        f"diff_thresh-{diff_threshold:.2f}",
                        f"comp_thresh-{completed_threshold:.2f}",
                        f"comp_only-{completed_only:.0f}",
                    )
                    os.makedirs(outpath, exist_ok=True)

                    outfiles.append(
                        open(
                            os.path.join(outpath, "out.txt"),
                            "w",
                        )
                    )

                    if fitness_type == "max_fit":
                        source_var = "Source"
                        target_var = "Target"
                        fitness_var = "Max_fitness"
                        proportion_var = "Proportion_max_fit"
                    elif fitness_type == "most_common":
                        source_var = "Source"
                        target_var = "Target"
                        fitness_var = "most_common_phenotype_fitness"
                        proportion_var = "most_common_phenotype_proportion"
                    elif fitness_type == "mrca":
                        source_var = "Source"
                        target_var = "Target"
                        fitness_var = "mrca_fitness"
                        proportion_var = None
                    elif fitness_type == "mono":
                        source_var = "Source_index"
                        target_var = "Target_index"
                        fitness_var = "Fitness_source"
                        proportion_var = None

                    # Make command string
                    command_str = [
                        "nice",
                        niceness,
                        python,
                        executable,
                        "--nomono_and_poly",
                        "--verbose",
                        "--path",
                        path,
                        "--poly_path",
                        poly_path,
                        "--mono_path",
                        mono_path,
                        "--outpath",
                        str(outpath),
                        "--fitness_var",
                        fitness_var,
                        "--source_var",
                        str(source_var),
                        "--target_var",
                        str(target_var),
                        "--by",
                        str(by),
                        "--by_mono",
                        str(by),
                        "--lengths",
                        args.lengths,
                        "--mono_pop_sizes",
                        args.mono_pop_sizes,
                    ]

                    if fitness_type in ["max_fit", "most_common"]:
                        command_str += [
                            "--diff_threshold",
                            f"{diff_threshold:.2f}",
                            "--completed_threshold",
                            f"{completed_threshold:.2f}",
                            "--proportion_var",
                            f"{proportion_var}",
                        ]
                    else:
                        command_str += []

                    if completed_only == 1:
                        command_str += [
                            "--completed_only",
                            "--completed_only_mono",
                        ]
                    else:
                        command_str += [
                            "--nocompleted_only",
                            "--nocompleted_only_mono",
                        ]

                    command_str += [
                        "--threshold_stride",
                        str(threshold_stride),
                    ]

                    print(command_str)

                    # Add to process list
                    procs.append(
                        subprocess.Popen(
                            command_str,
                            stdout=outfiles[-1],
                        )
                    )

            # Loop until complete
            for proc in procs:
                proc.communicate()
