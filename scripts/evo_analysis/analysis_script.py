#!/usr/bin/env python
# coding: utf-8

"""Analysis script for evolutionary simulations."""
import argparse
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
from analysis import (
    get_diffs_and_completed,
    get_nav,
    get_stats,
    get_timestamp,
    plot_nav_and_aborted,
    plot_nav_and_aborted_by_target,
    print_args,
)


parser = argparse.ArgumentParser(allow_abbrev=False)
parser.add_argument("--lengths", default="20,30,40", type=str)
parser.add_argument("--figsize", default="5,4", type=str)
parser.add_argument(
    "--completed_only", dest="completed_only", action="store_true"
)
parser.add_argument(
    "--nocompleted_only", dest="completed_only", action="store_false"
)
parser.set_defaults(completed_only=True)
parser.add_argument(
    "--completed_only_mono", dest="completed_only_mono", action="store_true"
)
parser.add_argument(
    "--nocompleted_only_mono", dest="completed_only_mono", action="store_false"
)
parser.set_defaults(completed_only_mono=True)
parser.add_argument("--testing", dest="testing", action="store_true")
parser.add_argument("--notesting", dest="testing", action="store_false")
parser.set_defaults(testing=False)
parser.add_argument(
    "--mono_and_poly", dest="mono_and_poly", action="store_true"
)
parser.add_argument(
    "--nomono_and_poly", dest="mono_and_poly", action="store_false"
)
parser.set_defaults(mono_and_poly=True)
parser.add_argument("--mono_pop_sizes", default="10,1000", type=str)
parser.add_argument("--poly_pop_sizes", default="100", type=str)
parser.add_argument("--fitness_measures", default="random,hamming", type=str)
parser.add_argument(
    "--path", default="../20211223_outs_rna_polymorphic/random/l20/", type=str
)
parser.add_argument(
    "--mono_path", default="../20211122_outs_rna_mono/", type=str
)
parser.add_argument(
    "--poly_path", default="../20211223_outs_rna_polymorphic/", type=str
)
parser.add_argument("--source_var", default="Source", type=str)
parser.add_argument("--source_var_mono", default="Source_index", type=str)
parser.add_argument("--target_var", default="Target", type=str)
parser.add_argument("--target_var_mono", default="Target_index", type=str)
parser.add_argument(
    "--fitness_var",
    default="Max_fitness",
    type=str,
    help="'Max_fitness', 'most_common_phenotype_fitness', 'mrca_fitness'",
)
parser.add_argument(
    "--fitness_var_mono",
    default="Fitness_source",
    type=str,
    help="'Max_fitness', 'most_common_phenotype_fitness', 'mrca_fitness'",
)
parser.add_argument(
    "--proportion_var",
    default=None,
    type=str,
    help="None, 'Proportion_max_fit', 'most_common_phenotype_proportion'",
)
parser.add_argument(
    "--proportion_var_mono",
    default=None,
    type=str,
    help="None, 'Proportion_max_fit', 'most_common_phenotype_proportion'",
)
parser.add_argument(
    "--by", default="source-target", type=str, help="'source-target', 'target'"
)
parser.add_argument(
    "--by_mono",
    default="source-target",
    type=str,
    help="'source-target', 'target'",
)
parser.add_argument("--most_common_threshold", default=0.5, type=float)
parser.add_argument("--diff_threshold", default=0.0, type=float)
parser.add_argument("--completed_threshold", default=0.0, type=float)
parser.add_argument("--most_common_threshold_mono", default=0.5, type=float)
parser.add_argument("--diff_threshold_mono", default=0.0, type=float)
parser.add_argument("--completed_threshold_mono", default=0.0, type=float)
parser.add_argument("--outpath", default="scratch/", type=str)
parser.add_argument("--color1", default="dodgerblue", type=str)
parser.add_argument("--color2", default="firebrick", type=str)
parser.add_argument("--alpha", default=0.3, type=float)
parser.add_argument("--threshold_start", default=0.0, type=float)
parser.add_argument("--threshold_end", default=1.01, type=float)
parser.add_argument("--threshold_stride", default=0.01, type=float)
parser.add_argument("--dpi", default=150, type=int)
args, unknown = parser.parse_known_args()

if os.path.dirname(args.outpath) != "":
    os.makedirs(os.path.dirname(args.outpath), exist_ok=True)

print_args(args, width=50)


# Index for identifying run
run_idx = ["file_name", args.target_var, args.source_var]

# Get timestamp
ts = get_timestamp()

# Get generation level required statistics added on to stats
configs, stats, summary, mrca = get_stats(
    path=args.path,
    source_var=args.source_var,
    target_var=args.target_var,
    fitness_var=args.fitness_var,
    extra_cols=[
        "Generation",
        "most_common_phenotype_count",
        "Proportion_max_fit",
        "Max_fitness",
        "most_common_phenotype_fitness",
        "Target",
    ],
    testing=args.testing,
)

# Get the required diffs
stats = get_diffs_and_completed(
    stats,
    fitness_var=args.fitness_var,
    source_var=args.source_var,
    target_var=args.target_var,
    proportion_var=args.proportion_var,
    diff_threshold=args.diff_threshold,
    completed_threshold=args.completed_threshold,
)
# Get the aggregates
navs = get_nav(
    stats,
    target_var=args.target_var,
    source_var=args.source_var,
    threshold_start=args.threshold_start,
    threshold_end=args.threshold_end,
    threshold_stride=args.threshold_stride,
    completed_only=args.completed_only,
    by=args.by,
)


# Write to file
fname = os.path.join(args.outpath, f"single_{args.by}.csv")

if args.by == "source-target":
    navs.round(6).to_csv(fname, index=None)
else:
    rename_cols = (
        {"Target": "target_phenotype"}
        if "target_phenotype" not in stats.columns
        else {}
    )
    (
        navs.reset_index()
        .merge(
            stats.rename(columns=rename_cols)[
                ["file_name", args.target_var, "target_phenotype"]
            ].drop_duplicates(),
            how="left",
            on=["file_name", args.target_var],
        )
        .set_index(["file_name", args.target_var, "target_phenotype"])
        .round(6)
        .to_csv(fname)
    )


# Plot single case and write to file
fig, ax = plt.subplots(
    1, 1, figsize=tuple(map(float, args.figsize.split(",")))
)

if args.by == "source-target":
    ax, ax2 = plot_nav_and_aborted(
        navs=navs,
        ax=ax,
        color1=args.color1,
        color2=args.color2,
        alpha=args.alpha,
    )
elif args.by == "target":
    ax, ax2 = plot_nav_and_aborted_by_target(
        navs=navs,
        ax=ax,
        color1=args.color1,
        color2=args.color2,
        alpha=0.5,
    )
else:
    raise ValueError(f"Not implemented for '{args.by}'.")

ax.set_title(
    f"FV:{args.fitness_var}; PV:{args.proportion_var}; "
    f"diff_thresh: {args.diff_threshold}\n"
    f"completed_thresh:{args.completed_threshold}; completed_only:{args.completed_only}"
)

if args.by == "source-target":
    ax.set_ylim(-0.01, 1.01)
else:
    ax.set_ylim(
        -0.01, stats.drop_duplicates(subset=["Target", "file_name"]).shape[0]
    )

# Filename
fname = os.path.join(args.outpath, f"single_{args.by}")
plt.savefig(f"{fname}.pdf", transparent=True, bbox_inches="tight")
plt.savefig(
    f"{fname}.png", transparent=True, bbox_inches="tight", dpi=args.dpi
)
plt.show()


# Exit script if only plotting single case
if not args.mono_and_poly:
    sys.exit()


# Part 2: big panel plots for paper comparing monomorphic and polymorphic regime
# Get monomorhic results
fitness_measures = args.fitness_measures.split(",")
lengths = list(map(int, args.lengths.split(",")))
mono_pop_sizes = list(map(int, args.mono_pop_sizes.split(",")))

navs_list = []
for length in list(dict.fromkeys(lengths)):
    for fitness in list(dict.fromkeys(fitness_measures)):
        for pop_size in list(dict.fromkeys(mono_pop_sizes)):
            for mono in [True]:
                print(length, fitness, pop_size, mono)
                path = os.path.join(
                    args.mono_path, fitness, f"l{length}", f"n{pop_size}"
                )
                # Get generation level required statistics added on to stats
                configs, stats, summary, mrca = get_stats(
                    path=path,
                    source_var=args.source_var_mono,
                    target_var=args.target_var_mono,
                    fitness_var=args.fitness_var_mono,
                    extra_cols=[
                        "most_common_phenotype_count",
                        "Proportion_max_fit",
                    ],
                    testing=args.testing,
                )

                # Get the required diffs
                stats = get_diffs_and_completed(
                    stats,
                    fitness_var=args.fitness_var_mono,
                    source_var=args.source_var_mono,
                    target_var=args.target_var_mono,
                    proportion_var=args.proportion_var_mono,
                    diff_threshold=args.diff_threshold_mono,
                    completed_threshold=args.completed_threshold_mono,
                )

                # Get the aggregates
                navs = get_nav(
                    stats,
                    target_var=args.target_var_mono,
                    source_var=args.source_var_mono,
                    threshold_start=args.threshold_start,
                    threshold_end=args.threshold_end,
                    threshold_stride=args.threshold_stride,
                    completed_only=args.completed_only_mono,
                    by=args.by_mono,
                )
                navs["length"] = length
                navs["fitness"] = fitness
                navs["pop_size"] = pop_size
                navs["mono"] = mono
                navs_list.append(navs)
navs_mono = pd.concat(navs_list, axis=0).reset_index(drop=True)


# Get polymorphic results
fitness_measures = args.fitness_measures.split(",")
lengths = list(map(int, args.lengths.split(",")))
poly_pop_sizes = list(map(int, args.poly_pop_sizes.split(",")))
navs_list_poly = []
most_common_threshold = args.most_common_threshold
for length in list(dict.fromkeys(lengths)):
    for fitness in list(dict.fromkeys(fitness_measures)):
        for pop_size in list(dict.fromkeys(poly_pop_sizes)):
            for mono in [False]:
                print(length, fitness, pop_size, mono)
                path = os.path.join(args.poly_path, fitness, f"l{length}")

                # Get generation level required statistics added on to stats
                configs, stats, summary, mrca = get_stats(
                    path=path,
                    source_var=args.source_var,
                    target_var=args.target_var,
                    fitness_var=args.fitness_var,
                    extra_cols=[
                        "most_common_phenotype_count",
                        "Proportion_max_fit",
                    ],
                    testing=args.testing,
                )

                # Get the required diffs
                stats = get_diffs_and_completed(
                    stats,
                    fitness_var=args.fitness_var,
                    source_var=args.source_var,
                    target_var=args.target_var,
                    proportion_var=args.proportion_var,
                    diff_threshold=args.diff_threshold,
                    completed_threshold=args.completed_threshold,
                )

                # Get the aggregates
                navs = get_nav(
                    stats,
                    target_var=args.target_var,
                    source_var=args.source_var,
                    threshold_start=args.threshold_start,
                    threshold_end=args.threshold_end,
                    threshold_stride=args.threshold_stride,
                    completed_only=args.completed_only,
                    by=args.by,
                )
                navs["length"] = length
                navs["fitness"] = fitness
                navs["pop_size"] = pop_size
                navs["mono"] = mono

                # Append to large list
                navs_list_poly.append(navs)

navs_poly = pd.concat(navs_list_poly, axis=0).reset_index(drop=True)


# Plotting mono and poly combined
fitness_measures = args.fitness_measures.split(",")
lengths = list(map(int, args.lengths.split(",")))
pop_size = list(map(int, args.mono_pop_sizes.split(",")))[0]
max_y = 2
min_y = 0.0
# plot_completed = False
plot_completed = True

ylim = (0.0, 1.05)
ylim_by_target = (
    0.0,
    stats.drop_duplicates(subset=["Target", "file_name"]).shape[0],
)
ylim2 = (0.0, 1.05)
xlim = (0, 1.0)
fig, axs = plt.subplots(3, 3, figsize=(12, 9), sharex=False, sharey=True)
axs[0, 0].set_title(f"Monomorphic, n={pop_size:,.0f}")
for i, length in enumerate(lengths):
    ax = axs[i, 0]
    ax2 = ax.twinx()
    for j, fitness in enumerate(["random", "hamming"]):
        df = navs_mono[
            (navs_mono["length"] == length)
            & (navs_mono["fitness"] == fitness)
            & (navs_mono["pop_size"] == pop_size)
        ]

        # Labelling
        if i == 1:
            ylabel = True
        else:
            ylabel = False
        ylabel2 = False
        xlabel = False

        # Plot nav against threshold
        if args.by_mono == "source-target":
            ax, ax2 = plot_nav_and_aborted(
                df,
                ax,
                ax2,
                color1=f"C{j}",
                color2=f"C{j}",
                color_labels=False,
                ylabel=ylabel,
                ylabel2=ylabel2,
                xlabel=xlabel,
                ylim=ylim,
                ylim2=ylim2,
            )
        elif args.by == "target":
            ax, ax2 = plot_nav_and_aborted_by_target(
                df,
                ax,
                ax2,
                color1=f"C{j}",
                color2=f"C{j}",
                color_labels=False,
                ylabel=ylabel,
                ylabel2=ylabel2,
                xlabel=xlabel,
                ylim=ylim_by_target,
                ylim2=ylim2,
            )
        else:
            raise ValueError(f"Not implemented for '{args.by_mono}'.")

pop_size = list(map(int, args.mono_pop_sizes.split(",")))[1]
axs[0, 1].set_title(f"Monomorphic, n={pop_size:,.0f}")
for i, length in enumerate(lengths):
    ax = axs[i, 1]
    ax2 = ax.twinx()
    for j, fitness in enumerate(fitness_measures):
        df = navs_mono[
            (navs_mono["length"] == length)
            & (navs_mono["fitness"] == fitness)
            & (navs_mono["pop_size"] == pop_size)
        ]
        # Labelling
        if i == 2:
            xlabel = True
        else:
            xlabel = False
        ylabel2 = False
        ylabel = False

        # Plot nav against threshold
        if args.by_mono == "source-target":
            ax, ax2 = plot_nav_and_aborted(
                df,
                ax,
                ax2,
                color1=f"C{j}",
                color2=f"C{j}",
                color_labels=False,
                ylabel=ylabel,
                ylabel2=ylabel2,
                xlabel=xlabel,
                ylim=ylim,
                ylim2=ylim2,
            )
        elif args.by == "target":
            ax, ax2 = plot_nav_and_aborted_by_target(
                df,
                ax,
                ax2,
                color1=f"C{j}",
                color2=f"C{j}",
                color_labels=False,
                ylabel=ylabel,
                ylabel2=ylabel2,
                xlabel=xlabel,
                ylim=ylim_by_target,
                ylim2=ylim2,
            )
        else:
            raise ValueError(f"Not implemented for '{args.by_mono}'.")

pop_size = 100
axs[0, 2].set_title(f"Polymorphic, n={pop_size:,.0f}")
for i, length in enumerate(lengths):
    ax = axs[i, 2]
    ax2 = ax.twinx()
    for j, fitness in enumerate(fitness_measures):
        df = navs_poly[
            (navs_poly["length"] == length)
            & (navs_poly["fitness"] == fitness)
            # &(navs_mono["pop_size"]==pop_size)
        ]
        # Labelling
        if i == 1:
            ylabel2 = True
        else:
            ylabel2 = False
        ylabel = False
        xlabel = False

        # Plot nav against threshold
        if args.by == "source-target":
            ax, ax2 = plot_nav_and_aborted(
                df,
                ax,
                ax2,
                color1=f"C{j}",
                color2=f"C{j}",
                color_labels=False,
                ylabel=ylabel,
                ylabel2=ylabel2,
                xlabel=xlabel,
                ylim=ylim,
                ylim2=ylim2,
                label=fitness,
                label2=f"aborted ({fitness})",
            )
        elif args.by == "target":
            ax, ax2 = plot_nav_and_aborted_by_target(
                df,
                ax,
                ax2,
                color1=f"C{j}",
                color2=f"C{j}",
                color_labels=False,
                ylabel=ylabel,
                ylabel2=ylabel2,
                xlabel=xlabel,
                ylim=ylim_by_target,
                ylim2=ylim2,
                label=fitness,
                label2=f"aborted ({fitness})",
            )
        else:
            raise ValueError(f"Not implemented for '{args.by}'.")

        if i == 0:
            lines, labels = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax.legend(
                lines + lines2,
                labels + labels2,
                loc="upper left",
                bbox_to_anchor=(1.1, 1.0),
                prop={"size": "small"},
            )

# Add L=x labels
for i, length in enumerate(lengths):
    axs[i, 0].text(
        -0.3,
        1.05,
        f"L={length}",
        transform=axs[i, 0].transAxes,
        ha="left",
        weight="bold",
    )

# Set xlims
[ax.set_xlim(*xlim) for ax in axs[:, 0]]
[ax.set_xlim(*xlim) for ax in axs[:, 1]]
[ax.set_xlim(*xlim) for ax in axs[:, 2]]

# Set xticks
[ax.set_xticklabels([]) for ax in axs[0, :]]
[ax.set_xticklabels([]) for ax in axs[1, :]]

# Set ylims
# axs[0, 0].set_ylim(-0.01, 1.05)
# ax.set_xscale("log")

# Fig title
fig.suptitle(
    f'Fitness: "{args.fitness_var}";  '
    f'Proportion: "{args.proportion_var}";  '
    f"diff_thresh: {args.diff_threshold};  "
    f"completed_thresh: {args.completed_threshold};  "
    f"completed_only: {args.completed_only}",
    y=0.95,
    weight="bold",
    va="top",
    ha="center",
    size="medium",
)

# Filename
fname = os.path.join(args.outpath, f"panel_{args.by}")
plt.savefig(f"{fname}.pdf", transparent=True, bbox_inches="tight")
plt.savefig(
    f"{fname}.png", transparent=True, bbox_inches="tight", dpi=args.dpi
)
plt.show()
