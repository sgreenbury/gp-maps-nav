"""Helper functions for evolution analysis"""
import os
import re
import yaml
from collections import OrderedDict
import numpy as np
import datetime
import glob
import pandas as pd


def print_args(args, na_value=np.nan, width=40):
    """Print arguments from argparser."""
    for key, value in vars(args).items():
        try:
            print(f"{key:>{width}} : {value : >{width}}")
        except Exception:
            print(f"{key:>{width}} : {na_value : >{width}}")


def get_file_props(fname):
    file_category, file_name, file_type = re.findall(
        r".*/(.*)_(.*)\.(.*)$", fname
    )[0]
    return file_category, file_name, file_type


def load_files(
    files,
    fitness_var="Fitness_source",
    source_var="Source_index",
    target_var="Target_index",
    extra_cols=[],
    verbose=False,
    testing=False,
):
    if verbose:
        print(files[0])
    file_category, file_name, file_type = get_file_props(files[0])

    # Subset of cols to load
    if extra_cols is not None:
        usecols = [source_var, target_var, fitness_var] + extra_cols
        if file_category in ["stats", "mrca"]:
            usecols += ["Generation"]
        usecols = list(OrderedDict.fromkeys(usecols))
        # Only keep subset that exist in files
        cols_in_file = pd.read_csv(
            files[0], nrows=0, sep="\t"
        ).columns.to_list()
        usecols = [col for col in usecols if col in cols_in_file]
    else:
        usecols = None

    # Set-up frames list and testing condition
    dfs = []
    if testing:
        print("Testing...")
        files = sorted(files)[:2]
    # Reading yamls or txts
    if file_type in ["yaml", "yml"]:
        for file in files:
            with open(file, "r") as f:
                dfs.append(pd.json_normalize(yaml.safe_load(f)))
    else:
        for file in files:
            df = pd.read_csv(file, sep="\t", usecols=usecols)
            file_category, file_name, file_type = get_file_props(file)
            df["file_name"] = file_name
            dfs.append(df)
    return pd.concat(dfs, axis=0).reset_index(drop=True)


def get_stats(
    path=None,
    source_var="Source_index",
    target_var="Target_index",
    fitness_var="Fitness_source",
    extra_cols=[],
    verbose=False,
    testing=False,
):
    """
    Loads the stats and config dataframes from a path with unique filenames.
    """
    # Load configs
    configs = load_files(
        glob.glob(os.path.join(path, "config_*.yaml")),
        fitness_var=fitness_var,
        source_var=source_var,
        verbose=verbose,
        testing=testing,
    )

    # Load main per generation stats
    stats = load_files(
        glob.glob(os.path.join(path, "stats_*.txt")),
        fitness_var=fitness_var,
        source_var=source_var,
        target_var=target_var,
        extra_cols=extra_cols,
        verbose=verbose,
        testing=testing,
    )

    run_idx = ["file_name", target_var, source_var]
    stats = stats.sort_values(
        ["file_name", target_var, source_var, "Generation"]
    )

    # Load summary data if it exists
    summary_files = glob.glob(os.path.join(path, "summary_*.txt"))
    if len(summary_files) > 0:
        summary = load_files(
            summary_files,
            fitness_var=fitness_var,
            source_var=source_var,
            target_var=target_var,
            extra_cols=None,
            verbose=verbose,
            testing=testing,
        )
        stats = stats.merge(summary, on=run_idx, how="left")
    else:
        summary = None

    # Load mrca data if it exists
    mrca_files = glob.glob(os.path.join(path, "mrca_*.txt"))
    if len(mrca_files) > 0:
        run_idx_and_gen = run_idx + ["Generation"]
        mrca = load_files(
            mrca_files,
            fitness_var=fitness_var,
            source_var=source_var,
            target_var=target_var,
            extra_cols=None,
            verbose=verbose,
            testing=testing,
        )
        stats = stats.merge(mrca, on=run_idx_and_gen, how="left")

        # Fill forwards with LOCF for MRCA cols
        mrca_cols = mrca.drop(
            mrca.filter(run_idx_and_gen).columns, axis=1
        ).columns.to_list()
        stats.loc[:, mrca_cols] = stats.groupby(run_idx)[mrca_cols].fillna(
            method="ffill"
        )
    else:
        mrca = None

    # Make new summary variable for first generation at max fitness
    stats = stats.merge(
        stats[stats[fitness_var] == 1]
        .groupby(run_idx)["Generation"]
        .min()
        .rename("first_generation_at_max"),
        on=run_idx,
        how="left",
    )

    if "most_common_phenotype_count" in stats.columns:
        pop_sizes = configs["population_size"].unique()
        if len(pop_sizes) == 1:
            pop_size = pop_sizes[0]
        else:
            raise ValueError(f"{len(pop_sizes)} population_sizes present.")

        stats.loc[:, "most_common_phenotype_proportion"] = (
            stats["most_common_phenotype_count"].div(pop_size).round(5)
        )

    return configs, stats, summary, mrca


def get_diffs_and_completed(
    stats,
    source_var="Source_index",
    target_var="Target_index",
    fitness_var="Fitness_source",
    proportion_var=None,  # "most_common_phenotype_proportion",
    diff_threshold=None,  # 0.0,
    completed_threshold=None,  # 0.0,
):
    """
    Takes a full generation or transition by transtion dataframe and returns
    the same dataframe back but with fitness differences and completion stats
    calculated.
    """
    # Get copy
    stats = stats.copy()

    # Index for identifying run
    run_idx = ["file_name", target_var, source_var]

    # For a required "diff_threshold", calculate the delta fitness to previous
    # generation that satisfies condition
    if proportion_var is not None:
        diffs = (
            stats[stats[proportion_var] > diff_threshold]
            .groupby(run_idx)[fitness_var]
            .diff()
            .rename("delta_fitness")
        )
    else:
        diffs = (
            stats.groupby(run_idx)[fitness_var].diff().rename("delta_fitness")
        )
    stats.loc[diffs.index, "delta_fitness"] = diffs

    # Calculate cumulative minimum (per source-target)
    stats["cummin_delta_fitness"] = (
        stats.groupby(run_idx)["delta_fitness"]
        .expanding()
        .min()
        .reset_index(level=list(range(len(run_idx))), drop=True)
    )

    # Calculate cumulative maximum fitness (per source-target)
    if proportion_var is not None:
        stats.loc[
            stats[stats[proportion_var] > completed_threshold].index,
            f"{fitness_var}_completed_threshold",
        ] = stats[fitness_var].copy()
    else:
        stats.loc[stats.index, f"{fitness_var}_completed_threshold"] = stats[
            fitness_var
        ].copy()

    # Get max so far across valid fitnesses to consider
    stats[f"cummax_{fitness_var}"] = (
        stats.groupby(run_idx)[f"{fitness_var}_completed_threshold"]
        .expanding()
        .max()
        .reset_index(level=list(range(len(run_idx))), drop=True)
    )

    # For a required "compeleted_threshold", has it reached the target
    stats["completed"] = stats[f"cummax_{fitness_var}"].eq(1).mul(1)

    return stats


def get_nav(
    df,
    threshold_start=-1,
    threshold_end=1.01,
    threshold_stride=0.01,
    completed_only=False,
    rounding=3,
    source_var="Source_index",
    target_var="Target_index",
    by="source-target",
    verbose=False,
):
    """
    Takes a full generation or transition by transtion dataframe and returns
    useful statistics by "source-target" or "target" measuring navigability
    at different thresholds.
    """
    if verbose:
        print(df.head())

    # Index for identifying run
    run_idx = ["file_name", target_var, source_var]

    thresholds = np.arange(threshold_start, threshold_end, threshold_stride)

    # Get navigability
    # NB. Commented out np.isclose condition as slow and does not affect output
    navs = pd.DataFrame(
        (
            (
                df[["cummin_delta_fitness"]].to_numpy()
                >= -thresholds[..., np.newaxis].T
            )
            # | np.isclose(
            #     df[["cummin_delta_fitness"]].to_numpy(),
            #     -thresholds[..., np.newaxis].T,
            # )
            | np.isnan(df[["cummin_delta_fitness"]].to_numpy())
        )
        * df[["completed"]].to_numpy(),
        index=df.index,
        columns=[round(el, rounding) for el in thresholds],
    )
    navs.loc[:, [source_var, target_var, "file_name", "completed"]] = df[
        [source_var, target_var, "file_name", "completed"]
    ]

    # Get nav summary
    # This merge is a bottleneck but there are no apparent simple solutions
    navs_summary = (
        navs[run_idx]
        .drop_duplicates(subset=run_idx)
        .merge(
            navs[navs["completed"].eq(1)].drop_duplicates(
                subset=run_idx, keep="first"
            ),
            on=run_idx,
            how="left",
        )
        .fillna(0)
        .sort_values(run_idx)
        .reset_index(drop=True)
    )

    # Whether or not to exclude runs that will definitely fail in aborted calc
    if completed_only:
        aborted = pd.DataFrame(
            np.tile(
                df[["completed"]].eq(0).mul(1).to_numpy(), len(thresholds)
            ),
            index=df.index,
            columns=[round(el, rounding) for el in thresholds],
        )
    else:
        # NB. Commented out np.isclose condition as slow and does not affect output
        aborted = pd.DataFrame(
            (
                (
                    df[["cummin_delta_fitness"]].to_numpy()
                    >= -thresholds[..., np.newaxis].T
                )
                # | np.isclose(
                #     df[["cummin_delta_fitness"]].to_numpy(),
                #     -thresholds[..., np.newaxis].T,
                # )
                | np.isnan(df[["cummin_delta_fitness"]].to_numpy())
            )
            * (df[["completed"]].eq(0).mul(1)).to_numpy(),
            index=df.index,
            columns=[round(el, rounding) for el in thresholds],
        )
    aborted.loc[:, [source_var, target_var, "file_name", "completed"]] = df[
        [source_var, target_var, "file_name", "completed"]
    ]

    last_gen_aborted = (
        aborted.drop_duplicates(subset=run_idx, keep="last")
        .sort_values(run_idx)
        .reset_index(drop=True)
    )

    # Mask out aborted cases
    mask = last_gen_aborted.filter(regex=r"\d+\.\d+").eq(0)
    vals = navs_summary.filter(regex=r"\d+\.\d+")
    navs_summary.loc[:, vals.columns] = vals.where(mask, other=np.nan)

    # Return the per target distribution if by=="target"
    navs_by_threshold_by_target = (
        navs_summary.groupby(["file_name", target_var])
        .agg("mean")
        .drop(columns=source_var)
    )
    if by == "target":
        return navs_by_threshold_by_target

    # Aggregate nav and aborted
    navs_out = (
        navs_summary.agg(["mean", "sem"])
        .filter(regex=r"\d+\.\d+")
        .T.add_prefix("nav_")
        .rename_axis(index="threshold")
    ).join(
        last_gen_aborted.agg(["mean", "sem"])
        .filter(regex=r"\d+\.\d+")
        .T.add_prefix("aborted_")
        .rename_axis(index="threshold")
    )

    # Add completed aggregates
    navs_out["completed_mean"] = (
        navs_summary["completed"].agg(["mean", "sem"]).loc["mean"]
    )
    navs_out["completed_sem"] = (
        navs_summary["completed"].agg(["mean", "sem"]).loc["sem"]
    )

    # Bring threshold out as a variable
    navs_by_threshold = navs_out.reset_index()

    return navs_by_threshold


def get_timestamp():
    return datetime.datetime.utcnow().strftime("%Y-%m-%d-%H-%M-%S")


def plot_nav_and_aborted(
    navs,
    ax,
    ax2=None,
    color1="dodgerblue",
    color2="firebrick",
    alpha=0.3,
    xticklabels=True,
    yticklabels=True,
    yticklabels2=True,
    xlabel=True,
    ylabel=True,
    ylabel2=True,
    ylim=(0, 1),
    ylim2=(0, 1),
    color_labels=True,
    label=None,
    label2=None,
):
    if ax2 is None:
        ax2 = ax.twinx()

    ax.plot(
        navs["threshold"],
        navs["nav_mean"],
        color=color1,
        zorder=1,
        label=label,
    )
    ax.fill_between(
        navs["threshold"],
        navs["nav_mean"] - navs["nav_sem"],
        navs["nav_mean"] + navs["nav_sem"],
        lw=0,
        color=color1,
        alpha=alpha,
    )
    ax.set_ylim(*ylim)
    if xlabel:
        ax.set_xlabel("Maximum fitness decrease")
    if ylabel:
        ax.set_ylabel(r"Navigability, $\left<\psi\right>$")
    if color_labels:
        ax.yaxis.label.set_color(color1)
        ax.tick_params(axis="y", colors=color1)

    ax2.plot(
        navs["threshold"],
        navs["aborted_mean"],
        color=color2,
        ls=":",
        label=label2,
    )
    ax2.set_ylim(*ylim2)
    if ylabel2:
        ax2.set_ylabel(r"Aborted fraction, $\alpha$")
    if color_labels:
        ax2.yaxis.label.set_color(color2)
        ax2.tick_params(axis="y", colors=color2)

    return ax, ax2


def plot_nav_and_aborted_by_target(
    navs,
    ax,
    ax2=None,
    color1="dodgerblue",
    color2="firebrick",
    alpha=0.5,
    xticklabels=True,
    yticklabels=True,
    yticklabels2=True,
    xlabel=True,
    ylabel=True,
    ylabel2=True,
    ylim=(0, 50),
    ylim2=(0, 1),
    color_labels=True,
    label=None,
    label2=None,
):
    if ax2 is None:
        ax2 = ax.twinx()

    ax.set_ylim(*ylim)
    bins = np.arange(0, 1.01, 0.1)
    ax.hist(
        navs[0].dropna(),
        color=color1,
        zorder=1,
        alpha=alpha,
        label=label,
        bins=bins,
    )
    if xlabel:
        ax.set_xlabel(r"Navigability, $\left<\psi\right>$")
    if ylabel:
        ax.set_ylabel(r"Count")

    return ax, ax2
