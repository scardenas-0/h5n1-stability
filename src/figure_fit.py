#!/usr/bin/env python3

import argparse
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from grizzlyplot.scales import ScaleXCategorical
from matplotlib.ticker import ScalarFormatter

import analyze as ana
import plotting as plot


def main(
    data_path: str,
    titer_mcmc_path: str,
    halflife_mcmc_path: str,
    output_path: str,
    separator: str = "\t",
    prior_annotate: bool = True,
) -> None:
    """
    Create the main text display figure,
    with one panel showing inferred halflives
    and another showing model fits to the raw data

    Parameters
    ----------
    data_path : str
        Path to the data used to fit the model,
        as a delimited text file
        (default .tsv: tab-delimited, change this
        with the separator argument)

    titer_mcmc_path : str
        Path to the MCMC output for individual
        titer inference, saved as a .pickle archive.

    halflife_mcmc_path : str
        Path to the MCMC output for virus half-life
        inference, saved as a .pickle archive.

    output_path : str
        Path to which to save the figure.

    separator : str
        Delimiter for the delimited text
        file specified in data_path. Default
        `\t` (tab-delimited).

    prior_annotate : bool
       Annotate the plot with key prior values?
       Boolean, default True.
    """
    print(
        f"Creating figure {os.path.basename(output_path)}..."
    )

    hl_model = ana.load_mcmc(halflife_mcmc_path)[
        0
    ].run_model

    tidy_results = ana.get_tidy_results(
        data_path,
        titer_mcmc_path,
        halflife_mcmc_path,
    )
    titers = tidy_results["titers"].with_columns(
        display_titer=pl.when(pl.col("detected"))
        .then(10 ** pl.col("log_titer"))
        .otherwise(10 ** pl.col("log10_approx_lod"))
    ).with_columns(
        pl.col('sample_id').str.extract(r'rep(\d+)', 1).cast(pl.Int64).alias('rep_number')
    )
    
    hls = tidy_results["halflives"].with_columns(
        halflife_days=pl.col("halflife")
    )
    hls_int = tidy_results["halflives_with_intercepts"]

    hls_reg = ana.downsample_draws(
        hls_int, 10, id_column="sample_id"
    ).with_columns(
        initial_titer=10 ** pl.col("log_titer_intercept")
    )
    
    raw_milk(titers, hls, hls_reg, hl_model, output_path)
    surfaces(titers, hls, hls_reg, hl_model, output_path)
    water(titers, hls, hls_reg, hl_model, output_path)
    
def raw_milk(titers, 
             hls,
             hls_reg, 
             hl_model, 
             output_path):
    
    hls = hls.filter(medium_name="milk")
    titers = titers.filter(medium_name="milk")
    hls_reg = hls_reg.filter(medium_name="milk")
    
    reg_plot = plot.titer_regression(
        titers,
        hls_reg,
        facet={
            "col": "temperature_celsius",
            "sharex": False,
            "label_cols": False,
        },
    )

    hl_plot = plot.halflife_violins(
        hls,
        x_column="condition_id",
        halflife_column="halflife_days",
        additional_mappings=dict(
            fillcolor="condition_id",
            markerfacecolor="condition_id",
        ),
        scales=dict(
            fillcolor=plot.condition_color_scale,
            markerfacecolor=plot.condition_color_scale,
            x=ScaleXCategorical(),
        ),
        facet=dict(
            col="temperature_celsius",
            sharex=False,
            label_cols=False,
        ),
        markeredgewidth=3,
    )

    fig, ax = plt.subplots(
        2, 2, figsize=[10, 8], sharex=None, sharey="row"
    )

    reg_plot.render(fig=fig, ax=ax[0, ::])
    hl_plot.render(fig=fig, ax=ax[1, ::])
    fig.supxlabel(None)
    fig.supylabel(None)

    title_milk_4C = "4C"
    title_milk_22C = "22C"

    ax[0, 0].set_title(title_milk_4C)
    ax[0, 1].set_title(title_milk_22C)
    ax[0, 0].set_ylim([1e-1, 1e8])
    ax[1, 0].set_ylim([0, 4.5])
    ax[0, 0].set_xlabel("Time (days)", x=1)
    ax[0, 0].set_ylabel("Virus titer (TCID$_{50}$/mL)")

    ax[1, 0].yaxis.set_major_formatter(ScalarFormatter())
    ax[1, 0].set_ylabel("Half-life (days)")

    if prior_annotate:
        ax[1, 0].set_xlabel(
            plot.get_annotation_string(hl_model)
        )

    ax[1, 0].set_xticks((0,), labels = ["4C"])
    ax[1, 1].set_xticks((1,), labels = ["22C"])
    
    output_path = output_path + "-milk.pdf"
    
    print(f"Saving figure to {output_path}...")
    fig.savefig(output_path)

def surfaces(titers, hls, hls_reg, hl_model, output_path):
    titers = titers.filter((pl.col("medium_name")=="steel") | (pl.col("medium_name")=="polypropylen"))
    hls_reg = hls_reg.filter((pl.col("medium_name")=="steel") | (pl.col("medium_name")=="polypropylen"))
    hls = hls.filter((pl.col("medium_name")=="steel") | (pl.col("medium_name")=="polypropylen"))
    
    reg_plot = plot.titer_regression(
        titers,
        hls_reg,
        facet={
            "row": "medium_name",
            "col": "temperature_celsius",
            "sharex": False,
            "label_cols": False,

            "label_rows": False
        },
    )
    
    hl_plot = plot.halflife_violins(hls,
        x_column= "medium_name",
        halflife_column="halflife_days",
        additional_mappings=dict(
            fillcolor="condition_id",
            markerfacecolor="condition_id",
        ),
        scales=dict(
            fillcolor=plot.condition_color_scale,
            markerfacecolor=plot.condition_color_scale,
            x=ScaleXCategorical(),
        ),
        facet=dict(
            col="temperature_celsius",
            sharex=False,
            label_cols=False,
        ),
        markeredgewidth=3,
    )

    fig, ax = plt.subplots(
        3, 2, figsize=[10, 12], sharex=None, sharey='row'
    )

    reg_plot.render(fig=fig, ax=ax[:2,::]) 
    hl_plot.render(fig=fig, ax=ax[2, ::])
    fig.supxlabel(None)
    fig.supylabel(None)

    title_4C = "4C"
    title_22C = "22C"

    ax[0, 0].set_title(title_4C)
    ax[0, 1].set_title(title_22C)
    ax[0, 0].set_ylim([1e-1, 1e8])
    ax[1, 0].set_ylim([1e-1, 1e8])

    # ax[2, 0].set_yscale('log')
    ax[2, 0].set_ylim([0, 3])
    ax[2, 0].set_xticks((0, 1), labels = ["Polypropylene", "Steel"])
    ax[2, 1].set_xticks((0, 1), labels = ["Polypropylene", "Steel"])
    
    ax[1, 0].set_xlabel("Time (days)", x=1)
    ax[0, 0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    ax[1, 0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    ax0_secondary = ax[0,1].twinx()
    ax1_secondary = ax[1,1].twinx()
    ax0_secondary.set_ylabel("Polypropylene", rotation=270, labelpad=15)
    ax1_secondary.set_ylabel("Steel", rotation=270, labelpad=15)
    ax0_secondary.set_yticklabels("")
    ax1_secondary.set_yticklabels("")
    ax0_secondary.grid(visible=False)
    ax1_secondary.grid(visible=False)
    
    ax[0, 0].set_xticklabels("")
    ax[0, 1].set_xticklabels("")

    ax[2, 0].yaxis.set_major_formatter(ScalarFormatter())
    ax[2, 0].set_ylabel("Half-life (days)")


    ax[2, 0].grid(visible = True, which = 'major', axis = 'both')
    ax[2, 1].grid(visible = True, which = 'major', axis = 'both')
    
    if prior_annotate:
        ax[2, 0].set_xlabel(
            plot.get_annotation_string(hl_model)
        )

    output_path = output_path + "-surfaces.pdf"
    print(f"Saving figure to {output_path}...")
    fig.savefig(output_path)

def surfaces2(titers, hls, hls_reg, hl_model, output_path):
    titers = titers.filter((pl.col("medium_name")=="steel") | (pl.col("medium_name")=="polypropylen"))
    hls_reg = hls_reg.filter((pl.col("medium_name")=="steel") | (pl.col("medium_name")=="polypropylen"))
    hls = hls.filter((pl.col("medium_name")=="steel") | (pl.col("medium_name")=="polypropylen"))
    
    reg_plot = plot.titer_regression(
        titers,
        hls_reg,
        facet={
            "row": "medium_name",
            # "col": "temperature_celsius",
            "sharex": False,
            "label_cols": False,
            # "color": "temperature_celsius" #?
        },
    )
    
    hl_plot = plot.halflife_violins(hls,
        x_column= "medium_name",
        halflife_column="halflife_days",
        additional_mappings=dict(
            fillcolor="condition_id",
            markerfacecolor="condition_id",
        ),
        scales=dict(
            fillcolor=plot.condition_color_scale,
            markerfacecolor=plot.condition_color_scale,
            x=ScaleXCategorical(),
        ),
        facet=dict(
            # col="temperature_celsius",
            sharex=False,
            label_cols=False,
            color = "temperature_celsius"
        ),
        markeredgewidth=3,
    )

    fig, ax = plt.subplots(
        3, 1, figsize=[10, 12], sharex=None, sharey='row'
    )

    reg_plot.render(fig=fig, ax=ax[:2]) 
    hl_plot.render(fig=fig, ax=ax[2])
    fig.supxlabel(None)
    fig.supylabel(None)

    title_surfaces = "Surfaces"

    ax[0].set_title(title_surfaces)
    ax[0].set_ylim([1e-1, 1e8])
    ax[1].set_ylim([1e-1, 1e8])

    # ax[2].set_yscale('log')
    ax[2].set_ylim([0, 3])
    ax[2].set_xticks((0, 1), labels = ["polypropylene", "steel"])
    # ax[2, 1].set_xticks((0, 1), labels = ["polypropylene", "steel"])
    
    ax[1].set_xlabel("Time (days)", x=0.5)
    ax[0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    ax[1].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    ax[0].set_xticklabels("")

    ax[2].yaxis.set_major_formatter(ScalarFormatter())
    ax[2].set_ylabel("Half-life (days)")


    # ax[2].set_yticks((0.1, 0.3, 1, 3))
    ax[2].grid(visible = True, which = 'major', axis = 'both')
    
    if prior_annotate:
        ax[2].set_xlabel(
            plot.get_annotation_string(hl_model)
        )

    output_path = output_path + "-surfaces-2.pdf"
    print(f"Saving figure to {output_path}...")
    fig.savefig(output_path)

def water(titers, 
             hls,
             hls_reg, 
             hl_model, 
             output_path):
    

    titers = titers.filter(pl.col("medium_name")=="wastewater")
    hls_reg = hls_reg.filter(pl.col("medium_name")=="wastewater")
    hls = hls.filter(pl.col("medium_name")=="wastewater")
    
    reg_plot = plot.titer_regression(
        titers,
        hls_reg,
        facet={
            "col": "medium_name",
            "sharex": False,
            "label_cols": False,
        },
    )

    hl_plot = plot.halflife_violins(
        hls,
        x_column="condition_id",
        halflife_column="halflife_days",
        additional_mappings=dict(
            fillcolor="condition_id",
            markerfacecolor="condition_id",
        ),
        scales=dict(
            fillcolor=plot.condition_color_scale,
            markerfacecolor=plot.condition_color_scale,
            x=ScaleXCategorical(),
        ),
        facet=dict(
            col="medium_name",
            sharex=False,
            label_cols=False,
        ),
        markeredgewidth=3,
    )

    fig, ax = plt.subplots(

        2, 1, figsize=[8, 8], sharex=None, sharey="row"
    )

    reg_plot.render(fig=fig, ax=ax[0])
    hl_plot.render(fig=fig, ax=ax[1])
    fig.supxlabel(None)
    fig.supylabel(None)

    title_milk_wwater = "Wastewater"

    # title_milk_DI = "DI water"

    ax[0].set_title(title_milk_wwater)
    # ax[0, 0].set_title(title_milk_DI)
    ax[0].set_ylim([1e-2, 1e5])
    ax[1].set_ylim([0, 0.8])
    ax[0].set_xlabel("Time (days)", x=0.5)
    ax[0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    ax[1].yaxis.set_major_formatter(ScalarFormatter())
    ax[1].set_ylabel("Half-life (days)")

    if prior_annotate:
        ax[1].set_xlabel(
            plot.get_annotation_string(hl_model)
        )

    # ax[1,0].set_xticks((0,), labels = ["DI"])
    ax[1].set_xticks((0,), labels = ["Wastewater"])
    
    output_path = output_path + "-water.pdf"
    print(f"Saving figure to {output_path}...")
    fig.savefig(output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Read in MCMC output and produce a main "
            "text figure showing fits to the data "
            "and inferred halflives"
        )
    )
    parser.add_argument(
        "data_path",
        type=str,
        help=(
            "Path to the data used for fitting, formatted as "
            "a delimited text file"
        ),
    )
    parser.add_argument(
        "titer_mcmc_path",
        type=str,
        help=(
            "Path to the MCMC output for individual "
            "titer inference, saved as a .pickle archive."
        ),
    )
    parser.add_argument(
        "halflife_mcmc_path",
        type=str,
        help=(
            "Path to the MCMC output for virus "
            "half-life inference, saved as a .pickle archive."
        ),
    )
    parser.add_argument(
        "output_path",
        type=str,
        help=("Path to save the generated figure."),
    )
    parser.add_argument(
        "-s",
        "--separator",
        type=str,
        help=(
            "Separator for the delimited text file containing "
            "the data (specified in data_path)"
        ),
        default="\t",
    )
    parsed = vars(parser.parse_args())
    # set seed for reproducibility
    # (since we use random draws)
    np.random.seed(52367)
    # do not annotate main text figure
    prior_annotate = "default" not in parsed["output_path"]
    main(
        parsed["data_path"],
        parsed["titer_mcmc_path"],
        parsed["halflife_mcmc_path"],
        parsed["output_path"],
        separator=parsed["separator"],
        prior_annotate=prior_annotate,
    )
