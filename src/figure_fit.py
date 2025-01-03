#!/usr/bin/env python3

import argparse
import os
import re

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

    old_titers_liquid_plot(titers, hls, hls_reg, output_path)
    old_titers_surface_plot(titers, hls, hls_reg, output_path)
    old_hl_plot(hls, output_path)
    
    bulk_milk_plot(titers, hls, hls_reg, output_path)
    surface_plot(titers, hls, hls_reg, output_path)
    hl_plot(hls, output_path)

def old_titers_liquid_plot(
        titers,
        hls,
        hls_reg,
        output_path
    ):
    media_names = ["raw", "wastewater"]
    
    titers = (
        titers.filter(pl.col("medium_name").is_in(media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(media_names)))
            .sort("medium_name")
    )

    hls_reg = (
        hls_reg.filter(pl.col("medium_name").is_in(media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(media_names)))
            .sort("medium_name")
    )

    hls = (
        hls.filter(pl.col("medium_name").is_in(media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(media_names)))
            .sort("medium_name")
    )
    
    reg_plot = plot.titer_regression(
        titers,
        hls_reg,
        facet = {
            "col": "medium_name",
            "row": "temperature_celsius",
            "sharex": True,
            "sharey": True,
            "label_cols": False,
            "label_rows": False,
        },
    )

    fig, ax = plt.subplots(
        2, 2, figsize=[11, 8], sharex='all', sharey='all'
    )

    reg_plot.render(fig=fig, ax=ax) 
    fig.supxlabel(None)
    fig.supylabel(None)

    
    ax[0,0].set_title("Raw milk")
    ax[1,1].set_title("Wastewater")

    ax[0,0].set_ylim([1e-1, 1e8])
    
    ax[1,0].set_xlabel("Time (days)", x=0.5)
    ax[1,1].set_xlabel("Time (days)", x=0.5)

    ax[0,0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    ax[1,0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    legend_elements = [Line2D([0], [0], color='orange', lw=8, label='22C'),
                   Line2D([0], [0], color='blue', lw = 8, label='4C')]

    ax[0,1].legend(handles=legend_elements, loc='center', prop = {"size": 20})
    ax[0,1].grid(visible = False, which = 'major', axis = 'both')

    output_path = output_path + "-titers-liquid.pdf"
    print(f"Saving figure to {output_path}...")
    fig.savefig(output_path)

def old_titers_surface_plot(
        titers,
        hls,
        hls_reg,
        output_path
    ):
    media_names = ["polypropylen", "steel", "rubber"]
    
    titers = (
        titers.filter(pl.col("medium_name").is_in(media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(media_names)))
            .sort("medium_name")
    )

    hls_reg = (
        hls_reg.filter(pl.col("medium_name").is_in(media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(media_names)))
            .sort("medium_name")
    )

    hls = (
        hls.filter(pl.col("medium_name").is_in(media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(media_names)))
            .sort("medium_name")
    )
    
    reg_plot = plot.titer_regression(
        titers,
        hls_reg,
        facet = {
            "col": "medium_name",
            "row": "temperature_celsius",
            "sharex": True,
            "sharey": True,
            "label_cols": False,
            "label_rows": False,
        },
    )

    fig, ax = plt.subplots(
        2, 3, figsize=[16, 8], sharex='all', sharey='all'
    )

    reg_plot.render(fig=fig, ax=ax) 
    fig.supxlabel(None)
    fig.supylabel(None)

    ax[0,0].set_title("Polypropylene plastic")
    ax[0,1].set_title("Steel")
    ax[0,2].set_title("Rubber")

    ax[0,0].set_ylim([1e-1, 1e8])
    
    ax[1,0].set_xlabel("Time (days)", x=0.5)
    ax[1,1].set_xlabel("Time (days)", x=0.5)
    ax[1,2].set_xlabel("Time (days)", x=0.5)

    ax[0,0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    ax[1,0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    legend_elements = [Line2D([0], [0], color='orange', lw = 4, label='22C'),
                   Line2D([0], [0], color='blue', lw = 4, label='4C')]

    ax[0,2].legend(handles=legend_elements, loc='upper right', prop = {"size": 12})
    # ax[0,2].grid(visible = False, which = 'major', axis = 'both')

    output_path = output_path + "-titers-surface.pdf"
    print(f"Saving figure to {output_path}...")
    fig.savefig(output_path)

def old_hl_plot(
        hls,
        output_path
    ):
    """
    Generates and saves a violin plot of half-lives for the different conditions.
    Parameters:
    hls (DataFrame): DataFrame containing half-life data.
    hl_model (Model): Model used for half-life calculations.
    output_path (str): Path to save the generated plot.
    Returns:
    None
    """
    half_lives_ordering = ["polypropylen", "steel", "rubber", "raw", "wastewater"]
    
    hls = (
        hls.filter(pl.col("medium_name").is_in(half_lives_ordering))
            .with_columns(pl.col("medium_name").cast(pl.Enum(half_lives_ordering)))
            .sort("medium_name")
    )
    
    hl_violins = plot.halflife_violins(
        hls,
        x_column="medium_name",
        halflife_column="halflife_days",
        additional_mappings=dict(
            fillcolor="temperature_celsius",
            markerfacecolor="temperature_celsius",
        ),
        scales=dict(
            fillcolor=plot.condition_color_scale,
            markerfacecolor=plot.condition_color_scale,
            x=ScaleXCategorical(),
        ),
        facet=dict(
            sharex=False,
            label_cols=False,
            color="temperature_celsius"
        ),
        markeredgewidth=3,
    )

    fig, ax = plt.subplots(
        1, 1, figsize = [12, 6], 
        sharex = None, sharey = None
        )
    
    hl_violins.render(fig=fig, ax=ax)
    fig.supxlabel(None)
    fig.supylabel(None)

    title_hls = "Half-lives"
    ax.set_title(title_hls)
    
    ax.set_ylim([0, 4.5])
    ax.set_xlabel("Medium", x=0.5)

    ax.set_xticks(
        ticks = [0, 1, 2, 3, 4],
        labels = ["Polypropylene", "Steel", "Rubber", "Raw milk", "Wastewater"]
    )

    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.set_ylabel("Half-life (days)")

    ax.grid(visible = True, which = 'major', axis = 'both')
    legend_elements = [Line2D([0], [0], color='orange', lw=8, label='22C'),
                   Line2D([0], [0], color='blue', lw = 8, label='4C')]

    ax.legend(handles=legend_elements, loc='upper right', prop = {"size": 15})
    
    output_path = output_path + "-original-halflives.pdf"
    print(f"Saving figure to {output_path}...")
    fig.savefig(output_path)

def bulk_milk_plot(
        titers,
        hls,
        hls_reg,
        output_path
    ):
    milk_media_names = ["raw", "skim", "fat", "whole", "wastewater"]
    
    titers = (
        titers.filter(pl.col("medium_name").is_in(milk_media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(milk_media_names)))
            .sort("medium_name")
    )

    hls_reg = (
        hls_reg.filter(pl.col("medium_name").is_in(milk_media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(milk_media_names)))
            .sort("medium_name")
    )

    hls = (
        hls.filter(pl.col("medium_name").is_in(milk_media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(milk_media_names)))
            .sort("medium_name")
    )
    
    reg_plot = plot.titer_regression(
        titers,
        hls_reg,
        facet = {
            "col": "medium_name",
            "row": "temperature_celsius",
            "sharex": True,
            "sharey": True,
            "label_cols": False,
            "label_rows": False,
        },
    )

    fig, ax = plt.subplots(
        2, 5, figsize=[16, 8], sharex='all', sharey='all'
    )

    reg_plot.render(fig=fig, ax=ax) 
    fig.supxlabel(None)
    fig.supylabel(None)

    ax[0,0].set_title("Raw")
    ax[0,1].set_title("Skim")
    ax[0,2].set_title("2%")
    ax[0,3].set_title("Whole")
    ax[1,4].set_title("Wastewater")

    ax[0,0].set_ylim([1e-1, 1e8])
    
    ax[1,0].set_xlabel("Time (days)", x=0.5)
    ax[1,1].set_xlabel("Time (days)", x=0.5)
    ax[1,2].set_xlabel("Time (days)", x=0.5)
    ax[1,3].set_xlabel("Time (days)", x=0.5)
    ax[1,4].set_xlabel("Time (days)", x=0.5)

    ax[0,0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    ax[1,0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    legend_elements = [Line2D([0], [0], color='orange', lw=8, label='22C'),
                   Line2D([0], [0], color='blue', lw = 8, label='4C')]

    ax[0,4].legend(handles=legend_elements, loc='center', prop = {"size": 20})
    ax[0,4].grid(visible = False, which = 'major', axis = 'both')

    output_path = output_path + "-milk-titers.pdf"
    print(f"Saving figure to {output_path}...")
    fig.savefig(output_path)

def surface_plot(
        titers,
        hls,
        hls_reg,
        output_path
    ):
    surface_media_names = ["polypropylen", "steel", "rubber"]
    
    titers = (
        titers.filter(pl.col("medium_name").is_in(surface_media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(surface_media_names)))
            .sort("medium_name")
    )

    hls_reg = (
        hls_reg.filter(pl.col("medium_name").is_in(surface_media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(surface_media_names)))
            .sort("medium_name")
    )

    hls = (
        hls.filter(pl.col("medium_name").is_in(surface_media_names))
            .with_columns(pl.col("medium_name").cast(pl.Enum(surface_media_names)))
            .sort("medium_name")
    )
    
    reg_plot = plot.titer_regression(
        titers,
        hls_reg,
        facet = {
            "col": "medium_name",
            "row": "temperature_celsius",
            "sharex": True,
            "sharey": True,
            "label_cols": False,
            "label_rows": False,
        },
    )

    fig, ax = plt.subplots(
        2, 3, figsize=[16, 8], sharex='all', sharey='all'
    )

    reg_plot.render(fig=fig, ax=ax) 
    fig.supxlabel(None)
    fig.supylabel(None)

    ax[0,0].set_title("Polypropylene plastic")
    ax[0,1].set_title("Steel")
    ax[0,2].set_title("Rubber")

    ax[0,0].set_ylim([1e-1, 1e8])
    
    ax[1,0].set_xlabel("Time (days)", x=0.5)
    ax[1,1].set_xlabel("Time (days)", x=0.5)
    ax[1,2].set_xlabel("Time (days)", x=0.5)

    ax[0,0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    ax[1,0].set_ylabel("Virus titer (TCID$_{50}$/mL)")
    legend_elements = [Line2D([0], [0], color='orange', lw=8, label='22C'),
                   Line2D([0], [0], color='blue', lw = 8, label='4C')]

    ax[1,2].legend(handles=legend_elements, loc='upper right', prop = {"size": 15})

    output_path = output_path + "-surface-titers.pdf"
    print(f"Saving figure to {output_path}...")
    fig.savefig(output_path)


def hl_plot(
        hls,
        output_path
    ):
    """
    Generates and saves a violin plot of half-lives for the different conditions.
    Parameters:
    hls (DataFrame): DataFrame containing half-life data.
    hl_model (Model): Model used for half-life calculations.
    output_path (str): Path to save the generated plot.
    Returns:
    None
    """
    half_lives_ordering = ["polypropylen", "steel", "rubber", "raw", "skim", "fat", "whole", "wastewater"]
    
    hls = (
        hls.filter(pl.col("medium_name").is_in(half_lives_ordering))
            .with_columns(pl.col("medium_name").cast(pl.Enum(half_lives_ordering)))
            .sort("medium_name")
    )
    
    hl_violins = plot.halflife_violins(
        hls,
        x_column="medium_name",
        halflife_column="halflife_days",
        additional_mappings=dict(
            fillcolor="temperature_celsius",
            markerfacecolor="temperature_celsius",
        ),
        scales=dict(
            fillcolor=plot.condition_color_scale,
            markerfacecolor=plot.condition_color_scale,
            x=ScaleXCategorical(),
        ),
        facet=dict(
            sharex=False,
            label_cols=False,
            color="temperature_celsius"
        ),
        markeredgewidth=3,
    )

    fig, ax = plt.subplots(
        1, 1, figsize = [12, 6], 
        sharex = None, sharey = None
        )
    
    hl_violins.render(fig=fig, ax=ax)
    fig.supxlabel(None)
    fig.supylabel(None)

    title_hls = "Half-lives"
    ax.set_title(title_hls)
    
    ax.set_ylim([0, 4.5])
    ax.set_xlabel("Medium", x=0.5)

    ax.set_xticks(
        ticks = [0, 1, 2, 3, 4, 5, 6, 7],
        labels = ["PP", "Steel", "Rubber", "Raw milk", "Skim", "2%", "Whole", "Wastewater"]
    )

    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.set_ylabel("Half-life (days)")

    ax.grid(visible = True, which = 'major', axis = 'both')
    legend_elements = [Line2D([0], [0], color='orange', lw=8, label='22C'),
                   Line2D([0], [0], color='blue', lw = 8, label='4C')]

    ax.legend(handles=legend_elements, loc='upper right', prop = {"size": 15})
    
    output_path = output_path + "-halflives.pdf"
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
    parsed = vars(parser.parse_args())

    # set seed for reproducibility
    # (since we use random draws)
    np.random.seed(52367)
    main(
        parsed["data_path"],
        parsed["titer_mcmc_path"],
        parsed["halflife_mcmc_path"],
        parsed["output_path"],
    )
