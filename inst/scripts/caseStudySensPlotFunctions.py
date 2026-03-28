import numpy as np
import pandas as pd
import shutil
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D


import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# Enable LaTeX rendering
USE_TEX = shutil.which("latex") is not None
plt.rc("text", usetex=USE_TEX)
plt.rc("font", family="serif")

if USE_TEX:
    print("Matplotlib LaTeX rendering enabled.")
else:
    print("LaTeX not found; using Matplotlib's built-in mathtext.")

def buildColorMapping(scenarios):
    """Assign a unique color to each scenario using a colorblind-friendly palette."""

    palette = sns.color_palette("tab20", n_colors=max(3, len(scenarios)))
    return {scenario: palette[idx % len(palette)] for idx, scenario in enumerate(scenarios)}


def createCaseStudySensPlot(
    summaryDf,
    outputFile,
    width = 16.0,
    height = 12.0,
    dpi = 150,
):
    """Create the sensitivity plot summarizing multiplicative RA changes.

    Parameters
    ----------
    summary_df:
        A tidy dataframe with one row per time point that contains the
        columns ``scenario``, ``taxa``, ``exposure``, ``time``, ``mean``,
        ``lower`` and ``upper``.
    taxa:
        Ordered collection of taxa (e.g., ["Bacilli", ...]) that define the rows.
    exposures:
        Ordered collection of exposure labels (e.g., breast milk contrasts).
    scenarios:
        Ordered scenario names used for color assignment and the legend.
    outputFile:
        Optional path to save the resulting figure. When omitted, the Matplotlib
        figure object is returned instead.
    width, height:
        Figure size in inches.
    dpi:
        Resolution for the output image when ``outputFile`` is provided.
    """

    if summaryDf is None:
        raise ValueError("summary_df cannot be None")

    if not isinstance(summaryDf, pd.DataFrame):
        summaryDf = pd.DataFrame(summaryDf)

    taxaList = sorted(summaryDf["taxa"].unique())
    exposureList = list(summaryDf["exposure"].unique())
    scenarioList = list(summaryDf["scenario"].unique())
    labelsList = list(summaryDf["label"].unique())

    plotOrder = scenarioList + [scenarioList[0]] # plot default again to have it on top
    legendOrder = labelsList

    summaryDf = summaryDf.copy()
    summaryDf["time"] = summaryDf["time"].astype(float)
    summaryDf["mean"] = summaryDf["mean"].astype(float)
    summaryDf["lower"] = summaryDf["lower"].astype(float)
    summaryDf["upper"] = summaryDf["upper"].astype(float)

    colorMap = buildColorMapping(scenarioList)

    fig, axes = plt.subplots(
        len(taxaList),
        len(exposureList),
        figsize=(width, height),
        sharex=True,
        sharey=True,
        squeeze=False,
        constrained_layout=True,
    )

    fontsize = 10

    for rowIdx, taxon in enumerate(taxaList):
        for colIdx, exposure in enumerate(exposureList):
            ax = axes[rowIdx, colIdx]
            subset = summaryDf[
                (summaryDf["taxa"] == taxon)
                & (summaryDf["exposure"] == exposure)
            ]
            if subset.empty:
                raise ValueError(
                    f"No rows available for taxa '{taxon}' and exposure '{exposure}'."
                )

            subset = subset.sort_values("time")
            for scenario in plotOrder:
                scenarioSlice = subset[subset["scenario"] == scenario]
                if scenarioSlice.empty:
                    continue
                color = colorMap[scenario]
                sns.lineplot(
                    data=scenarioSlice,
                    x="time",
                    y="mean",
                    ax=ax,
                    color=color,
                    linewidth=1.5,
                    legend=False,
                )
                sns.lineplot(
                    data=scenarioSlice,
                    x="time",
                    y="lower",
                    ax=ax,
                    color=color,
                    linewidth=0.9,
                    linestyle="--",
                    alpha=0.75,
                    legend=False,
                )
                sns.lineplot(
                    data=scenarioSlice,
                    x="time",
                    y="upper",
                    ax=ax,
                    color=color,
                    linewidth=0.9,
                    linestyle="--",
                    alpha=0.75,
                    legend=False,
                )

            ax.axhline(1.0, color="black", linestyle=":", linewidth=0.9)
            ax.set_title(f"{taxon}: {exposure}", fontsize=11)
            if rowIdx == len(taxaList) - 1:
                ax.set_xlabel("Days after birth", fontsize=fontsize)
            if colIdx == 0:
                ax.set_ylabel("Multiplicative change in RA", fontsize=fontsize)
            ax.tick_params(axis='x', labelsize=fontsize, rotation=0)
            ax.tick_params(axis='y', labelsize=fontsize, rotation=0)

    legendHandles = [
        Line2D([0], [0], color=colorMap[scenario], label=scenario, linewidth=1.5)
        for scenario in legendOrder
    ]
    fig.legend(
        handles=legendHandles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=min(len(legendOrder), 4),
        fontsize=10,
        frameon=True,
    )
    fig.suptitle(
        "Sensitivity of Multiplicative Change in Relative Abundance",
        fontsize=14,
        y=0.99,
    )
    plt.tight_layout(rect=[0, 0.05, 1, 0.95])

    if outputFile is None:
        return fig

    path = outputFile
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path
