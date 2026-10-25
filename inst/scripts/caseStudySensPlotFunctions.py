import numpy as np
import pandas as pd
import shutil
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import os

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# Enable LaTeX rendering
os.environ["PATH"] = "/curc/sw/texlive/2021/bin/x86_64-linux" + ":" + os.environ["PATH"]
USE_TEX = shutil.which("latex") is not None
plt.rc("text", usetex=USE_TEX)
plt.rc("font", family="serif")

if USE_TEX:
    print("Matplotlib LaTeX rendering enabled.")
else:
    print("LaTeX not found; using Matplotlib's built-in mathtext.")

def buildColorMapping(scenarios):
    """Assign a unique color to each scenario using a colorblind-friendly palette."""

    palette = sns.color_palette("tab20", n_colors=max(3, len(scenarios) - 1))
    color_mapping = {}

    # First scenario is black
    color_mapping[scenarios[0]] = (0, 0, 0)
    
    # Remaining scenarios chosen from palette colors
    for idx, scenario in enumerate(scenarios[1:]):
        color_mapping[scenario] = palette[idx % len(palette)]

    return color_mapping


def createCaseStudySensPlot(
    summaryDf,
    outputFile = None,
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
            #ax.set_facecolor("white")
            #ax.grid(False)
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
                    alpha=0.5 if scenario != "CS_Orig" else 1,
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
                    alpha=0.25 if scenario != "CS_Orig" else .5,
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
                    alpha=0.25 if scenario != "CS_Orig" else .5,
                    legend=False,
                )

            ax.axhline(1.0, color="black", linestyle=":", linewidth=0.9)
            expTitle = exposure
            if "> 50" in exposure:
                expTitle = "Breast milk $> 50\\%$ vs $< 10\\%$"
            elif "10-50" in exposure:
                expTitle = "Breast milk $10-50\\%$ vs $< 10\\%$"
            ax.set_title(f"{taxon}: {expTitle}", fontsize=11)
            ax.tick_params(axis='both', which='both', labelsize=fontsize)
            ax.xaxis.set_tick_params(labelbottom=True)
            ax.yaxis.set_tick_params(labelleft=True)
            if rowIdx == len(taxaList) - 1:
                ax.set_xlabel("Days after birth", fontsize=fontsize)
            if colIdx == 0:
                ax.set_ylabel("$\Delta_v$RA$_{cp}[t, \mathbf{x}(t)]$", fontsize=fontsize)

    legendHandles = [
        Line2D([0], [0], color=colorMap[plotOrder[s]], label=legendOrder[s], linewidth=1.5) for s in range(len(legendOrder))
    ]
    fig.legend(
        handles=legendHandles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.03),
        ncol=min(len(legendOrder), 4),
        fontsize=10,
        frameon=True,
    )
    fig.suptitle(
        "Sensitivity of Multiplicative Difference in Relative Abundance",
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