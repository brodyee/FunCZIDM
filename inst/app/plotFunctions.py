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

def ciAxLabels(ax, title, ylabel, xlabel, fontsize):
    ax.set_title(title, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.tick_params(axis='x', labelsize=fontsize, rotation=0)
    ax.tick_params(axis='y', labelsize=fontsize, rotation=0)
    #ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=fontsize)
    #ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=fontsize)

def makeCIPlot(data, ax, cat, titlePrefix, ylabel, xLabel, isRA=True, 
               labelX=True, changeIn=True, by=0.5, fontsize=10, 
               lineWidth=1, yMin=None, yMax=None):
    sns.lineplot(data=data, x="time", y="lowerCI", color="gray", linestyle="--",
                 ax=ax, linewidth=lineWidth)
    sns.lineplot(data=data, x="time", y="upperCI", color="gray", linestyle="--",
                 ax=ax, linewidth=lineWidth)
    sns.lineplot(data=data, x="time", y="mean",   color="blue", ax=ax, 
                 linewidth=lineWidth)

    sameYaxis = True
    if yMin is None or yMax is None:
        yMin, yMax = min(data["lowerCI"]) - .03, max(data["upperCI"]) + .03
        sameYaxis = False
    if changeIn:
      noEffect = 1 if isRA else 0
      if (noEffect >= yMin) or (noEffect <= yMax):
        ax.axhline(y=noEffect, color="red", linestyle="-.", linewidth=lineWidth)
    
    if isRA and not changeIn and (yMin < 0):
        yMin = 0
    if isRA and not changeIn and (yMax > 1):
        yMax = 1
    if sameYaxis:
        y_ticks = np.linspace(yMin, yMax, by)
    else:
        y_ticks = np.arange(yMin, yMax, by)
    
    roundTo = 2
    while len(np.unique(y_ticks)) != len(np.unique(np.round(y_ticks, roundTo))):
        roundTo += 1
    y_ticks = np.round(y_ticks, roundTo)
    ax.set_yticks(y_ticks)
    
    xlabel = xLabel if labelX else ""
    full_title = f"{titlePrefix}{cat}"
    ciAxLabels(ax, full_title, ylabel, xlabel, fontsize)

def curvePlots(dataDict, suptitle, raTitlePrefix, betaTitlePrefix, raYlabel, 
               betaYlabel, xLabel, changeIn, plotBetas, lineWidth=1, byRA=0.05,
               byBeta=0.5, fontsize=10, titleSize=17, figsize=(15, 10), 
               xAxisRange=None, figRowCol=None, legendLoc=(0.46, -0.01), 
               suptitleLoc=(.45, .77), sameYaxis=False):
    n_cats = len(dataDict)
    if plotBetas:
        n_rows = n_cats
        n_cols = 2
    else:
        n_rows = figRowCol[0] if figRowCol[0] != -1 else \
                     n_cats//figRowCol[1] + 1 if n_cats%figRowCol[1] != 0 else \
                     n_cats//figRowCol[1]
        n_cols = figRowCol[1] if figRowCol[1] != -1 else \
                                 n_cats//n_rows + 1 if n_cats%n_rows != 0 else \
                                 n_cats//n_rows 
        n_rows, n_cols = (int(n_rows), int(n_cols)) 
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False, 
                             constrained_layout=True)
    axes = axes.flatten()
    
    if sameYaxis:
        yMinRA = min(min(subdict["RA"]["lowerCI"]) for subdict in
                     dataDict.values()) - .03
        yMaxRA = max(max(subdict["RA"]["upperCI"]) for subdict in 
                     dataDict.values()) + .03

        if plotBetas:
            yMinBeta = min(min(subdict["beta"]["lowerCI"]) for subdict in 
                           dataDict.values()) - .03
            yMaxBeta = max(max(subdict["beta"]["upperCI"]) for subdict in 
                           dataDict.values()) + .03
            byBeta = byRA
        else:
            yMinBeta, yMaxBeta = None, None
    else:
        yMinRA, yMaxRA, yMinBeta, yMaxBeta = None, None, None, None

    for idx, (cat, subdict) in enumerate(dataDict.items()):
        # Left column: RA‐change
        ax_ra = axes[2*idx] if plotBetas else axes[idx]
        df_ra = pd.DataFrame(subdict["RA"])
        df_ra.columns = ["lowerCI", "upperCI", "mean", "time"]
        if xAxisRange is not None:
            df_ra = df_ra[(df_ra["time"] >= xAxisRange[0]) 
                          & (df_ra["time"] <= xAxisRange[1])]
        labelX_ra = (idx == n_cats - 1) if plotBetas else \
                                                 (n_rows*n_cols - n_cols <= idx)
        ylabel = raYlabel if plotBetas else raYlabel if (idx%n_cols == 0) else \
                                                                              ""
        makeCIPlot(data=df_ra, ax=ax_ra, cat=cat, titlePrefix=raTitlePrefix, 
                   ylabel=ylabel, xLabel=xLabel, isRA=True, labelX=labelX_ra, 
                   changeIn=changeIn, by=byRA, fontsize=fontsize, 
                   lineWidth=lineWidth, yMin=yMinRA, yMax=yMaxRA)
        
        # Right column: beta
        if not plotBetas:
            continue
        ax_beta = axes[2*idx + 1]
        df_beta = pd.DataFrame(subdict["beta"])
        df_beta.columns = ["lowerCI", "upperCI", "mean", "time"]
        if xAxisRange is not None:
            df_beta = df_beta[(df_beta["time"] >= xAxisRange[0]) 
                              & (df_beta["time"] <= xAxisRange[1])]
        labelX_beta = (idx == n_cats - 1)
        makeCIPlot(data=df_beta, ax=ax_beta, cat=cat, 
                   titlePrefix=betaTitlePrefix, ylabel=betaYlabel, 
                   xLabel=xLabel, isRA=False, labelX=labelX_beta, 
                   changeIn=changeIn, by=byBeta, fontsize=fontsize, 
                   lineWidth=lineWidth, yMin=yMinBeta, yMax=yMaxBeta)

    # Shared legend at bottom center
    legend_handles = [Line2D([0], [0], color="gray", linestyle="--", 
                             label="Lower/Upper 95\% CI"),
                      Line2D([0], [0], color="blue", label="Mean"),
                      Line2D([0], [0], color="red", linestyle="-.", 
                             label="No Effect")]
    if not changeIn:
        legend_handles = legend_handles[:-1]  # Remove "No Effect" 
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
               bbox_to_anchor=legendLoc, fontsize=fontsize)

    fig.suptitle(suptitle, fontsize=titleSize, y=suptitleLoc[1],
                 x=suptitleLoc[0])
    plt.close(fig)
    
    return fig

def divCurvePlots(dataDict, suptitle, yLabel, xLabel, changeIn, 
                  multiChange = False, lineWidth = 1, byY = 0.05, fontsize = 10,
                  titleSize = 17, figsize = (15, 10), xAxisRange = None, 
                  legendLoc = (0.46, -0.01), suptitleLoc = (.45, .77)):
    fig, axes = plt.subplots(1, 1, figsize=figsize, squeeze=False, 
                             constrained_layout=True)
    axes = axes.flatten()
    ax = axes[0]
    df = pd.DataFrame(dataDict["div"])
    df.columns = ["lowerCI", "upperCI", "mean", "time"]
    if xAxisRange is not None:
        df = df[(df["time"] >= xAxisRange[0]) & (df["time"] <= xAxisRange[1])]

    makeCIPlot(data=df, ax=ax, cat="", titlePrefix="", ylabel=yLabel, 
               xLabel=xLabel, isRA=multiChange, labelX=True, changeIn=changeIn,
               by=byY, fontsize=fontsize, lineWidth=lineWidth)

    # Shared legend at bottom center
    legend_handles = [
        Line2D([0], [0], color="gray", linestyle="--", 
               label="Lower/Upper 95\% CI"),
        Line2D([0], [0], color="blue", label="Mean"),
        Line2D([0], [0], color="red", linestyle="-.", label="No Effect")
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=legendLoc,
        fontsize=fontsize
    )
    
    fig.suptitle(suptitle, fontsize=titleSize, y=suptitleLoc[1], 
                 x=suptitleLoc[0])
    plt.close(fig)
    
    return fig

def heatmapAxLabels(ax, cat, fontsize, xlabel, ylabel, xPos, yPos, xLabels, 
                    yLabels):
    ax.set_title(cat, fontsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_xticks(xPos)
    ax.set_yticks(yPos)
    ax.set_yticklabels(yLabels, rotation=0, fontsize=fontsize)
    ax.set_xticklabels(xLabels, rotation=0, fontsize=fontsize)

def makeHeatmapPlot(ax, heatmapDat, norm, cmap, cat, fontsize, xlabel, ylabel, 
                    cbar = False, 
                    cbarLabel = r"$\Delta_{v} \mathrm{RA}_{cp}(t)$", 
                    cbarTicks = None, xPos = None, yPos = None, xLabels = None, 
                    yLabels = None):
    # 3. Draw the seaborn heatmap
    hm = sns.heatmap(heatmapDat, cmap=cmap, norm=norm, ax=ax, cbar=cbar,
                     cbar_kws={"label": cbarLabel, "ticks": cbarTicks} if cbar 
                                                                       else {})
    if cbar:
        cbar_obj = hm.collections[0].colorbar   # get the matplotlib colorbar
        cbar_obj.ax.tick_params(labelsize=fontsize)  # tick labels
        cbar_obj.ax.yaxis.label.set_size(fontsize)   # cbar label
        
    # 4. Apply titles and axis labels via helper
    heatmapAxLabels(ax, cat, fontsize, xlabel, ylabel, xPos, yPos, xLabels,
                    yLabels)

    return hm

def heatmapPlots(dataDict, suptitle, covChange, inCI, figsize = (15, 6), 
                 fontsize = 12, xlabel = "", ylabel = "", cmap = 'PuOr', 
                 vmin = 0, vcenter = 1, vmax = 2, 
                 cbarLabel = r"$\Delta_{v} \mathrm{RA}_{cp}(t)$", 
                 suptitleLoc = (1, 1.025), titleSize = 17, numXticks = 10,
                 numYticks = 10):
    vmax = max(np.max(df["means"]) for df in dataDict.values())
    vmin = min(np.min(df["means"]) for df in dataDict.values())
    maxDist = max(vmax - vcenter, vcenter - vmin)
    if maxDist > 1:
        norm = TwoSlopeNorm(vmin=0, vcenter=vcenter, vmax=vmax)
        cbarTicks = np.linspace(0, vmax, 10)
    else:
        norm = TwoSlopeNorm(vmin=vcenter - maxDist, vcenter=vcenter,
                            vmax=vcenter + maxDist)
        cbarTicks = np.linspace(vcenter - maxDist, vcenter + maxDist, 10)
    roundTo = 2
    while len(np.unique(cbarTicks)) != \
                                   len(np.unique(np.round(cbarTicks, roundTo))):
        roundTo += 1
    cbarTicks = np.round(cbarTicks, roundTo)
    
    # 2. Create a row of subplots, one per category
    n_cats = len(dataDict)
    fig, axes = plt.subplots(1, n_cats, figsize=figsize, squeeze=False, 
                             constrained_layout=True)
    axes = axes.flatten()
    for ax, (cat, dataList) in zip(axes, dataDict.items()):
        heatmapDat = pd.DataFrame(dataList["means"])
        heatmapDat.columns = np.round(dataList["testpoints"], 2)
        heatmapDat.index = np.round(covChange, 2)
        if inCI:
            heatmapDat = heatmapDat.mul(1 - dataList["nullInCI"])
            heatmapDat = heatmapDat.replace(0, 1)
        # X-axis ticks
        xPos = np.linspace(0, heatmapDat.shape[1]-1, numXticks, dtype=int)
        xLabels = [heatmapDat.columns[i] for i in xPos]

        # Y-axis ticks
        yPos = np.linspace(0, heatmapDat.shape[0]-1, numYticks, dtype=int)
        yLabels = [heatmapDat.index[i] for i in yPos]

        hm = makeHeatmapPlot(ax=ax, heatmapDat=heatmapDat, norm=norm, cmap=cmap,
                             cat=cat, fontsize=fontsize, xlabel=xlabel, 
                             ylabel="", xPos=xPos, yPos=yPos, xLabels=xLabels,
                             yLabels=yLabels)
    
    # 4. Adjust layout, add suptitle, save
    pos = axes[0].get_position()
    fig.supylabel(ylabel, fontsize=fontsize, x=-0.025)
    pos = axes[-1].get_position()
    cbar_ax = fig.add_axes([pos.x1 + 0.01, pos.y0, 0.02, pos.height])
    cbar = fig.colorbar(hm.collections[0], cax=cbar_ax, ticks=cbarTicks)
    cbar.set_label(cbarLabel, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    for spine in cbar_ax.spines.values():
        spine.set_visible(False)
    fig.suptitle(suptitle, fontsize=titleSize, y=suptitleLoc[1],
                 x=suptitleLoc[0])
    plt.close(fig)

    return fig

def divHeatmapPlots(dataDict, suptitle, covChange, inCI, figsize=(15, 6), 
                    fontsize = 12, xlabel = "", ylabel = "", cmap = 'PuOr', 
                    vmin = 0, vcenter = 1, vmax = 2, 
                    cbarLabel = r"$\Delta_{v} \mathrm{RA}_{cp}(t)$", 
                    suptitleLoc = (1, 1.025), titleSize = 17, numXticks = 10, 
                    numYticks = 10):
    vmax = max(np.max(df["means"]) for df in dataDict.values())
    vmin = min(np.min(df["means"]) for df in dataDict.values())
    maxDist = max(vmax - vcenter, vcenter - vmin)
    if maxDist > 1:
        norm = TwoSlopeNorm(vmin=0, vcenter=vcenter, vmax=vmax)
        cbarTicks = np.linspace(0, vmax, 10)
    else:
        norm = TwoSlopeNorm(vmin=vcenter - maxDist, vcenter=vcenter,
                            vmax=vcenter + maxDist)
        cbarTicks = np.linspace(vcenter - maxDist, vcenter + maxDist, 10)
    roundTo = 2
    while len(np.unique(cbarTicks)) != len(np.unique(np.round(cbarTicks,
                                                              roundTo))):
        roundTo += 1
    cbarTicks = np.round(cbarTicks, roundTo)
    
    # 2. Create a row of subplots, one per category
    fig, axes = plt.subplots(1, 1, figsize=figsize, squeeze=False, 
                             constrained_layout=True)
    axes = axes.flatten()
    for ax, (_, dataList) in zip(axes, dataDict.items()):
        heatmapDat = pd.DataFrame(dataList["means"])
        heatmapDat.columns = np.round(dataList["testpoints"], 2)
        heatmapDat.index = np.round(covChange, 2)
        if inCI:
            heatmapDat = heatmapDat.mul(1 - dataList["nullInCI"])
            heatmapDat = heatmapDat.replace(0, 1)
        # X-axis ticks
        xPos = np.linspace(0, heatmapDat.shape[1]-1, numXticks, dtype=int)
        xLabels = [heatmapDat.columns[i] for i in xPos]

        # Y-axis ticks
        yPos = np.linspace(0, heatmapDat.shape[0]-1, numYticks, dtype=int)
        yLabels = [heatmapDat.index[i] for i in yPos]
        makeHeatmapPlot(ax=ax, heatmapDat=heatmapDat, norm=norm, cmap=cmap, 
                        cat="", fontsize=fontsize, xlabel=xlabel, ylabel=ylabel,
                        cbar=True, cbarLabel=cbarLabel, cbarTicks=cbarTicks, 
                        xPos=xPos, yPos=yPos, xLabels=xLabels, yLabels=yLabels)
    
    # 4. Adjust layout, add suptitle, save
    fig.suptitle(suptitle, fontsize=titleSize, y=suptitleLoc[1], 
                 x=suptitleLoc[0])
    plt.close(fig)

    return fig