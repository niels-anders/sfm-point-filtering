#!/usr/bin/env python
"""
This file includes methods for accuracy assessment based on ground truth data
"""

from scipy.interpolate import interp2d
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from termcolor import colored


def acc(grid, ref_data):
    """ check grid with reference data"""
    xi = grid.grid_x[0, :]
    yi = grid.grid_y[:, 0]
    dtm = interp2d(xi, yi, np.flipud(grid.dtm), kind='linear')
    dsm = interp2d(xi, yi, np.flipud(grid.dsm), kind='linear')
    compare = np.zeros((len(ref_data.point_z), 5))
    for i in range(0, len(ref_data.point_z)):
        compare[i, 0] = ref_data.point_z[i]  # GPS
        compare[i, 1] = dtm(ref_data.point_x[i], ref_data.point_y[i])[0]  # DTM
        compare[i, 3] = dsm(ref_data.point_x[i], ref_data.point_y[i])[0]  # DSM
    compare[:, 2] = compare[:, 1] - compare[:, 0]  # DTM_DIFF
    compare[:, 4] = compare[:, 3] - compare[:, 0]  # DSM_DIFF
    return compare


def boxplot(grid, gps_data):
    if len(gps_data.classification) == 0:
        colored("WARN: not able to assess accuracy: no reference points inside bounding box", "yellow")
        return

    compare = acc(grid, gps_data)

    df = pd.DataFrame(data=compare, index=gps_data.classification, columns=list(["GPS", "DTM", "DTM_DIFF", "DSM", "DSM_DIFF"]))
    groups = df.groupby(df.index)
    fig, ax = plt.subplots()

    data = []
    labels = []
    for label, group in groups:
        data.append(group.DTM_DIFF)
        labels.append(label)
    ax.boxplot(data, labels=labels)
    ax.set_title("filtered: " + grid.filter_method)
    ax.set_ylim([-1, 10])
    fig.savefig("out/plots/boxplot" + grid.filter_method + ".png", dpi=150)


def compare_as_df(grid, gps_data):
    compare = acc(grid, gps_data)
    df = pd.DataFrame(data=compare, index=gps_data.classification,
                      columns=list(["GPS", "DTM", "DTM_DIFF", "DSM", "DSM_DIFF"]))
    return df


def categorized_boxplot(df):
    titles = ['Bare', 'Shrubs', 'Trees']
    names = df.columns
    fig, axes = plt.subplots(1, 3)
    i = 0
    for title in titles:
        dfs = df.loc[title]
        dfs.plot(kind='box', ax=axes[i])
        axes[i].set_xticklabels(names, rotation=90)
        axes[i].spines['right'].set_visible(False)  # Zet deze aan voor hele box
        axes[i].spines['top'].set_visible(False)  # Zet deze aan voor hele box
        axes[i].set_title(title)
        i += 1
    axes[0].set_ylabel('Elevation Difference (m)')
    fig.set_size_inches(8.27, 3.5)  # hele breedte a4, ~1/3 hoogte
    fig.tight_layout()
    #plt.show()
    plt.savefig('out/plots/boxplots_comparison.png', dpi=300)

    # Maak een DF met alle MSE waardes
    dfnames = [x for x in df.columns if not x == 'category']
    namedict = dict(zip(dfnames, names))
    errordf = pd.DataFrame(index=names, columns=titles)
    calcrmse = lambda x: np.sqrt(np.sum(np.power(x, 2)) / len(x))
    for filters in dfnames:
        for title in titles:
            mse = calcrmse(df.loc[title][filters])
            errordf.loc[namedict[filters], title] = mse
    errordf.to_csv('RMSE.csv')


