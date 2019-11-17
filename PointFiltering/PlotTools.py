#!/usr/bin/env python
"""
Set of tools to optimize plotting
"""
import numpy
import matplotlib.pyplot as plt
from . import hillshade

def stretch(values):
    """ Calculate stretched values based on mean and 2 standard deviations """
    values = values[numpy.isnan(values) == 0].flatten()
    mean = values.mean()
    stdev = values.std()
    if mean-2*stdev < values.min():
        low = values.min()
    else:
        low = mean-2*stdev
    if mean+2*stdev > values.max():
        high = values.max()
    else:
        high = mean+2*stdev
    return low, high


def scatter(x, y, c, clim, fn, title=None):
    fig, ax = plt.subplots()
    ax.scatter(x, y, c=c, s=3, vmin=clim[0], vmax=clim[1])
    ax.axis("equal")
    ax.set_title(title)
    fig.savefig(fn, dpi=150)


def image(grid, fn, do_shade=False, res=None, title=None):
    fig, ax = plt.subplots()
    ax.imshow(grid, clim=stretch(grid))
    if do_shade:
        shade = hillshade.hillshade(grid, res, res)
        ax.imshow(shade, clim=stretch(shade), cmap='gray', alpha=0.5)
    ax.axis("equal")
    ax.set_title(title)
    fig.savefig(fn, dpi=150)

