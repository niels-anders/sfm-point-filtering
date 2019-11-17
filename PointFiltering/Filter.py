#!/usr/bin/env python
"""
This file includes several filter methods.
"""
import numpy as np
from scipy.interpolate import interp2d

from . import VegetationIndices
from . import Grid
from .Helpers import make_grid


FILTERS = ('isl', 'isl_vi', 'vi', 'tin', 'ps', 'low')


def filter_by_exg(points, threshold):
    """ Use ExG to classify ground points """
    vegindex = VegetationIndices.exg(points)
    points.ground = vegindex < threshold
    points.filter_method += "_exg"


def filter_by_vvi(points, threshold):
    """ Use VVI to classify ground points """
    vegindex = VegetationIndices.vvi(points)
    points.ground = vegindex > threshold
    points.filter_method += "_vvi"


def filter_by_tin(points):
    """ Copy LAStools classification """
    points.ground = points.classification == 2
    points.filter_method += "_tin"


def filter_by_ps(points):
    """ Copy PhotoScan classification """
    points.ground = points.classification == 2
    points.filter_method += "_ps"


def filter_by_lowest(points):
    points.filter_method += "_low"


def filter_by_dtm(points):
    """
    checks whether ground points are on top of bumps
    """
    # get terrain height
    res = 2
    [x_grid, y_grid] = make_grid(points, res)
    xi = x_grid[0, :]
    yi = y_grid[:, 0]
    neighbors = 5

    max_iterations = 3
    min_ground_points = 10
    threshold_above_surface = 0.1  # m
    it = 0
    prev_nr_ground_points = -1
    while it < max_iterations:
        if (sum(points.ground) < min_ground_points) & \
            (sum(points.ground) != prev_nr_ground_points):
            break
        try:
            prev_nr_ground_points = sum(points.ground)

            # make grid from all ground points
            dtm = Grid.idw(
                points.point_x[points.ground],
                points.point_y[points.ground],
                points.point_z[points.ground],
                x_grid,
                y_grid,
                neighbors
            ).reshape(x_grid.shape)

            # get grid elevation at point location
            dtm[np.isnan(dtm) == 1] = 999
            f = interp2d(xi, yi, dtm, kind='linear')
            zi = np.zeros(points.point_x.shape)
            for i in np.arange(len(zi)):
                zi[i] = f(points.point_x[i], points.point_y[i])

            # determine vegetation points point elevation is higher than grid + threshold
            veg = (points.point_z > zi + threshold_above_surface)

            # remove ground label from vegetation points
            points.ground[points.ground & veg] = False

            # prepare for next iteration
            print('iteration {}: nr of ground points from {} to {}'
                  .format(it, prev_nr_ground_points, sum(points.ground)))
            it += 1

        except:
            raise RuntimeError("filter by dtm failed")
    points.filter_method += "_dtm"




