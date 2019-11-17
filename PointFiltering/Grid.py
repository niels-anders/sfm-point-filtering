#!/usr/bin/env python
"""
This script creates a grid out of points.
"""
try:
    import cPickle as pickle
except ImportError:
    import pickle
import numpy as np
import gdal

from . import invdisttree
from .import VegetationIndices


class Grid:
    """ A regularly spaced grid, structured as flattened arrays """

    def __init__(self, las, res=0.1, low=False):
        """ Initialization """
        self.res = res
        self.point_x = np.array([])
        self.point_y = np.array([])
        self.point_r = np.array([])
        self.point_g = np.array([])
        self.point_b = np.array([])
        self.point_dsm = np.array([])
        self.point_dtm = np.array([])
        self.point_vi = np.array([])
        self.point_gvi = np.array([])
        self.grid_x = None
        self.grid_y = None
        self.dsm = None
        self.dtm = None
        self.rgb = None
        self.raster_vi = None
        self.pointdensity = np.array([])
        self.ground_pointdensity = np.array([])
        self.low = low
        make_grid(self, las=las)
        average(self, res, las)
        self.extent = Extent(las)
        self.mask = None
        self.filter_method = None


    def create_dtm(self):
        """
        Rasterize the grid and fill gaps with idw interpolation.
        Returns a 2D matrix with z values
        """
        nnear = 5
        values = np.isnan(self.point_dtm) == 0
        self.dtm = np.flipud(idw(self.point_x[values], self.point_y[values], self.point_dtm[values],
                                 self.grid_x, self.grid_y, nnear).reshape(self.grid_x.shape))

    def create_dsm(self):
        """
        Rasterize the grid and fill gaps with idw interpolation.
        Returns a 2D matrix with z values
        """
        nnear = 5
        values = np.isnan(self.point_dsm) == 0
        self.dsm = np.flipud(idw(self.point_x[values], self.point_y[values], self.point_dsm[values],
                                 self.grid_x, self.grid_y, nnear).reshape(self.grid_x.shape))

    def create_rgb(self):
        """
        Rasterize the grid and fill gaps with idw interpolation.
        Returns a 2D matrix with z values
        """
        nnear = 5
        values = np.isnan(self.point_r) == 0
        red = np.flipud(idw(self.point_x[values], self.point_y[values], self.point_r[values], self.grid_x, self.grid_y,
                            nnear).reshape(self.grid_x.shape))
        green = np.flipud(idw(self.point_x[values], self.point_y[values], self.point_g[values], self.grid_x,
                              self.grid_y, nnear).reshape(self.grid_x.shape))

        blue = np.flipud(idw(self.point_x[values], self.point_y[values], self.point_b[values], self.grid_x,
                             self.grid_y, nnear).reshape(self.grid_x.shape))

        [red, green, blue] = VegetationIndices.convert_colors(red, green, blue)

        rgb = np.array([red, green, blue])
        rgb = np.rollaxis(rgb.T, 2).T
        self.rgb = rgb

    def create_vi(self):
        """
        Rasterize the grid and fill gaps with idw interpolation.
        Returns a 2D matrix with z values
        """
        nnear = 5
        self.raster_vi = np.flipud(idw(
            self.point_x, self.point_y, self.point_vi, self.grid_x, self.grid_y, nnear).reshape(self.grid_x.shape))

    def merge_grids(self, subgrid):
        """ Merge two grids """
        # append points
        self.point_x = np.append(self.point_x, subgrid.point_x)
        self.point_y = np.append(self.point_y, subgrid.point_y)
        self.point_r = np.append(self.point_r, subgrid.point_r)
        self.point_g = np.append(self.point_g, subgrid.point_g)
        self.point_b = np.append(self.point_b, subgrid.point_b)
        self.point_dsm = np.append(self.point_dsm, subgrid.point_dsm)
        self.point_dtm = np.append(self.point_dtm, subgrid.point_dtm)
        self.point_vi = np.append(self.point_vi, subgrid.point_vi)
        self.point_gvi = np.append(self.point_gvi, subgrid.point_gvi)
        self.pointdensity = np.append(self.pointdensity, subgrid.pointdensity)
        self.ground_pointdensity = np.append(self.ground_pointdensity, subgrid.ground_pointdensity)
        # update grid & extent
        make_grid(self, las=None)
        self.extent = Extent(self)

    def export_tiff(self, raster, tiff, mask=None):
        """
        Export to geotiff
        """
        raster[np.isnan(raster)] = -999
        driver = gdal.GetDriverByName('GTiff')
        rows, cols = self.grid_x.shape
        if len(raster.shape) == 2:
            nr_of_bands = 1
        elif len(raster.shape) == 3:
            nr_of_bands = raster.shape[-1]
        out = driver.Create(tiff, cols, rows, nr_of_bands, gdal.GDT_Float32)
        # get band to write to
        try:
            for band_nr in np.arange(nr_of_bands):
                band = out.GetRasterBand(band_nr+1)
                # write array
                if nr_of_bands > 1:
                    band.WriteArray(raster[:, :, band_nr], 0, 0)
                else:
                    band.WriteArray(raster, 0, 0)
                band.FlushCache()
                # set no data
                band.SetNoDataValue(-999)
                # calculate statistics
                band.GetStatistics(0, 1)
            # geo_transform = (left, cell_size, 0.0, top, 0.0, -cell_size)
            geo_transform = (self.grid_x.min() - self.res / 2, self.res, 0.0,
                             self.grid_y.max() + self.res / 2, 0.0, -self.res)
            # georeference
            out.SetGeoTransform(geo_transform)
            out.SetProjection("")

            if mask is not None:
                alpha = out.GetRasterBand(nr_of_bands+1)
                alpha.WriteArray(mask, 0, 0)

            # empty memory
            out = None
        except:
            raise ("can't write to " + tiff)

    def in_grid(self, points_x, points_y):
        """ return true when point falls in bounding box of the grid """
        return np.where(
            (points_x > self.extent.xmin) & (points_x < self.extent.xmax) & (points_y > self.extent.ymin) & (
                points_y < self.extent.ymax))

    def serialize(self, pickle):
        """ Write to .pickle file"""

        with open(pickle, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)


class Extent:
    """ Min/max values of point coordinates """

    def __init__(self, data):
        """ initialization """
        self.xmin = np.nanmin(data.point_x)
        self.xmax = np.nanmax(data.point_x)
        self.ymin = np.nanmin(data.point_y)
        self.ymax = np.nanmax(data.point_y)


def make_grid(grid, las):
    """ align grid """
    if las is not None:
        range_x = (las.point_x.min(), las.point_x.max())
        range_y = (las.point_y.min(), las.point_y.max())
    else:
        range_x = (np.nanmin(grid.point_x), np.nanmax(grid.point_x))
        range_y = (np.nanmin(grid.point_y), np.nanmax(grid.point_y))

    left = np.round(range_x[0] - np.mod(range_x[0], grid.res), 2)
    right = np.round(range_x[1] + grid.res - np.mod(range_x[1] + grid.res, grid.res), 2)
    bottom = np.round(range_y[0] - np.mod(range_y[0], grid.res), 2)
    top = np.round(range_y[1] + grid.res - np.mod(range_y[1] + grid.res, grid.res), 2)

    x_centers = np.arange(left + grid.res / 2, right, grid.res)
    y_centers = np.arange(bottom + grid.res / 2, top, grid.res)

    grid_x, grid_y = np.meshgrid(x_centers, y_centers, indexing="xy")
    grid.grid_x = grid_x
    grid.grid_y = grid_y


def average(grid, res, las):
    """ Calculate average values on (x,y) centroid locations """
    # TODO: unit test
    # determine cell borders
    x_borders = np.arange(
        grid.grid_x.min() - res / 2, grid.grid_x.max() + res, res)
    y_borders = np.arange(
        grid.grid_y.min() - res / 2, grid.grid_y.max() + res, res)
    nr_cells = len(grid.grid_x.flatten())

    # initialize arrays
    grid.point_x = np.zeros(nr_cells) * np.NaN
    grid.point_y = np.zeros(nr_cells) * np.NaN
    grid.point_r = np.zeros(nr_cells) * np.NaN
    grid.point_g = np.zeros(nr_cells) * np.NaN
    grid.point_b = np.zeros(nr_cells) * np.NaN
    grid.point_dsm = np.zeros(nr_cells) * np.NaN
    grid.point_dtm = np.zeros(nr_cells) * np.NaN
    grid.point_vi = np.zeros(nr_cells) * np.NaN
    grid.point_gvi = np.zeros(nr_cells) * np.NaN
    grid.pointdensity = np.zeros(nr_cells)
    grid.ground_pointdensity = np.zeros(nr_cells)
    k = 0
    # fill arrays
    for i in np.arange(len(x_borders) - 1):
        query = (las.point_x > x_borders[i]) & (las.point_x < x_borders[i + 1])
        col_x = np.compress(query, las.point_x)
        col_y = np.compress(query, las.point_y)
        col_z = np.compress(query, las.point_z)
        col_v = np.compress(query, las.vegindex)
        col_r = np.compress(query, las.red)
        col_g = np.compress(query, las.green)
        col_b = np.compress(query, las.blue)
        col_ground = np.compress(query, las.ground)
        for j in np.arange(len(y_borders) - 1):
            query = (col_y > y_borders[j]) & (col_y < y_borders[j + 1])
            cell_x = np.compress(query, col_x)
            cell_y = np.compress(query, col_y)
            cell_z = np.compress(query, col_z)
            cell_v = np.compress(query, col_v)
            cell_r = np.compress(query, col_r)
            cell_g = np.compress(query, col_g)
            cell_b = np.compress(query, col_b)
            ground = np.compress(query, col_ground)
            if len(cell_v) > 0:
                grid.point_x[k] = cell_x.mean()
                grid.point_y[k] = cell_y.mean()
                grid.point_dsm[k] = cell_z.mean()
                grid.point_vi[k] = cell_v.mean()
                grid.point_r[k] = cell_r.mean()
                grid.point_g[k] = cell_g.mean()
                grid.point_b[k] = cell_b.mean()
                grid.pointdensity[k] = len(cell_z)
                if len(cell_z[np.argwhere(ground == 1)]) > 0:
                    # grid.point_x[k] = cell_x[np.argwhere(ground == 1)].mean()
                    # grid.point_y[k] = cell_y[np.argwhere(ground == 1)].mean()
                    # grid.point_dtm[k] = np.percentile(cell_z[np.argwhere(ground == 1)], 25)
                    if grid.low:
                        grid.point_dtm[k] = cell_z[np.argwhere(ground == 1)].min()
                    else:
                        grid.point_dtm[k] = cell_z[np.argwhere(ground == 1)].mean()
                    grid.point_gvi[k] = cell_v[np.argwhere(ground == 1)].mean()
                    grid.ground_pointdensity[k] = len(cell_z[ground == 1])
            k += 1
    values = np.argwhere(grid.pointdensity > 0)
    # clip out nans
    # NB: point_dtm is still nan at zero ground points
    grid.point_x = grid.point_x[values]
    grid.point_y = grid.point_y[values]
    grid.point_r = grid.point_r[values]
    grid.point_g = grid.point_g[values]
    grid.point_b = grid.point_b[values]
    grid.point_dsm = grid.point_dsm[values]
    grid.point_dtm = grid.point_dtm[values]
    grid.point_vi = grid.point_vi[values]
    grid.point_gvi = grid.point_gvi[values]
    grid.extent = Extent(las)


def alt_idw(point_x, point_y, point_z, grid_x, grid_y, nnear, power):
    """ alternative idw """
    values = np.zeros(grid_x.shape)
    for i in range(0, len(grid_x)):
        dist_x = point_x - grid_x[i]
        dist_y = point_y - grid_y[i]
        dist = np.sqrt(dist_x ** 2 + dist_y ** 2)
        k = dist.argsort()[:nnear]
        dist_i = dist[k]
        weight = 1 / (dist_i ** power)
        values[i] = sum(point_z[k] * weight) / sum(weight)
    return values


def idw(point_x, point_y, point_z, xi, yi, nnear):
    """ Inverse distance weighting """
    # TODO: unit test
    # flatten variables
    point_x = point_x.flatten()
    point_y = point_y.flatten()
    point_z = point_z.flatten()
    xi = xi.flatten()
    yi = yi.flatten()
    # perform idw
    xxyy = np.append([point_x], [point_y], axis=0).transpose()
    xiyi = np.append([xi], [yi], axis=0).transpose()
    tree = invdisttree.Invdisttree(xxyy, point_z, leafsize=10, stat=1)
    weights = None  # np.ones(point_z.shape)
    return tree(xiyi, nnear=nnear, eps=0.5, p=2)


def nearest(point_x, point_y, grid_x, grid_y, grid_z):
    """ nearest neighbor """
    point_z = np.zeros(point_x.shape)
    for pnt in np.arange(len(point_x)):
        left = np.argwhere(point_x[pnt] < grid_x[0])
        print (pnt, point_x[pnt], point_y[pnt], left, grid_x)
        if len(left) > 1:
            if (point_x[pnt] - grid_x[0][left[0]][0]) < (grid_x[0][left[1]][0]) - point_x[pnt]:
                xii = left[0][0]
            else:
                xii = left[1][0]
        if len(left) == 1:
            xii = left[0][0]

        bottom = np.argwhere(point_y[pnt] < grid_y[:, 0])
        if len(bottom) > 1:
            if (point_y[pnt] - grid_y[:, 0][bottom[0]][0]) < (grid_y[:, 0][bottom[1]][0] - point_y[pnt]):
                yii = bottom[0][0]
            else:
                yii = bottom[1][0]
        if len(bottom) == 1:
            yii = bottom[0][0]
        if (len(bottom) > 0) & (len(left) > 0):
            point_z[pnt] = grid_z[xii][yii]
    return point_z


def create_mask(grid):
    """ Create mask of data points in grid """
    # TODO: create unit test
    return (grid.pointdensity > 0).reshape(grid.grid_x.shape);
