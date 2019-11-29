#!/usr/bin/env python
"""
This file includes several methods for general use.
"""
import numpy as np
import ogr
from .Reference import ReferenceData
from . import Grid
import time
from scipy import interpolate
from osgeo import gdal


def openshp(fname, c=False):
    """ Read shapefile and extract x,y points from the geometry and z from the attribute table """
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.Open(fname, 0)
    lyr = ds.GetLayer(0)
    proj = lyr.GetSpatialRef().ExportToProj4()
    n_shapes = lyr.GetFeatureCount()

    reference = ReferenceData(n_shapes, proj)
    if c:
        reference.classification = []

    i = 0
    for shape in lyr:
        reference.point_x[i] = shape.geometry().GetX()
        reference.point_y[i] = shape.geometry().GetY()
        reference.point_z[i] = shape.GetField("Z")
        if c:
            reference.classification.append(shape.GetField("Class"))
        i += 1
    return reference


def make_grid(points, res):
    """ align grid """
    range_x = (points.point_x.min(), points.point_x.max())
    range_y = (points.point_y.min(), points.point_y.max())

    left = np.round(range_x[0] - np.mod(range_x[0], res), 2)
    right = np.round(range_x[1] + res - np.mod(range_x[1] + res, res), 2)
    bottom = np.round(range_y[0] - np.mod(range_y[0], res), 2)
    top = np.round(range_y[1] + res - np.mod(range_y[1] + res, res), 2)

    x_centers = np.arange(left + res / 2, right, res)
    y_centers = np.arange(bottom + res / 2, top, res)

    grid_x, grid_y = np.meshgrid(x_centers, y_centers, indexing="xy")
    return [grid_x, grid_y]


def correction(grid, gps_data, vi_threshold, s_err = None):
    in_grid = grid.in_grid(gps_data.point_x, gps_data.point_y)[0]
    values = np.isnan(grid.point_dtm) == False
    zi = Grid.idw(grid.point_x[values], grid.point_y[values], grid.point_dtm[values], gps_data.point_x, gps_data.point_y, 3)

    dtm_diff_before = (gps_data.point_z - zi)[in_grid]

    # GPS ground points
    gps_data.point_vi = Grid.idw(grid.point_x, grid.point_y, grid.point_vi, gps_data.point_x, gps_data.point_y, 5)
    gps_data.ground = gps_data.point_vi > vi_threshold

    t0 = time.time()
    if len(gps_data.ground) > 0:
        # correction model
        if s_err is None:
            spline = interpolate.SmoothBivariateSpline(gps_data.point_x[gps_data.ground], gps_data.point_y[gps_data.ground], dtm_diff_before[gps_data.ground], kx=1, ky=1)
            corr_plane = np.flipud(spline.ev(grid.grid_x, grid.grid_y))
            # correction
            corr = Grid.idw(grid.grid_x, grid.grid_y, corr_plane, grid.point_x, grid.point_y, 3)
        else:
            corr = s_err
        grid.point_dtm = grid.point_dtm + corr
        grid.point_dsm = grid.point_dsm + corr
    else:
        print("No ground-based GPS data found within grid boundaries.")
    if len(gps_data.point_x) > 0:
        zi = Grid.idw(grid.point_x[values], grid.point_y[values], grid.point_dtm[values], gps_data.point_x, gps_data.point_y, 3)
        dtm_diff_after = gps_data.point_z - zi
        print("avg error before: %.2f, after: %.2f" % (np.mean(np.abs(dtm_diff_before)),
                                                        np.mean(np.abs(dtm_diff_after))))
    else:
        print("No GPS data found within grid boundaries.")
    print("Correcting with with gps data took %d seconds" % (time.time() - t0))
    return grid, dtm_diff_before, dtm_diff_after


def savergb(grid, raster, filename, wkt=''):
    geotransform = (grid.grid_x.min() - grid.res / 2, grid.res, 0.0, grid.grid_y.max() + grid.res / 2, 0.0, -grid.res)
    rows, cols = grid.grid_x.shape

    # change NaN to -999 no data value
    raster[np.isnan(raster)] = -999

    # get driver
    driver = gdal.GetDriverByName('HFA')
    # create file
    out_data = driver.Create(filename, cols, rows, 3, gdal.GDT_Float32)

    # get band 1 to write to
    out_band = out_data.GetRasterBand(1)
    out_band.WriteArray(raster[:, :, 0], 0, 0)

    # get band 2 to write to
    out_band = out_data.GetRasterBand(2)
    out_band.WriteArray(raster[:, :, 1], 0, 0)

    # get band 3 to write to
    out_band = out_data.GetRasterBand(3)
    out_band.WriteArray(raster[:, :, 2], 0, 0)

    out_band.FlushCache()
    # set no_data
    out_band.SetNoDataValue(-999)
    # calculate statistics
    out_band.GetStatistics(0, 1)
    # georeference
    out_data.SetGeoTransform(geotransform)
    out_data.SetProjection(wkt)
    out_band.FlushCache()
    out_data = None


def saveimg(grid, raster, filename, wkt=''):
    geotransform = (grid.grid_x.min() - grid.res / 2, grid.res, 0.0, grid.grid_y.max() + grid.res / 2, 0.0, -grid.res)
    rows, cols = grid.grid_x.shape

    # change NaN to -999 nodata value
    raster[np.isnan(raster)] = -999
    # get driver
    driver = gdal.GetDriverByName('GTiff')
    # create file
    out_data = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)
    # get band to write to
    out_band = out_data.GetRasterBand(1)
    # write array
    out_band.WriteArray(raster, 0, 0)
    out_band.FlushCache()
    # set nodata
    out_band.SetNoDataValue(-999)
    # calculate statistics
    out_band.GetStatistics(0,1)
    # georeference

    out_data.SetGeoTransform(geotransform)
    out_data.SetProjection(wkt)
    out_data = None

