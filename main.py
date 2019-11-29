#!/usr/bin/env python3
"""
SfM point filtering
"""
import sys
import time
import numpy
import glob
from termcolor import colored
import pylab as plt

from PointFiltering import LazData, VegetationIndices, Filter, Grid, PlotTools
from PointFiltering.Helpers import saveimg, savergb

def run(filter, **kwargs):
    """ Run point filter
    Args:
        filter (string): options are 'isl', 'vi', 'isl_vi'
        **kwargs (dict): optional keyword arguments
    Raises:
        ValueError: when invalid filter or VI is provided
    """
    print('Run {} filter'.format(filter))

    FOLDER = kwargs.get('folder', 'example_data/')
    RES = kwargs.get('res', 0.25)
    VI = kwargs.get('vi', 'exg')
    PLOT = kwargs.get('plot', False)
    wkt = kwargs.get('wkt', '')

    if VI == 'exg':
        THRESHOLD = kwargs.get('threshold', 0.01)
    elif VI == 'vvi':
        THRESHOLD = kwargs.get('threshold', 5e-5)
    else:
        raise ValueError('invalid VI')

    OK = "green"

    if not filter or filter not in Filter.FILTERS:
        raise ValueError('invalid filter')

    nr_points = 0
    nr_ground_points = 0

    i = 0
    LISTDIR = glob.glob(FOLDER + "*.las")
    LISTDIR.sort()
    if len(LISTDIR) == 0:
        raise RuntimeError('No las files found in folder')
    print("Processing tiles")

    # iterate over tiles
    for LASFILE in LISTDIR:
        print(LASFILE+"..."),
        points = LazData.createlaz(LASFILE)
        nr_points += len(points.point_x)
        if len(points.ground) == 0:
            print ("no points in tile, skipping"),
        else:
            try:
                # In the first iteration all points are considered ground points.
                # The ground point labels are stored as an array in points.ground as bool
                # each filter sets the label to False when points are considered not-ground

                print ("total: %6d" % len(points.ground)),
                # first filter based on vegetation index
                if VI == "vvi":
                    points.vegindex = VegetationIndices.vvi(points)
                    if filter in ['vi', 'isl_vi']:
                        Filter.filter_by_vvi(points, THRESHOLD)
                if VI == "exg":
                    points.vegindex = VegetationIndices.exg(points)
                    if filter in ['vi', 'isl_vi']:
                        Filter.filter_by_exg(points, THRESHOLD)
                # then filter by ISL
                if filter in ['isl', 'isl_vi']:
                    Filter.filter_by_dtm(points)
                print ("ground: %d" % sum(points.ground)),
                nr_ground_points += sum(points.ground)

            except Exception as err:
                print (err)
            finally:
                if "GRID" not in locals():
                    GRID = Grid.Grid(points, res=RES)
                    GRID.filter_method = points.filter_method
                else:
                    GRID_PART = Grid.Grid(points, res=RES)
                    GRID.merge_grids(GRID_PART)
        i += 1
        print(colored("done", OK)),

    print('total nr points: {}; nr ground points: {}'.format(
        nr_points, nr_ground_points))
    # create raster files
    t0 = time.time()
    print("rasterizing %d points into a (%d, %d) grid..." %
           (len(GRID.point_x), GRID.grid_x.shape[0], GRID.grid_x.shape[1])),
    GRID.create_dtm()
    GRID.create_dsm()
    GRID.create_rgb()
    GRID.create_vi()
    print(colored("done", OK)),

    chm = GRID.dsm - GRID.dtm
    chm[chm < 0] = 0

    print(colored("done", OK))

    if PLOT:
        print ("plotting..."),
        # points vi
        values = numpy.isnan(GRID.point_vi) == 0
        data = GRID.point_vi[values]
        clim = PlotTools.stretch(data)
        filename = "out/plots/points_vi.png"
        PlotTools.scatter(
            GRID.point_x[values],
            GRID.point_y[values],
            data,
            clim,
            filename,
            title='VI points'
        )
        # points gvi
        values = numpy.isnan(GRID.point_dtm) == 0
        data = GRID.point_gvi[values]
        # use clim from vi
        filename = "out/plots/points_vi-ground" + GRID.filter_method + ".png"
        PlotTools.scatter(
            GRID.point_x[values],
            GRID.point_y[values],
            data,
            clim,
            filename,
            title='VI (ground only)'
        )
        # points dsm -----------------------------
        values = numpy.isnan(GRID.point_dsm) == 0
        data = GRID.point_dsm[values]
        clim = PlotTools.stretch(data)
        fn = "out/plots/points_dsm" + GRID.filter_method + ".png"
        PlotTools.scatter(
            GRID.point_x[values],
            GRID.point_y[values],
            data,
            clim,
            fn,
            title='DSM points'
        )

        # ground points z ------------------------
        values = numpy.isnan(GRID.point_dtm) == 0
        PlotTools.scatter(
            GRID.point_x[values],
            GRID.point_y[values],
            GRID.point_dtm[values],
            clim,
            "out/plots/points_dtm" + GRID.filter_method + ".png",
            title='DTM points'
        )
        # dsm
        PlotTools.image(
            GRID.dsm,
            "out/plots/grid_dsm.png",
            do_shade=True,
            res=RES,
            title='DSM grid'
        )
        # dtm
        PlotTools.image(
            GRID.dtm,
            "out/plots/grid_dtm" + GRID.filter_method + ".png",
            do_shade=True,
            res=RES,
            title='DTM grid'
        )
        # chm
        PlotTools.image(
            chm,
            "out/plots/grid_chm" + GRID.filter_method + ".png",
            title='CHM grid'
        )
        print(colored("done", OK))

    # save images
    saveimg(GRID, GRID.dtm, 'out/dtm.tif', wkt=wkt)
    saveimg(GRID, GRID.dsm, 'out/dsm.tif', wkt=wkt)
    saveimg(GRID, GRID.dsm-GRID.dtm, 'out/chm.tif', wkt=wkt)
    saveimg(GRID, GRID.raster_vi, 'out/vi.tif', wkt=wkt)
    savergb(GRID, GRID.rgb, 'out/rgb.tif', wkt=wkt)


if __name__ == "__main__":
    args = sys.argv

    kwargs = {}
    # get filter from command line argument
    if len(args) > 1:
        for arg in args:
            if '=' not in arg:
                continue
            key = arg.split('=')[0]
            value = arg.split('=')[1]
            if 'filter' in arg:
                filter = value.lower()
            else:
                kwargs[key] = value
    # use default if not provided
    if 'filter' not in locals():
        filter = 'isl_vi'
    wkt = 'PROJCS["WGS 84 / UTM zone 32N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",9],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32632"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'
    run(filter=filter, plot=True, wkt=wkt, **kwargs)

