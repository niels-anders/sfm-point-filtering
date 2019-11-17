#!/usr/bin/env python
"""
Documentation
"""

import numpy


class ReferenceData:
    """ Reference (x, y, z) point data and proj4 info """
    def __init__(self, nrshapes, proj):
        self.point_x = numpy.zeros(nrshapes)
        self.point_y = numpy.zeros(nrshapes)
        self.point_z = numpy.zeros(nrshapes)
        self.point_vi = None
        self.ground = None
        self.classification = None
        self.proj = proj

    def remove_outside(self, grid):
        """ Remove points outside grid """
        inside = numpy.where((self.point_x > grid.extent.xmin) & (self.point_x < grid.extent.xmax) & (self.point_y > grid.extent.ymin) & (self.point_y < grid.extent.ymax))
        self.point_x = self.point_x[inside]
        self.point_y = self.point_y[inside]
        self.point_z = self.point_z[inside]
        self.classification = numpy.array(self.classification)[inside]



