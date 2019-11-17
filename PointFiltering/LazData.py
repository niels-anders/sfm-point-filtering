#!/usr/bin/env python
"""
This script reads laz/las files, creates a dataset from the points
and reads/writes to and from pickle files.
"""
try:
    import cPickle
except ImportError:
    import pickle
import laspy
import numpy


class Lazdata:
    # pylint: disable=C0103, C1001, R0902
    """ Dataset object with x, y, z coordinates, color values """

    def __init__(self, las):
        self.point_x = las.x
        self.point_y = las.y
        self.point_z = las.z
        self.red = las.red
        self.green = las.green
        self.blue = las.blue
        self.classification = las.Classification
        self.vegindex = numpy.ones(len(las.x))
        self.ground = numpy.ones(len(las.x)).astype("bool")
        self.filter_method = ""

    def serialize(self, pickle):
        """ Write to .pickle file"""
        with open(pickle, 'wb') as output:
            cPickle.dump(self, output, cPickle.HIGHEST_PROTOCOL)

    def merge_points(self, other):
        """ Merge two Datasets """
        self.point_x = numpy.append(self.point_x, other.point_x)
        self.point_y = numpy.append(self.point_y, other.point_y)
        self.point_z = numpy.append(self.point_z, other.point_z)
        self.red = numpy.append(self.red, other.red)
        self.green = numpy.append(self.green, other.green)
        self.blue = numpy.append(self.blue, other.blue)
        self.classification = numpy.append(
            self.classification, other.classification)


def deserialize(pickled):
    """ Read .pickle file and return Dataset object"""
    with open(pickled, "rb") as pickle_object:
        lazdata = cPickle.load(pickle_object)
    return lazdata


def createlaz(lazfile):
    """ Import LAZ/LAS files using laspy library"""
    laz = laspy.file.File(lazfile, mode='r')
    return Lazdata(laz)


def get_boundingbox(lazfile):
    """ Print bounding box """
    laz = createlaz(lazfile)
    print ("%s: [%.2f, %.2f, %.2f, %.2f]" % (lazfile, laz.point_x.min(), laz.point_x.max(), laz.point_y.min(), laz.point_y.max()))

"""
if __name__ == "__main__":
    IN = sys.argv[1]
    BASENAME = IN.split(".")[0]
    OUT = BASENAME + ".pickle"

    POINTS = createlaz(IN)
    POINTS.serialize(OUT)
"""
