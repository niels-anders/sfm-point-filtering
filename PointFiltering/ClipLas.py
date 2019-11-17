import laspy
from osgeo import gdal


class BoundingBox:
    def __init__(self, x_min, y_min, x_max, y_max):
        self.x_min = x_min
        self.y_min = y_min
        self.x_max = x_max
        self.y_max = y_max


def clip(in_file, out_file, area_file):
    """
    Clip las points inside bounding box and write to new las file
    :param in_file: input las filename
    :param out_file: output las filename
    :param area_file: input tiff file to get bounding box
    :return: None
    """
    bb = load_bb(area_file)
    in_las = laspy.file.File(in_file, mode='r')
    out_las = laspy.file.File(out_file, mode='w', header=in_las.header)

    inside = (in_las.x > bb.x_min) & (in_las.x < bb.x_max) & (in_las.y > bb.y_min) & (in_las.y < bb.y_max)
    out_las.points = in_las.points[inside]
    out_las.close()


def load_bb(filename):
    """
    Loads Get the bounding box from a geotiff
    :param filename: input tiff filename
    :return: bounding box
    """
    in_data = gdal.Open(filename, 0)
    geotransform = in_data.GetGeoTransform()
    nx = in_data.RasterXSize
    ny = in_data.RasterYSize
    return geotransform2bb(geotransform, nx, ny)


def geotransform2bb(geotransform, nx, ny):
    """
    Transforms the geotiff geotransform into a bounding box
    :param geotransform: (left, cellsize, 0.0, top, 0.0, -cellsize)
    :param nx: nr of x cells
    :param ny: nr of y cells
    :return: bounding box object
    """
    res = geotransform[1]
    x_min = geotransform[0] - res/2
    y_max = geotransform[3] + res/2
    x_max = x_min + nx*res + res
    y_min = y_max - ny*res - res
    return BoundingBox(x_min, y_min, x_max, y_max)


if __name__ == "__main__":
    img = "/home/niels/code/git.pi/pointfiltering/python/out/dtm.tif"
    inLas = "/home/niels/Downloads/PNOA_2009_Lote3_MUR_602-4184_ORT-CLA-COL.las"
    outLas = "/home/niels/Downloads/PNOA_2009_Lote3_MUR_602-4184_ORT-CLA-COL_clip.las"

    clip(inLas, outLas, img)







