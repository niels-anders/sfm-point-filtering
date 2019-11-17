#!/usr/bin/env python
"""
Set of vegetation indices
"""

import numpy


def exg(points):
    """
    Excessive Greenness (ExG)
    """
    [red, green, blue] = convert_colors(points.red, points.green, points.blue)
    normed_r = red / (red+green+blue)
    normed_g = green / (red+green+blue)
    normed_b = blue / (red+green+blue)
    return normed_g*2 - normed_r - normed_b


def vvi(points):
    """
    Visual Vegetation Index (VVI)
    """
    rgb0 = (40, 60, 10)
    [red, green, blue] = convert_colors(points.red, points.green, points.blue)
    weight = 1
    normed_r = 1-numpy.abs((red-rgb0[0])/(red+rgb0[0]))
    normed_g = 1-numpy.abs((green-rgb0[1])/(green+rgb0[1]))
    normed_b = 1-numpy.abs((blue-rgb0[2])/(blue+rgb0[2]))
    return (normed_r * normed_g * normed_b)**(1/weight)


def convert_colors(red16bit, green16bit, blue16bit):
    """
    Convert color values from 16 bit signed integers
    to 8 bit unsigned integers
    """
    # from 16 bit signed integers to [-0.5 - 0.5]
    red8bit = red16bit/float(2**16)
    green8bit = green16bit/float(2**16)
    blue8bit = blue16bit/float(2**16)
    # from signed values to unsigned values
    red8bit[red8bit < 0] = red8bit[red8bit < 0] + 1
    green8bit[green8bit < 0] = green8bit[green8bit < 0] + 1
    blue8bit[blue8bit < 0] = blue8bit[blue8bit < 0] + 1
    return red8bit, green8bit, blue8bit
