""" dummy functions will be here"""
import numpy as np

def rescale(x, scale_min, scale_max):
    """ Using the law of re-scalling this function will take an array and
    scale the minimum and maximum to a specified limit.

    Input
    ----
    x (array): data you want to re-scale
    scale_min(float): minimum of your scaled array
    scale_max(float): maximum of your scaled array

    Output
    ------
    x_scale (array): the re-scaled array
    """
    return (((scale_max-scale_min)*(x-np.max(x)))/(np.max(x) - np.min(x))) + scale_max
