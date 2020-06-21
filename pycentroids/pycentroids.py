import pandas as pd
import numpy as np
from _pycentroids import find_photons as _find_photons


def find_photons(images, filter=None,
                 threshold=200, box=2, search_box=1, pixel_photon=9,
                 pixel_bgnd=15, com_photon=9, overlap_max=0, sum_min=800,
                 sum_max=1250, pixel_lut=None, pixel_lut_range=None,
                 fit_pixels_2d=True, fit_pixels_1d_x=True,
                 fit_pixels_1d_y=True, fit_constraints=None,
                 fit_weights_2d=None, fit_weights_1d=None,
                 return_pixels='none', return_map=False, tag_pixels=False):
    """Find photons in CCD images and process for sub-pixel center.

    Parameters
    ----------
    images : numpy.ndarray
        Images to process of 3 dimensions
    filter : numpy.ndarray or None
        Filter image of 2d or 3d shape. Should be non-zero to filter pixels.
        If none, then no pixels are filtered.
    threshold : int
        Pixel value threshold for first search of images
    box : int
        Value for pixel box, Box size is
        (2 * box + 1) * (2 * box + 1)
    search_box : int
        Value for pixel seaarch box, Box size is
        (2 * box + 1) * (2 * box + 1)
        This value defines the box which is used to ensure that the
        center pixel is the highest. This should be set to the size of
        the charge cloud from the photon event.
    pixel_photon : int
        Number of pixels to include for photon intensity. After sorting
        from high to low values the first pixel_photon values are summed
        to get the intensity values.
    pixel_background : int
        Number of pixels to include for background determination. After sorting
        from high to low values the values from pixel_background to the end of
        the array are averaged to get the background value.
    com_photon : int
        Number of pixels to include for COM calculation. After sorting
        from high to low values the first com_photon values are used
        to determine the photons COM.
    overlap_max : int
        The maximum overlap to accept before rejecting the photon
        from the table.
    sum_min : float
        The minimum integrated intensity to filter the output table.
    sum_max : float
        The maximum integrated intensity to filter the output table.
    pixel_lut : np.array
        Lookup table for the pixel COM correction. 1D array of
        corrected position
    pixel_lut_range: tuple
        Tuple of the start and end coordinates for the LUT
    fit_pixels_2d : bool
        If true, fit the pixels from a photon with a 2D gaussian
    fit_pixels_1dx : bool
        If true, fit the pixels from a photon with a 1D gaussian integrating
        along the y-axis to get the x position
    fit_pixels_1dy : bool
        If true, fit the pixels from a photon with a 1D gaussian integrating
        along the x-axis to get the y position
    fit_constraints : dict
        Dictionary of constraints to use for fitting pixels
        Currently the dictionary keys can be:
            pos_range, pos_cent : Position (x) of photon range and center
            sigma_range, sigma_cent : Sigma of the gaussian range anc center
    fit_weights_2d : np.array
        array of size (2 * box + 1) * (2 * box + 1) of the weights to use for
        the pixel fit
    fit_weights_1d : np.array
        array of size (2 * box + 1) of the weights to use for the pixel fits
    return_pixels : 'none', 'sorted', 'unsorted'
        Option to return array of pixel values
    return_map : bool
        Option to return map of located photons.
    tag_pixels : bool
        Option to tag the found photons in the pixel map by setting the MSB

    Returns
    -------
    photon_table : Pandas DataFrame
        Table of photons and their parameters
    photon_map : nummpy.ndarray
        Map of located photons
    pixel_values : numpy.ndarray
        Pixel values from located photons

    """

    if filter is None:
        filter = np.zeros_like(images[0])

    if fit_constraints is None:
        fit_constraints = {}

    if fit_weights_1d is None:
        fit_weights_1d = np.ones((2 * box + 1))

    if fit_weights_2d is None:
        fit_weights_2d = np.ones((2 * box + 1, 2 * box + 1))

    if pixel_lut_range is None:
        pixel_lut_range = (-1, 1)

    if pixel_lut is None:
        pixel_lut = np.linspace(pixel_lut_range[0], pixel_lut_range[1], 2000)

    _rtn = _find_photons(images=images, filter=filter,
                         threshold=threshold, box=box,
                         search_box=search_box,
                         pixel_photon=pixel_photon, pixel_bgnd=pixel_bgnd,
                         com_photon=com_photon, overlap_max=overlap_max,
                         sum_min=sum_min, sum_max=sum_max,
                         pixel_lut=pixel_lut,
                         pixel_lut_range=pixel_lut_range,
                         fit_pixels_2d=fit_pixels_2d,
                         fit_pixels_1dx=fit_pixels_1d_x,
                         fit_pixels_1dy=fit_pixels_1d_y,
                         fit_constraints=fit_constraints,
                         fit_weights_2d=fit_weights_2d,
                         fit_weights_1d=fit_weights_1d,
                         return_pixels=return_pixels, return_map=return_map,
                         tag_pixels=tag_pixels)

    df = pd.DataFrame(_rtn[0], columns=_rtn[1])

    rtn = (df, ) + _rtn[2:]
    return rtn
