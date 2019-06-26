import pandas as pd
from _pycentroids import find_photons as _find_photons


def find_photons(images, threshold=200, box=2, pixel_photon=10, pixel_bgnd=15,
                 overlap_max=0, sum_min=800, sum_max=1250,
                 return_pixels='none', return_map=False):
    _rtn = _find_photons(images=images, threshold=threshold, box=box,
                         pixel_photon=pixel_photon, pixel_bgnd=pixel_bgnd,
                         overlap_max=overlap_max,
                         sum_min=sum_min, sum_max=sum_max,
                         return_pixels=return_pixels, return_map=return_map)
    """Find photons in CCD images and process for sub-pixel center.

    Parameters
    ----------
    images : numpy.ndarray
        Images to process of 3 dimensions
    threshold : int
        Pixel value threshold for first search of images
    box : int
        Value for pixel seaarch box, Box size is
        (2 * box + 1) * (2 * box + 1)
    pixel_photon : int
        Number of pixels to include for photon intensity. After sorting
        from high to low values the first pixel_photon values are summed
        to get the intensity values.
    pixel_background : int
        Number of pixels to include for background determination. After sorting
        from high to low values the values from pixel_background to the end of
        the array are averaged to get the background value.
    overlap_max : int
        The maximum overlap to accept before rejecting the photon
        from the table.
    sum_min : float
        The minimum integrated intensity to filter the output table.
    sum_max : float
        The maximum integrated intensity to filter the output table.
    return_pixels : 'none', 'sorted', 'unsorted'
        Option to return array of pixel values
    return_map : bool
        Option to return map of located photons.

    Returns
    -------
    photon_table : Pandas DataFrame
        Table of photons and their parameters
    photon_map : nummpy.ndarray
        Map of located photons
    pixel_values : numpy.ndarray
        Pixel values from located photons

    """

    df = pd.DataFrame(_rtn[0], columns=_rtn[1])

    rtn = (df, ) + _rtn[2:]
    return rtn
