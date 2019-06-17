import pandas as pd
from _pycentroids import find_photons as _find_photons


def find_photons(images, threshold=200, box=2, pixel_photon=10, overlap_max=0,
                 sum_min=800, sum_max=1250, return_pixels='none',
                 return_map=False):
    _rtn = _find_photons(images=images, threshold=threshold, box=box,
                         pixel_photon=pixel_photon, overlap_max=overlap_max,
                         sum_min=sum_min, sum_max=sum_max,
                         return_pixels=return_pixels, return_map=return_map)

    df = pd.DataFrame(_rtn[0], columns=_rtn[1])

    rtn = (df, ) + _rtn[2:]
    return rtn;
