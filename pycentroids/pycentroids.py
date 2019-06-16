import pandas as pd
from _pycentroids import find_photons as _find_photons


def find_photons(images, threshold=200, box=2, pixel_photon=10, overlap_max=0,
                 sum_min=800, sum_max=1250, store_pixels='none'):
    rtn = _find_photons(images=images, threshold=threshold, box=box,
                        pixel_photon=pixel_photon, overlap_max=overlap_max,
                        sum_min=sum_min, sum_max=sum_max,
                        store_pixels=store_pixels)

    df = pd.DataFrame(rtn[0], columns=rtn[2])

    return (df, rtn[1])
