import pytest
import numpy as np
from scipy.special import erf
# from numpy.testing import assert_array_equal, assert_array_almost_equal
from pycentroids import find_photons


@pytest.fixture
def dataframe():
    def _dataframe(size, offset, mag):
        data = np.random.rand(*size)
        data *= mag
        data += offset
        return data

    return _dataframe


@pytest.fixture
def gauss():
    def _gauss(box, x, y, amp, sigma):
        lin = np.linspace(-box, box, 2 * box + 1)
        xx, yy = np.meshgrid(lin, lin)

        xmin = (xx - 0.5 - x) / (sigma * np.sqrt(2))
        xmax = (xx + 0.5 - x) / (sigma * np.sqrt(2))
        ymin = (yy - 0.5 - y) / (sigma * np.sqrt(2))
        ymax = (yy + 0.5 - y) / (sigma * np.sqrt(2))

        out = erf(xmax) - erf(xmin)
        out *= erf(ymax) - erf(ymin)

        out = out / out.max()
        out *= amp
        return out

    return _gauss


def test_null(dataframe):
    data = dataframe((1, 1400, 1200), 150, 10)
    table, grid, photons = find_photons(data.astype(np.uint16),
                                        threshold=250, box=2)

    assert len(table) == 0


def test_find_photons(dataframe, gauss):
    x = 17
    y = 20
    cen_x = 0.4
    cen_y = 0.15
    sigma = 0.46
    bgnd = 150
    box = 3
    pixel_photon = 9
    pixel_bgnd = 12

    photon = gauss(box, cen_x, cen_y, 500, sigma)
    data = dataframe((1, 1400, 1200), bgnd, 5)

    data[0, y - box:y + box + 1, x - box:x + box + 1] += photon
    data = data.astype(np.uint16)

    int_photon = data[0, y - box:y + box + 1, x - box:x + box + 1]
    photon_sorted = np.sort(int_photon.ravel())[::-1]
    photon_bgnd = photon_sorted[pixel_bgnd:].mean()
    photon_int = (photon_sorted[:pixel_photon] - photon_bgnd).sum()
    photon_bgnd = photon_sorted[pixel_bgnd:].mean()

    table, grid, photons = find_photons(data.astype(np.uint16),
                                        threshold=250, box=box,
                                        sum_min=800, sum_max=1400,
                                        pixel_photon=pixel_photon,
                                        pixel_bgnd=pixel_bgnd)

    assert len(table) == 1
    assert table['Pixel X'][0] == x
    assert table['Pixel Y'][0] == y
    assert pytest.approx(table['Fit X'][0], 0.01) == x + cen_x
    assert pytest.approx(table['Fit Y'][0], 0.01) == y + cen_y
    assert pytest.approx(table['COM X'][0], 0.01) == x + cen_x
    assert pytest.approx(table['COM Y'][0], 0.01) == y + cen_y
    assert pytest.approx(table['Fit Sigma'][0], 0.01) == sigma
    assert pytest.approx(table['Bgnd'][0]) == photon_bgnd
    assert pytest.approx(table['Int'][0]) == photon_int
