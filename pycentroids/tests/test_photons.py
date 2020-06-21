import pytest
import numpy as np
import pandas as pd
from scipy.special import erf
from numpy.testing import assert_array_equal
from pycentroids import find_photons
from packaging import version

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

if version.parse(pd.__version__) < version.parse("1.0"):
    pd.set_option('display.max_colwidth', -1)
else:
    pd.set_option('display.max_colwidth', None)


@pytest.fixture
def dataframe():
    def _dataframe(size, offset, sigma):
        data = np.random.normal(offset, sigma, np.product(size))
        data = data.reshape(size)
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
    data = dataframe((1, 1400, 1200), 150, 1)
    table, grid, photons = find_photons(data.astype(np.uint16),
                                        None,
                                        threshold=250, box=2)

    assert len(table) == 0


def test_find_photons(dataframe, gauss):
    image_x = 1400
    image_y = 1200
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
    data = dataframe((1, image_x, image_y), bgnd, 1)

    data[0, y - box:y + box + 1, x - box:x + box + 1] += photon
    data = data.astype(np.uint16)

    int_photon = data[0, y - box:y + box + 1, x - box:x + box + 1]
    photon_sorted = np.sort(int_photon.ravel())[::-1]
    photon_bgnd = photon_sorted[pixel_bgnd:].mean()
    photon_int = (photon_sorted[:pixel_photon] - photon_bgnd).sum()
    photon_bgnd = photon_sorted[pixel_bgnd:].mean()

    fit_constraints = {
        "pos_range": 0.5,
        "pos_cent": 0,
        "sigma_range": 0.2,
        "sigma_cent": 0.45
    }

    table, grid, photons = find_photons(data.astype(np.uint16),
                                        None,
                                        threshold=250, box=box,
                                        search_box=box,
                                        sum_min=800, sum_max=1400,
                                        pixel_lut=np.linspace(-2, 2, 2000),
                                        pixel_lut_range=(-1, 1),
                                        pixel_photon=pixel_photon,
                                        pixel_bgnd=pixel_bgnd,
                                        fit_constraints=fit_constraints,
                                        return_map=True,
                                        return_pixels='unsorted')

    assert len(table) == 1
    assert photons.shape == (1, 2 * box + 1, 2 * box + 1)
    assert grid.shape == (1, image_x, image_y)

    # Check returned photons
    assert_array_equal(photons[0],
                       data[0, y - box:y + box + 1, x - box:x + box + 1])

    # Check photon mask
    mask = np.zeros_like(grid)
    mask[0, y - box:y + box + 1, x - box:x + box + 1] = 1
    assert_array_equal(mask, grid)

    # Check pixel values and fit
    print(table)
    assert table['Pixel X'][0] == x
    assert table['Pixel Y'][0] == y
    assert pytest.approx(table['COM X'][0], 0.01) == x + cen_x
    assert pytest.approx(table['COM Y'][0], 0.01) == y + cen_y
    assert pytest.approx(table['COR COM X'][0], 0.01) == x + 2 * cen_x
    assert pytest.approx(table['COR COM Y'][0], 0.01) == y + 2 * cen_y
    assert pytest.approx(table['Int'][0]) == photon_int
    assert pytest.approx(table['Bgnd'][0]) == photon_bgnd
    assert pytest.approx(table['Fit X'][0], 0.05) == x + cen_x
    assert pytest.approx(table['Fit Y'][0], 0.05) == y + cen_y
    assert pytest.approx(table['Fit Sigma'][0], 0.1) == sigma
    assert pytest.approx(table['Fit 1DX X'][0], 0.05) == x + cen_x
    assert pytest.approx(table['Fit 1DY Y'][0], 0.05) == y + cen_y
    assert pytest.approx(table['Fit 1DX Sigma'][0], 0.1) == sigma
    assert pytest.approx(table['Fit 1DY Sigma'][0], 0.1) == sigma


def test_mask(dataframe, gauss):
    image_x = 1400
    image_y = 1200
    x = np.array([17, 50, 89, 400, 17, 762, 90, 1100])
    y = np.array([20, 100, 400, 753, 20, 1000, 41, 200])
    cen_x = 0.4
    cen_y = 0.15
    sigma = 0.46
    bgnd = 150
    box = 3
    pixel_photon = 9
    pixel_bgnd = 12

    photon = gauss(box, cen_x, cen_y, 500, sigma)
    data = dataframe((2, image_x, image_y), bgnd, 5)

    for _x, _y in zip(x[:4], y[:4]):
        data[0, _y - box:_y + box + 1, _x - box:_x + box + 1] += photon
    for _x, _y in zip(x[4:], y[4:]):
        data[1, _y - box:_y + box + 1, _x - box:_x + box + 1] += photon

    data = data.astype(np.uint16)
    mask = np.zeros_like(data)

    table, grid, photons = find_photons(data,
                                        mask,
                                        threshold=250, box=box,
                                        search_box=box,
                                        sum_min=800, sum_max=1400,
                                        pixel_photon=pixel_photon,
                                        pixel_bgnd=pixel_bgnd,
                                        return_map=True,
                                        return_pixels='unsorted')

    assert len(table) == len(x)
    for _x, _y in zip(x, y):
        tx = table['Pixel X'] == _x
        ty = table['Pixel Y'] == _y
        assert (tx & ty).any()

    table, grid, photons = find_photons(data, mask[0],
                                        threshold=250, box=box,
                                        search_box=box,
                                        sum_min=800, sum_max=1400,
                                        pixel_photon=pixel_photon,
                                        pixel_bgnd=pixel_bgnd,
                                        return_map=True,
                                        return_pixels='unsorted')

    assert len(table) == len(x)
    for _x, _y in zip(x, y):
        tx = table['Pixel X'] == _x
        ty = table['Pixel Y'] == _y
        assert (tx & ty).any()

    test_mask = np.copy(mask[0])
    test_mask[y[0], x[0]] = 1

    table, grid, photons = find_photons(data, test_mask,
                                        threshold=250, box=box,
                                        search_box=box,
                                        sum_min=800, sum_max=1400,
                                        pixel_photon=pixel_photon,
                                        pixel_bgnd=pixel_bgnd,
                                        return_map=True,
                                        return_pixels='unsorted')

    assert len(table) == (len(x) - 2)
    for _x, _y in zip(np.concatenate((x[1:3], x[5:])),
                      np.concatenate((y[1:3], y[5:]))):
        tx = table['Pixel X'] == _x
        ty = table['Pixel Y'] == _y
        assert (tx & ty).any()

    # Check returned mask has MSB set
    assert_array_equal(grid[:, y[0], x[0]], np.ones(data.shape[0]) * 0x8000)

    test_mask = np.copy(mask)
    test_mask[0, y[0], x[0]] = 1

    table, grid, photons = find_photons(data, test_mask,
                                        threshold=250, box=box,
                                        search_box=box,
                                        sum_min=800, sum_max=1400,
                                        pixel_photon=pixel_photon,
                                        pixel_bgnd=pixel_bgnd,
                                        return_map=True,
                                        return_pixels='unsorted')

    assert len(table) == (len(x) - 1)
    for _x, _y in zip(x[1:], y[1:]):
        tx = table['Pixel X'] == _x
        ty = table['Pixel Y'] == _y
        assert (tx & ty).any()
