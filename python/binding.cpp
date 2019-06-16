//
// CENTROIDS : C++ implementation of single photon counting for CCDs
// Stuart B. Wilkins, Brookhaven National Laboratory
//
//
// BSD 3-Clause License
//
// Copyright (c) 2019, Brookhaven Science Associates
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
// THE POSSIBILITY OF SUCH DAMAGE.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <memory>

#include "centroids.h"

namespace py = pybind11;

extern const char* CENTROIDS_GIT_VERSION;

template <typename DT, typename OT>
std::vector<std::string> get_column_names(
        const centroid_params<DT, OT> &params) {
    std::vector<std::string> names;

    names.insert(names.end(),
            &centroids_photon_table_names[0],
            &centroids_photon_table_names[CENTROIDS_TABLE_COLS]);

    if (params.fit_pixels != CENTROIDS_FIT_NONE) {
        names.insert(names.end(),
                &centroids_photon_table_names[CENTROIDS_TABLE_COLS],
                &centroids_photon_table_names
                    [CENTROIDS_TABLE_COLS
                     + (2 * CENTROIDS_FIT_PARAMS_N)
                     + CENTROIDS_FIT_EXTRA_N]);
    }

    if (params.store_pixels != CENTROIDS_STORE_NONE) {
        for (int n = 0; n < params.box_t; n++) {
            names.push_back(std::string("Pixel ") + std::to_string(n));
        }
    }

    return names;
}

py::tuple _find_photons(py::array_t<uint16_t> images,
                        uint16_t threshold, int box, int pixel_photon,
                        int overlap_max, double sum_min, double sum_max,
                        const std::string &store_pixels) {
    // Do numpy array request to get size
    py::buffer_info buf1 = images.request();

    centroid_params<uint16_t, double> params;
    centroids_initialize_params<uint16_t, double>(&params);

    params.threshold = threshold;
    params.box = box;
    params.pixel_photon_num = pixel_photon;
    params.overlap_max = overlap_max;
    params.sum_min = sum_min;
    params.sum_max = sum_max;
    params.fit_pixels = CENTROIDS_FIT_LMMIN;
    params.x = buf1.shape[2];
    params.y = buf1.shape[1];
    params.n = buf1.shape[0];

    if (!store_pixels.compare("sorted")) {
        params.store_pixels = CENTROIDS_STORE_SORTED;
    } else if (!store_pixels.compare("unsorted")) {
        params.store_pixels = CENTROIDS_STORE_UNSORTED;
    } else if (!store_pixels.compare("none")) {
        params.store_pixels = CENTROIDS_STORE_NONE;
    } else {
        throw std::invalid_argument("Invalid pixel store option");
    }

    centroids_calculate_params<uint16_t>(&params);

    // Do numpy arrays
    py::array_t<uint16_t> result = py::array_t<uint16_t>(buf1.size);
    py::buffer_info buf2 = result.request();
    uint16_t *in_ptr = reinterpret_cast<uint16_t*>(buf1.ptr);
    uint16_t *out_ptr = reinterpret_cast<uint16_t*>(buf2.ptr);

    // Setup the photon table
    PhotonTable<double>* photon_table(new PhotonTable<double>);

    // Process photons
    size_t nphotons = centroids_process<uint16_t, double>(
            in_ptr, out_ptr, photon_table, params);

    size_t photon_table_cols = CENTROIDS_TABLE_COLS;
    if (params.store_pixels != CENTROIDS_STORE_NONE) {
        photon_table_cols += params.box_t;
    }

    if (params.fit_pixels != CENTROIDS_FIT_NONE) {
        photon_table_cols += 2 * CENTROIDS_FIT_PARAMS_N;
        photon_table_cols += CENTROIDS_FIT_EXTRA_N;
    }

    // The following is some jiggery-pokery so we dont
    // have to copy the vector....
    auto capsule = py::capsule(photon_table, [](void *(photon_table))
            { delete reinterpret_cast<std::vector<double>*>((photon_table)); });
    auto table = py::array({static_cast<int>(nphotons),
                            static_cast<int>(photon_table_cols)},
                            photon_table->data(), capsule);

    // Reshape the output array...
    result.resize({params.n, params.y, params.x});

    // Return also column names
    std::vector<std::string>names = get_column_names(params);
    py::tuple args = py::make_tuple(table, result, names);
    return args;
}

PYBIND11_MODULE(_pycentroids, m) {
     m.doc() = "Fast centroiding routines for CCD detectors";
     m.def("find_photons", &_find_photons,
           "Find photons",
           py::arg("images"),
           py::arg("threshold") = 200,
           py::arg("box") = 2,
           py::arg("pixel_photon") = 10,
           py::arg("overlap_max") = 0,
           py::arg("sum_min") = 800,
           py::arg("sum_max") = 1250,
           py::arg("store_pixels") = "none");
     m.attr("__version__") = CENTROIDS_GIT_VERSION;
}
