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
#ifdef _OPENMP
#include <omp.h>
#endif
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

    return names;
}

py::dict omp_info(void) {
#ifdef _OPENMP
    auto d1 = py::dict(py::arg("threads_max") =
            omp_get_max_threads());
    d1 = py::dict(py::arg("threads_limit") =
            omp_get_thread_limit(), **d1);
    d1 = py::dict(py::arg("num_procs") =
            omp_get_num_procs(), **d1);
    d1 = py::dict(py::arg("dynamic") =
            omp_get_dynamic(), **d1);
    return d1;
#else
    return pybind11::cast<pybind11::none>(Py_None);
#endif
}

py::tuple find_photons(py::array_t<uint16_t> images,
                       uint16_t threshold, int box,
                       int pixel_photon, int pixel_bgnd,
                       int overlap_max, double sum_min, double sum_max,
                       const std::string &return_pixels, bool return_map) {
    py::list out_list;

    py::buffer_info images_buffer = images.request();

    if (images_buffer.ndim != 3) {
        throw std::runtime_error("Number of dimensions must be 3");
    }

    uint16_t *images_ptr = reinterpret_cast<uint16_t*>(images_buffer.ptr);

    PhotonTable<double>* photon_table(new PhotonTable<double>);
    std::vector<uint16_t>* pixels = NULL;

    centroid_params<uint16_t, double> params;
    centroids_initialize_params<uint16_t, double>(&params);

    params.threshold = threshold;
    params.box = box;
    params.pixel_photon_num = pixel_photon;
    params.pixel_bgnd_num = pixel_bgnd;
    params.overlap_max = overlap_max;
    params.sum_min = sum_min;
    params.sum_max = sum_max;
    params.fit_pixels = CENTROIDS_FIT_LMMIN;
    params.x = images_buffer.shape[2];
    params.y = images_buffer.shape[1];
    params.n = images_buffer.shape[0];
    params.return_map = return_map;

    if (!return_pixels.compare("sorted")) {
        params.return_pixels = CENTROIDS_STORE_SORTED;
    } else if (!return_pixels.compare("unsorted")) {
        params.return_pixels = CENTROIDS_STORE_UNSORTED;
    } else if (!return_pixels.compare("none")) {
        params.return_pixels = CENTROIDS_STORE_NONE;
    } else {
        throw std::invalid_argument("Invalid pixel store option");
    }

    if (params.return_pixels != CENTROIDS_STORE_NONE) {
        pixels = new std::vector<uint16_t>;
    }


    if (centroids_calculate_params<uint16_t>(&params)
            != CENTROIDS_PARAMS_OK) {
        throw std::invalid_argument("Invalid parameter combination");
    }

    // Setup our array if needed
    uint16_t *out_ptr = NULL;
    py::array_t<uint16_t> result;
    if (return_map) {
        result = py::array_t<uint16_t>(images_buffer.size);
        py::buffer_info buf2 = result.request();
        out_ptr = reinterpret_cast<uint16_t*>(buf2.ptr);
    }

    pybind11::gil_scoped_release release;

    size_t nphotons = centroids_process<uint16_t, double>(
            images_ptr, out_ptr, photon_table, pixels, params);

    pybind11::gil_scoped_acquire acquire;

    size_t photon_table_cols = CENTROIDS_TABLE_COLS;

    if (params.fit_pixels != CENTROIDS_FIT_NONE) {
        photon_table_cols += 2 * CENTROIDS_FIT_PARAMS_N;
        photon_table_cols += CENTROIDS_FIT_EXTRA_N;
    }

    fprintf(stderr, "Finished centroids_process\n");

    // The following is some jiggery-pokery so we dont
    // have to copy the vector....
    auto capsule = py::capsule(photon_table, [](void *(photon_table))
            { delete reinterpret_cast<std::vector<double>*>((photon_table)); });
    auto table = py::array({static_cast<int>(nphotons),
                            static_cast<int>(photon_table_cols)},
                            photon_table->data(), capsule);
    out_list.append(table);

    // Return also column names
    std::vector<std::string>names = get_column_names(params);
    out_list.append(names);

    // Reshape the output array...
    if (return_map) {
        result.resize({params.n, params.y, params.x});
        out_list.append(result);
    } else {
        out_list.append(pybind11::cast<pybind11::none>(Py_None));
    }

    if (params.return_pixels != CENTROIDS_STORE_NONE) {
        auto pixels_capsule = py::capsule(pixels, [](void *(pixels))
                { delete reinterpret_cast<std::vector<uint16_t>*>((pixels)); });
        auto pixels_array = py::array({static_cast<int>(nphotons),
                static_cast<int>(params.box_t)},
                pixels->data(), pixels_capsule);

        out_list.append(pixels_array);
    } else {
        out_list.append(pybind11::cast<pybind11::none>(Py_None));
    }

    py::tuple out(out_list);
    return out;
}

PYBIND11_MODULE(_pycentroids, m) {
     m.doc() = "Fast centroiding routines for CCD detectors";

     m.def("find_photons", &find_photons,
           "Find photons",
           py::arg("images"),
           py::arg("threshold"),
           py::arg("box"),
           py::arg("pixel_photon"),
           py::arg("pixel_bgnd"),
           py::arg("overlap_max"),
           py::arg("sum_min"),
           py::arg("sum_max"),
           py::arg("return_pixels"),
           py::arg("return_map"));

     m.def("omp_info", &omp_info,
             "Return OpenMP info");

     m.attr("__version__") = CENTROIDS_GIT_VERSION;
}
