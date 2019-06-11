#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <memory>

#include "photons.h"
#include "version.h"

namespace py = pybind11;

py::tuple _find_photons(py::array_t<uint16_t> images, 
                        uint16_t threshold, int box, int pixel_photon,
                        int overlap_max, double sum_min, double sum_max)
{
    /* read input arrays buffer_info */
    py::buffer_info buf1 = images.request();

    /* allocate the output buffer */
    py::array_t<uint16_t> result = py::array_t<uint16_t>(buf1.size);
    py::buffer_info buf2 = result.request();

    // Allcoate the table buffer

    uint16_t *in_ptr = (uint16_t*)buf1.ptr;
    uint16_t *out_ptr = (uint16_t*)buf2.ptr;
    size_t X = buf1.shape[2];
    size_t Y = buf1.shape[1];
    size_t N = buf1.shape[0];

    centroid_params<uint16_t,double> params;
    centroids_initialize_params<uint16_t,double>(&params);

    params.threshold = threshold;
    params.box = box;
    params.pixel_photon_num = pixel_photon;
    params.overlap_max = overlap_max;
    params.sum_min = sum_min;
    params.sum_max = sum_max;
    params.store_pixels = false;
    params.x = buf1.shape[2];
    params.y = buf1.shape[1];
    params.n = buf1.shape[0];

    centroids_calculate_params<uint16_t>(&params);

    PhotonTablePtr<double> photon_table(new PhotonTable<double>);
    size_t nphotons = centroids_process<uint16_t, double>(in_ptr, out_ptr, photon_table, X, Y, N, params);

    size_t photon_table_cols = 9;
    if(params.store_pixels)
    {
        photon_table_cols += params.box_t;
    }

    // The following is some jiggery-pokery so we dont have to copy the vector.... 
    // https://stackoverflow.com/questions/54876346/pybind11-and-stdvector-how-to-free-data-using-capsules
    auto capsule = py::capsule(photon_table, [](void *(photon_table)) 
            { delete reinterpret_cast<std::vector<double>*>((photon_table)); });
    auto table = py::array({(int)nphotons, (int)photon_table_cols}, photon_table->data(), capsule);

    // Reshape the output array...
    result.resize({N, Y, X});

    py::tuple args = py::make_tuple(table, result);
    return args;
}

PYBIND11_MODULE(pycentroids, m) {
   m.doc() = "Fast centroiding routines for CCD detectors";
   m.def("find_photons", &_find_photons, 
         "Find photons", 
         py::arg("images"), 
         py::arg("threshold") = 200, 
         py::arg("box") = 2,
         py::arg("pixel_photon") = 10,
         py::arg("overlap_max") = 0,
         py::arg("sum_min") = 800,
         py::arg("sum_max") = 1250);
   m.attr("__version__") = GIT_VERSION;
}
