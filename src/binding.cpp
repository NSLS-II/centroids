#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <memory>

#include "photons.h"
#include "version.h"

namespace py = pybind11;

py::tuple _find_photons(py::array_t<uint16_t> images, uint16_t threshold, int box)
{
    centroid_params<uint16_t,double> params;
    centroids_initialize_params<uint16_t,double>(&params);

    params.box = box;
    params.threshold = threshold;
    params.pixel_photon_num = 10;
    params.overlap_max = 0;
    params.sum_min = 800;
    params.sum_max = 1210;

    centroids_calculate_params<uint16_t>(&params);

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

    PhotonTablePtr<double> photon_table(new PhotonTable<double>);

    size_t nphotons = centroids_process<uint16_t, double>(in_ptr, out_ptr, photon_table, X, Y, N, params);

    /* Reshape result to have same shape as input */
    result.resize({N, Y, X});
    
    // The following is some jiggery-pokery so we dont have to copy the vector.... 
    auto capsule = py::capsule(photon_table, [](void *(photon_table)) 
            { delete reinterpret_cast<std::vector<double>*>((photon_table)); });
    auto table = py::array(photon_table->size(), photon_table->data(), capsule);
    table.resize({(int)nphotons, 8 + params.box_t});

    py::tuple args = py::make_tuple(table, result);
    return args;
}

PYBIND11_MODULE(pycentroids, m) {
   m.doc() = "Fast centroiding routines for CCD detectors";
   m.def("find_photons", &_find_photons, "Find photons", 
           py::arg("images"), py::arg("threshold"), py::arg("box") = 2);
   m.attr("__version__") = GIT_VERSION;
}
