#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "photons.h"
#include "version.h"

int _debug_print_flag = 1;

namespace py = pybind11;

py::tuple _find_photons(py::array_t<uint16_t> images, uint16_t threshold, int box)
{
    centroid_params<uint16_t> params;
    centroids_initialize_params<uint16_t>(params);

    params.box = box;
    params.threshold = threshold;
    params.pixel_photon_num = 10;
    params.overlap_max = 0;
    params.sum_min = 800;
    params.sum_max = 1200;

    centroids_calculate_params<uint16_t>(params);

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
    fprintf(stderr, "nphotons = %ld\n", nphotons);

    py::array_t<double> table = py::array_t<double>(8 * nphotons);
    py::buffer_info buf3 = table.request();
    double *table_ptr = (double*)buf3.ptr;

    // This is not so good, we now do a copy into the memory`
    for(auto table = photon_table->begin(); 
        table != photon_table->end(); 
        table++)
    {
        *(table_ptr++) = table->x;
        *(table_ptr++) = table->y;
        *(table_ptr++) = table->com_x;
        *(table_ptr++) = table->com_y;
        *(table_ptr++) = table->sum;
        *(table_ptr++) = table->bgnd;
        *(table_ptr++) = table->box_sum;
        *(table_ptr++) = 0;
    }

    /* Reshape result to have same shape as input */
    result.resize({N, Y, X});
    table.resize({(int)nphotons, 8});

    py::tuple args = py::make_tuple(table, result);
    return args;
}

PYBIND11_MODULE(pycentroids, m) {
   m.doc() = "Fast centroiding routines for CCD detectors";
   m.def("find_photons", &_find_photons, "Find photons", 
           py::arg("images"), py::arg("threshold"), py::arg("box") = 2);
   m.attr("__version__") = GIT_VERSION;
}
