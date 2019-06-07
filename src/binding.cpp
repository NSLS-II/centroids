#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "photons.h"
#include "version.h"

int _debug_print_flag = 1;

namespace py = pybind11;

py::tuple _find_photons(py::array_t<uint16_t> images, uint16_t threshold, int box)
{
    centroid_params<uint16_t> params;
    params.box = box;
    params.threshold = threshold;

    centroids_initialize_params<uint16_t>(params);

    /* read input arrays buffer_info */
    py::buffer_info buf1 = images.request();

    /* allocate the output buffer */
    py::array_t<uint16_t> result = py::array_t<uint16_t>(buf1.size);
    py::buffer_info buf2 = result.request();

    // Allcoate the table buffer
   
    int box_total = (2 * box + 1) * (2 * box + 1);
    py::array_t<double> table = py::array_t<double>(100000 * (box_total + 7));
    py::buffer_info buf3 = table.request();

    /* allocate the bias buffer */
    py::array_t<double> bias = py::array_t<uint16_t>(buf1.size);
    py::buffer_info buf4 = bias.request();

    uint16_t *in_ptr = (uint16_t*)buf1.ptr;
    uint16_t *out_ptr = (uint16_t*)buf2.ptr;
    double *table_ptr = (double*)buf3.ptr;
    double *bias_ptr = (double*)buf4.ptr;
    size_t X = buf1.shape[1];
    size_t Y = buf1.shape[0];

    // Blank the table

    for(int i=0;i<100000*(7 + box_total);i++)
    {
        table_ptr[i] = 0.0;
    }

    int nphotons = process_image<uint16_t>(in_ptr, out_ptr, table_ptr, bias_ptr, X, Y, params);

    /* Reshape result to have same shape as input */
    result.resize({Y, X});
    bias.resize({Y, X});
    table.resize({nphotons, 7 + box_total});

    py::tuple args = py::make_tuple(table, result, bias);
    return args;
}

PYBIND11_MODULE(pycentroids, m) {
   m.doc() = "Fast centroiding routines for CCD detectors";
   m.def("find_photons", &_find_photons, "Find photons", 
           py::arg("images"), py::arg("threshold"), py::arg("box") = 1);
   m.attr("__version__") = GIT_VERSION;
}
