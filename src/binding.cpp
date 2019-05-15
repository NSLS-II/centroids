#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "photons.h"

int _debug_print_flag = 1;

namespace py = pybind11;

py::tuple _find_photons(py::array_t<uint16_t> images)
{
    /* read input arrays buffer_info */
    py::buffer_info buf1 = images.request();

    /* allocate the output buffer */
    py::array_t<uint16_t> result = py::array_t<uint16_t>(buf1.size);
    py::buffer_info buf2 = result.request();

    // Allcoate the table buffer
    
    py::array_t<double> table = py::array_t<double>(100000 * 6);
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

    for(int i=0;i<1000*6;i++)
    {
        table_ptr[i] = 0.0;
    }

    find_photons_uint16(in_ptr, out_ptr, table_ptr, bias_ptr, X, Y, 400);  

    /* Reshape result to have same shape as input */
    result.resize({Y, X});
    bias.resize({Y, X});
    table.resize({1000, 6});

    py::tuple args = py::make_tuple(table, result, bias);
    return args;
}

PYBIND11_MODULE(centroids, m) {
   m.doc() = "Fast centroiding routines for CCD detectors";
   m.def("find_photons", &_find_photons, "Find photons");
#ifdef VERSION_INFO
   m.attr("__version__") = VERSION_INFO;
#else
   m.attr("__version__") = "dev";
#endif
}
