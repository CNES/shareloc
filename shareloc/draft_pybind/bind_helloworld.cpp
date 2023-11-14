#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "rpc.cpp"

namespace py = pybind11;


PYBIND11_MODULE(pbrpc, m) {

    py::class_<GeoModelTemplate>(m, "GeoModelTemplate")
        .def(py::init<string>())
        .def("direct_loc_h", &GeoModelTemplate::direct_loc_h)
        .def("direct_loc_dtm", &GeoModelTemplate::direct_loc_dtm)
        .def("inverse_loc", &GeoModelTemplate::inverse_loc);

    py::class_<RPC,GeoModelTemplate>(m, "RPC")
        .def(py::init<string>())
        .def("direct_loc_h", &RPC::direct_loc_h)
        .def("direct_loc_grid_h", &RPC::direct_loc_grid_h)
        .def("direct_loc_dtm", &RPC::direct_loc_dtm)
        .def("filter_coordinates", &RPC::filter_coordinates)
        .def("compute_loc_inverse_derivates", &RPC::compute_loc_inverse_derivates)
        .def("direct_loc_inverse_iterative", &RPC::direct_loc_inverse_iterative)
        .def("get_alt_min_max", &RPC::get_alt_min_max)
        .def("los_extrema", &RPC::los_extrema);

    //m.doc() = "Pybind hello world"; // optional module docstring
    m.def("parse_coeff_line", &parse_coeff_line, "TODO: doc");
    m.def("polynomial_equation", &polynomial_equation, "TODO: doc");
    m.def("compute_rational_function_polynomial", &compute_rational_function_polynomial, "TODO: doc");
    m.def("derivative_polynomial_latitude", &derivative_polynomial_latitude, "TODO: doc");
    m.def("derivative_polynomial_longitude", &derivative_polynomial_longitude, "TODO: doc");
    m.def("compute_loc_inverse_derivates_numba", &compute_loc_inverse_derivates_numba, "TODO: doc");
}

//c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) bind_helloworld.cpp -o pbrpc$(python3-config --extension-suffix)