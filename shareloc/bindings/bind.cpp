/*
!/usr/bin/env python
coding: utf8

Copyright (c) 2023 Centre National d'Etudes Spatiales (CNES).

This file is part of shareloc
(see https://github.com/CNES/shareloc).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/**
Thes purpose of this module is only to "bind" cpp code.
It gives to the compiler the instructions to compile the usefull cpp code into an .so file
which is callable in a python code as a python module.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "rpc.cpp"

namespace py = pybind11;


PYBIND11_MODULE(rpc_c, m) {

    py::class_<GeoModelTemplate>(m, "GeoModelTemplate")
        .def(py::init<>())
        .def("direct_loc_h", &GeoModelTemplate::direct_loc_h)
        .def("direct_loc_dtm", &GeoModelTemplate::direct_loc_dtm)
        .def("inverse_loc", &GeoModelTemplate::inverse_loc);

    py::class_<RPC,GeoModelTemplate>(m, "RPC")
        .def(py::init<>())
        .def("direct_loc_h", &RPC::direct_loc_h)
        .def("direct_loc_grid_h", &RPC::direct_loc_grid_h)
        .def("direct_loc_dtm", &RPC::direct_loc_dtm)
        .def("filter_coordinates", &RPC::filter_coordinates)
        .def("compute_loc_inverse_derivates", &RPC::compute_loc_inverse_derivates)
        .def("direct_loc_inverse_iterative", &RPC::direct_loc_inverse_iterative)
        .def("get_alt_min_max", &RPC::get_alt_min_max)
        .def("los_extrema", &RPC::los_extrema);

    //m.doc() = "Pybind hello world"; // optional module docstring
    m.def("polynomial_equation", &polynomial_equation, "TODO: doc");
    m.def("compute_rational_function_polynomial", &compute_rational_function_polynomial,
        "TODO: doc");
    m.def("derivative_polynomial_latitude", &derivative_polynomial_latitude, "TODO: doc");
    m.def("derivative_polynomial_longitude", &derivative_polynomial_longitude, "TODO: doc");
    m.def("compute_loc_inverse_derivates_numba", &compute_loc_inverse_derivates_numba, "TODO: doc");
}

//c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes)
//bind.cpp -o pbrpc$(python3-config --extension-suffix)