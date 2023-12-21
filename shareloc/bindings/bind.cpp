/*
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
        .def(py::init<array<double, 20>,
        array<double, 20>,
        array<double, 20>,
        array<double, 20>,
        array<double, 10>>())
        .def("direct_loc_h", &RPC::direct_loc_h)
        .def("direct_loc_grid_h", &RPC::direct_loc_grid_h)
        .def("direct_loc_dtm", &RPC::direct_loc_dtm)
        .def("inverse_loc", &GeoModelTemplate::inverse_loc)
        .def("inverse_loc",\
        (tuple<double,double,double> (RPC::*)\
        (double,double,double)) &RPC::inverse_loc)

        .def("filter_coordinates", &RPC::filter_coordinates)
        .def("compute_loc_inverse_derivates",\
        (tuple<vector<double>,vector<double>,vector<double>,vector<double>> (RPC::*)\
        (vector<double>,vector<double>,vector<double>)) &RPC::compute_loc_inverse_derivates)

        .def("compute_loc_inverse_derivates",\
        (tuple<double,double,double,double> (RPC::*)\
        (double,double,double)) &RPC::compute_loc_inverse_derivates)

        .def("direct_loc_inverse_iterative", &RPC::direct_loc_inverse_iterative)
        .def("get_alt_min_max", &RPC::get_alt_min_max)
        .def("los_extrema", &RPC::los_extrema)
        .def("get_num_col", &RPC::get_num_col)
        .def("get_den_col", &RPC::get_den_col)
        .def("get_num_row", &RPC::get_num_row)
        .def("get_den_row", &RPC::get_den_row)
        .def("get_num_lon", &RPC::get_num_lon)
        .def("get_den_lon", &RPC::get_den_lon)
        .def("get_num_lat", &RPC::get_num_lat)
        .def("get_den_lat", &RPC::get_den_lat)
        .def("get_offset_row", &RPC::get_offset_row)
        .def("get_scale_row", &RPC::get_scale_row)
        .def("get_offset_col", &RPC::get_offset_col)
        .def("get_scale_col", &RPC::get_scale_col)
        .def("get_offset_alt", &RPC::get_offset_alt)
        .def("get_scale_alt", &RPC::get_scale_alt)
        .def("get_offset_lon", &RPC::get_offset_lon)
        .def("get_scale_lon", &RPC::get_scale_lon)
        .def("get_offset_lat", &RPC::get_offset_lat)
        .def("get_scale_lat", &RPC::get_scale_lat);

    //m.doc() = "Pybind hello world"; // optional module docstring
    m.def("polynomial_equation", &polynomial_equation, "Compute polynomial equation");
    m.def("compute_rational_function_polynomial_unitary", \
        &compute_rational_function_polynomial_unitary,
        "Compute rational function polynomial for only one point");
    m.def("compute_rational_function_polynomial", &compute_rational_function_polynomial,
        "Compute rational function polynomial. Useful to compute direct and inverse localization"
        "using direct or inverse RPC.");
    m.def("derivative_polynomial_latitude", &derivative_polynomial_latitude,
    "Compute latitude derivative polynomial equation");
    m.def("derivative_polynomial_longitude", &derivative_polynomial_longitude,
    "Compute longitude derivative polynomial equation");
    m.def("compute_loc_inverse_derivates_optimized_unitary",
    &compute_loc_inverse_derivates_optimized_unitary,
    "Analytically compute the partials derivatives of inverse localization for only one point");
    m.def("compute_loc_inverse_derivates_optimized",
    &compute_loc_inverse_derivates_optimized,
    "Analytically compute the partials derivatives of inverse localization");
}

//c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes)
//bind.cpp -o pbrpc$(python3-config --extension-suffix)