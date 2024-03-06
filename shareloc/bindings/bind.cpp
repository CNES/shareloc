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

#include <pybind11/stl.h>

#include "rpc.hpp"
#include "rectification.cpp"

namespace py = pybind11;


PYBIND11_MODULE(bindings_cpp, m) {

    py::class_<DTMIntersection>(m, "DTMIntersection")
        .def(py::init<>())
        .def(py::init<int,py::array_t<double, py::array::c_style | py::array::forcecast> \
                ,int,int,std::tuple<double,double,double,double,double,double>>())
        .def("eq_plan", &DTMIntersection::eq_plan)
        .def("ter_to_index", &DTMIntersection::ter_to_index)
        .def("index_to_ter", &DTMIntersection::index_to_ter)
        .def("get_footprint_corners", &DTMIntersection::get_footprint_corners)
        .def("interpolate", &DTMIntersection::interpolate)
        .def("intersect_dtm_cube", &DTMIntersection::intersect_dtm_cube)
        .def("intersection", &DTMIntersection::intersection)
        .def("intersection_n_los_dtm", &DTMIntersection::intersection_n_los_dtm)
        .def("get_alt_data", &DTMIntersection::get_alt_data)
        .def("get_alt_min", &DTMIntersection::get_alt_min)
        .def("get_alt_max", &DTMIntersection::get_alt_max)
        .def("get_plane_coef_a", &DTMIntersection::get_plane_coef_a)
        .def("get_plane_coef_b", &DTMIntersection::get_plane_coef_b)
        .def("get_plane_coef_c", &DTMIntersection::get_plane_coef_c)
        .def("get_plane_coef_d", &DTMIntersection::get_plane_coef_d)
        .def("get_alt_min_cell", &DTMIntersection::get_alt_min_cell)
        .def("get_alt_max_cell", &DTMIntersection::get_alt_max_cell)
        .def("get_tol_z", &DTMIntersection::get_tol_z)
        .def("get_epsg", &DTMIntersection::get_epsg)
        .def("get_plans", &DTMIntersection::get_plans)
        .def("get_trans_inv", &DTMIntersection::get_trans_inv)
        .def("get_transform", &DTMIntersection::get_transform)
        .def("get_nb_rows", &DTMIntersection::get_nb_rows)
        .def("get_nb_columns", &DTMIntersection::get_nb_columns)
        .def("set_alt_data", &DTMIntersection::set_alt_data)
        .def("set_alt_min", &DTMIntersection::set_alt_min)
        .def("set_alt_max", &DTMIntersection::set_alt_max)
        .def("set_plane_coef_a", &DTMIntersection::set_plane_coef_a)
        .def("set_plane_coef_b", &DTMIntersection::set_plane_coef_b)
        .def("set_plane_coef_c", &DTMIntersection::set_plane_coef_c)
        .def("set_plane_coef_d", &DTMIntersection::set_plane_coef_d)
        .def("set_alt_min_cell", &DTMIntersection::set_alt_min_cell)
        .def("set_alt_max_cell", &DTMIntersection::set_alt_max_cell)
        .def("set_tol_z", &DTMIntersection::set_tol_z)
        .def("set_epsg", &DTMIntersection::set_epsg)
        .def("set_plans", &DTMIntersection::set_plans)
        .def("set_trans_inv", &DTMIntersection::set_trans_inv)
        .def("set_transform", &DTMIntersection::set_transform)
        .def("set_nb_rows", &DTMIntersection::set_nb_rows)
        .def("set_nb_columns", &DTMIntersection::set_nb_columns)
        .def(py::pickle(
                [](const DTMIntersection &p) { // __getstate__
                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(
                        p.get_alt_data(),
                        p.get_alt_min(),
                        p.get_alt_max(),
                        p.get_plane_coef_a(),
                        p.get_plane_coef_b(),
                        p.get_plane_coef_c(),
                        p.get_plane_coef_d(),
                        p.get_alt_min_cell(),
                        p.get_alt_max_cell(),
                        p.get_tol_z(),
                        p.get_epsg(),
                        p.get_plans(),
                        p.get_trans_inv(),
                        p.get_transform(),
                        p.get_nb_rows(),
                        p.get_nb_columns());
                },
                [](py::tuple t) { // __setstate__
                if (t.size() != 16)
                        throw std::runtime_error("Invalid state!");

                /* Create a new C++ instance */
                DTMIntersection p;

                /* Assign any additional state */
                p.set_alt_data(t[0].cast<std::vector<double>>());
                p.set_alt_min(t[1].cast<double>());
                p.set_alt_max(t[2].cast<double>());
                p.set_plane_coef_a(t[3].cast<std::array<double,6>>());
                p.set_plane_coef_b(t[4].cast<std::array<double,6>>());
                p.set_plane_coef_c(t[5].cast<std::array<double,6>>());
                p.set_plane_coef_d(t[6].cast<std::array<double,6>>());
                p.set_alt_min_cell(t[7].cast<std::vector<double>>());
                p.set_alt_max_cell(t[8].cast<std::vector<double>>());
                p.set_tol_z(t[9].cast<double>());
                p.set_epsg(t[10].cast<int>());
                p.set_plans(t[11].cast<std::vector<double>>());
                p.set_trans_inv(t[12].cast<std::array<double,6>>());
                p.set_transform(t[13].cast<std::array<double,6>>());
                p.set_nb_rows(t[14].cast<int>());
                p.set_nb_columns(t[15].cast<int>());
                return p;
                }
        ));

    py::class_<GeoModelTemplate>(m, "GeoModelTemplate")
        .def(py::init<>())
        .def("direct_loc_h", py::overload_cast<double,
                                        double,
                                        double,
                                        bool,
                                        bool>(&GeoModelTemplate::direct_loc_h, py::const_))

        .def("direct_loc_h", py::overload_cast<std::vector<double> const&,
                                        std::vector<double> const&,
                                        std::vector<double> const&,
                                        bool,
                                        bool>(&GeoModelTemplate::direct_loc_h, py::const_))

        .def("direct_loc_dtm", py::overload_cast<double,
                                        double,
                                        DTMIntersection const&>
                                        (&GeoModelTemplate::direct_loc_dtm, py::const_))

        .def("direct_loc_dtm", py::overload_cast<std::vector<double> const&,
                                        std::vector<double> const&,
                                        DTMIntersection const&>
                                        (&GeoModelTemplate::direct_loc_dtm, py::const_))

        .def("inverse_loc",py::overload_cast<double,
                                        double,
                                        double>
                                        (&GeoModelTemplate::inverse_loc, py::const_))

        .def("inverse_loc",py::overload_cast<std::vector<double> const&,
                                        std::vector<double> const&,
                                        std::vector<double> const&>
                                        (&GeoModelTemplate::inverse_loc, py::const_));


    py::class_<RPC,GeoModelTemplate>(m, "RPC")
        .def(py::init<bool,
                bool,
                std::array<double, 20>,
                std::array<double, 20>,
                std::array<double, 20>,
                std::array<double, 20>,
                std::array<double, 20>,
                std::array<double, 20>,
                std::array<double, 20>,
                std::array<double, 20>,
                std::array<double, 10>>())
        .def("direct_loc_h", py::overload_cast<double,
                                                double,
                                                double,
                                                bool,
                                                bool>
                                                (&RPC::direct_loc_h, py::const_))

        .def("direct_loc_h", py::overload_cast<std::vector<double> const&,
                                                std::vector<double> const&,
                                                std::vector<double> const&,
                                                bool,
                                                bool>
                                                (&RPC::direct_loc_h, py::const_))

        .def("direct_loc_dtm", py::overload_cast<double,
                                                double,
                                                DTMIntersection const&>
                                                (&RPC::direct_loc_dtm, py::const_))

        .def("direct_loc_dtm", py::overload_cast<std::vector<double> const&,
                                                std::vector<double> const&,
                                                DTMIntersection const&>
                                                (&RPC::direct_loc_dtm, py::const_))

        .def("inverse_loc",py::overload_cast<double,
                                        double,
                                        double>
                                        (&RPC::inverse_loc, py::const_))

        .def("inverse_loc",py::overload_cast<std::vector<double> const&,
                                        std::vector<double> const&,
                                        std::vector<double> const&>
                                        (&RPC::inverse_loc, py::const_))

        .def("compute_loc_inverse_derivates", &RPC::compute_loc_inverse_derivates)

        .def("direct_loc_inverse_iterative",py::overload_cast<double,double,double,int,bool>\
                (&RPC::direct_loc_inverse_iterative, py::const_))

        .def("direct_loc_inverse_iterative",py::overload_cast<std::vector<double> const&,\
                std::vector<double> const&,std::vector<double> const&,int,bool>\
                (&RPC::direct_loc_inverse_iterative, py::const_))

        .def("get_alt_min_max", &RPC::get_alt_min_max)
        .def("los_extrema", &RPC::los_extrema)
        .def("compute_rational_function_polynomial_unitary",
                &RPC::compute_rational_function_polynomial_unitary)

        .def("compute_rational_function_polynomial", 
                &RPC::compute_rational_function_polynomial)
        .def("get_num_col", &RPC::get_num_col)
        .def("get_den_col", &RPC::get_den_col)
        .def("get_num_row", &RPC::get_num_row)
        .def("get_den_row", &RPC::get_den_row)
        .def("get_num_lon", &RPC::get_num_lon)
        .def("get_den_lon", &RPC::get_den_lon)
        .def("get_num_lat", &RPC::get_num_lat)
        .def("get_den_lat", &RPC::get_den_lat)
        .def("get_alt_minmax", &RPC::get_alt_minmax)
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

    m.def("polynomial_equation", &polynomial_equation, "Compute polynomial equation");

    m.def("derivative_polynomial_latitude", &derivative_polynomial_latitude,
            "Compute latitude derivative polynomial equation");

    m.def("derivative_polynomial_longitude", &derivative_polynomial_longitude,
    "Compute longitude derivative polynomial equation");

    m.def("init_min_max", &init_min_max,
    "init_min_max");
    
    m.def("compute_epipolar_angle", &compute_epipolar_angle,
            "compute epipolar angle");




m.def("coloc", py::overload_cast <GeoModelTemplate const&,
                                GeoModelTemplate const&,
                                vector<double> const&,
                                vector<double> const&,
                                DTMIntersection const&>(&coloc));

m.def("coloc", py::overload_cast<GeoModelTemplate const&,
                                GeoModelTemplate const&,
                                double,
                                double,
                                DTMIntersection const&>(&coloc));

m.def("coloc", py::overload_cast<GeoModelTemplate const&,
                                GeoModelTemplate const&,
                                vector<double> const&,
                                vector<double> const&,
                                vector<double> const&>(&coloc));

m.def("coloc", py::overload_cast<GeoModelTemplate const&,
                                GeoModelTemplate const&,
                                double,
                                double,
                                double>(&coloc));

m.def("moving_along_axis", &moving_along_axis<DTMIntersection const&>);
m.def("moving_along_axis", &moving_along_axis<double>);

m.def("compute_local_epipolar_line", &compute_local_epipolar_line<DTMIntersection const&>);
m.def("compute_local_epipolar_line", &compute_local_epipolar_line<double>);

m.def("compute_strip_of_epipolar_grid", &compute_strip_of_epipolar_grid<DTMIntersection const&>,
        "compute_strip_of_epipolar_grid");
m.def("compute_strip_of_epipolar_grid", &compute_strip_of_epipolar_grid<double>,
        "compute_strip_of_epipolar_grid");
}

// Commande to compile c++ bindings
// c++ -w -O3 -Wall -Wextra -shared -std=c++20 -march=native -fPIC
//$(python3 -m pybind11 --includes) bind.cpp -o bindings_cpp$(python3-config --extension-suffix)
