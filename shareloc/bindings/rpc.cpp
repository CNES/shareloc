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
Cpp copy of rpc.py
*/

#include "rpc.hpp"

//---- RPC methodes ----//

RPC::RPC(array<double, 20> num_col,
        array<double, 20> den_col,
        array<double, 20> num_row,
        array<double, 20> den_row,
        array<double, 10> norm_coeffs):GeoModelTemplate(){

    num_col = num_col;
    den_col = den_col;
    num_row = num_row;
    den_row = den_row;

    offset_lon = norm_coeffs[0];
    scale_lon = norm_coeffs[1];
    offset_lat = norm_coeffs[2];
    scale_lat = norm_coeffs[3];
    offset_alt = norm_coeffs[4];
    scale_alt = norm_coeffs[5];
    offset_col = norm_coeffs[6];
    scale_col = norm_coeffs[7];
    offset_row = norm_coeffs[8];
    scale_row = norm_coeffs[9];

}

vector<vector<double>> RPC::direct_loc_h(
    vector<double> row,
    vector<double> col,
    double alt,
    bool fill_nan){
    vector<vector<double>> vect;
    return vect;
}

tuple<vector<vector<double>>,vector<vector<double>>> RPC::direct_loc_grid_h(
    int row0,
    int col0,
    int steprow,
    int stepcol,
    int nbrow,
    int nbcol,
    double alt){
    tuple<vector<vector<double>>,vector<vector<double>>> res;
    return res;
}

vector<vector<double>> RPC::direct_loc_dtm(
    double row,
    double col,
    string dtm){// dtm intersection model ->python class
    vector<vector<double>> vect;
    return vect;
}

tuple<vector<double>,vector<double>,vector<double>> RPC::inverse_loc(
    vector<double> lon,
    vector<double> lat,
    double alt){
    tuple<vector<double>,vector<double>,vector<double>> vect;
    return vect;
    }


vector<vector<double>> RPC::filter_coordinates(
    vector<double> first_coord,
    vector<double> second_coord,
    bool fill_nan,
    string direction){
    vector<vector<double>> vect;
    return vect;
}

tuple<vector<double>,
vector<double>,
vector<double>,
vector<double>> RPC::compute_loc_inverse_derivates(
    vector<double> lon,
    vector<double> lat,
    vector<double> alt){
    tuple<vector<double>,vector<double>,vector<double>,vector<double>> res;
    return res;
}

vector<vector<double>> RPC::direct_loc_inverse_iterative(
    vector<double> row,
    vector<double> col,
    double alt,
    int nb_iter_max,
    bool fill_nan){
    vector<vector<double>> vect;
    return vect;
}

vector<double> RPC::get_alt_min_max(){
    vector<double> vect;
    return vect;
}

vector<vector<double>> RPC::los_extrema(
    double row,
    double col,
    double alt_min,
    double alt_max,
    bool fill_nan,
    int epsg){
    vector<vector<double>> vect;
    return vect;
}

//---- Functions ----//

double polynomial_equation(
    double xnorm,
    double ynorm,
    double znorm,
    vector<double> coeff){
    double res;
    return res;
}


tuple<vector<double>,vector<double>> compute_rational_function_polynomial(
    vector<double> lon_col_norm,
    vector<double> lat_row_norm,
    vector<double> alt_norm,
    vector<double> num_col,
    vector<double> den_col,
    vector<double> num_lin,
    vector<double> den_lin,
    double scale_col,
    double offset_col,
    double scale_lin,
    double offset_lin
){
    tuple<vector<double>,vector<double>> res;
    return res;
}


double derivative_polynomial_latitude(
    double lon_norm,
    double lat_norm,
    double alt_norm,
    vector<double> coeff){
    double res;
    return res;
}


double derivative_polynomial_longitude(
    double lon_norm,
    double lat_norm,
    double alt_norm,
    vector<double> coeff){
    double res;
    return res;
}


tuple<vector<double>,
vector<double>,
vector<double>,
vector<double>> compute_loc_inverse_derivates_numba(
    vector<double> lon_norm,
    vector<double> lat_norm,
    vector<double> alt_norm,
    vector<double> num_col,
    vector<double> den_col,
    vector<double> num_lin,
    vector<double> den_lin,
    double scale_col,
    double scale_lon,
    double scale_lin,
    double scale_lat
){
    tuple<vector<double>,vector<double>,vector<double>,vector<double>> res;
    return res;
}