/*
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


#include "GeoModelTemplate.hpp"
#include "GeoModelTemplate.cpp"

/**
Class RPC
Framework of the RPC python class.
*/

class RPC : public GeoModelTemplate{

private:

    string type;
    int epsg;
    string datum;
    map<string, double> rpc_params;//map<string, double> is a simple dict -> maybe inappropiate
    double lim_extrapol;

    vector<vector<double>> monomes;// to convert to array
    vector<vector<double>> monomes_deriv_1;
    vector<vector<double>> monomes_deriv_2;

    bool inverse_coefficient;
    bool direct_coefficient;

    array<double, 20> num_col;
    array<double, 20> den_col;
    array<double, 20> num_row;
    array<double, 20> den_row;

    array<double, 20> num_lon;
    array<double, 20> den_lon;
    array<double, 20> num_lat;
    array<double, 20> den_lat;

    array<double, 2> alt_minmax;

    double col0;
    double colmax;
    double row0;
    double rowmax;

    double offset_row;
    double scale_row;
    double offset_col;
    double scale_col;
    double offset_alt;
    double scale_alt;
    double offset_lon;//py: lon=x
    double scale_lon;//py: lon=x
    double offset_lat;//py: lat=y
    double scale_lat;//py: lat=y


public:
    using GeoModelTemplate::GeoModelTemplate;

    /**Constructor*/
    RPC(bool inverse_coefficient,
        bool direct_coefficient,
        array<double, 20> num_col,
        array<double, 20> den_col,
        array<double, 20> num_row,
        array<double, 20> den_row,
        array<double, 20> num_lon,
        array<double, 20> den_lon,
        array<double, 20> num_lat,
        array<double, 20> den_lat,
        array<double, 10> norm_coeffs);

    /**direct_loc_h*/
    tuple<vector<double>,vector<double>,vector<double>> direct_loc_h(
        vector<double> row,
        vector<double> col,
        vector<double> alt,
        bool fill_nan=false);//override

    /**direct_loc_grid_h*/
    tuple<vector<vector<double>>,vector<vector<double>>> direct_loc_grid_h(
        int row0,
        int col0,
        int steprow,
        int stepcol,
        int nbrow,
        int nbcol,
        double alt);

    /**direct_loc_dtm*/
    vector<vector<double>> direct_loc_dtm(
        double row,
        double col,
        string dtm);//override + dtm is a python class not a string

    /**inverse_loc unitary*/
    tuple<double,double,double> inverse_loc(
        double lon,
        double lat,
        double alt);

    /**inverse_loc*/
    tuple<vector<double>,vector<double>,vector<double>> inverse_loc(
        vector<double> lon,
        vector<double> lat,
        vector<double> alt) override;

    /**filter_coordinates*/
    tuple<vector<bool>,vector<double>,vector<double>> filter_coordinates(
        vector<double> first_coord,
        vector<double> second_coord,
        bool fill_nan=false,
        string direction="direct");


    /**compute_loc_inverse_derivates unitary*/
    tuple<double,
    double,
    double,
    double> compute_loc_inverse_derivates(
        double lon,
        double lat,
        double alt);

    /**direct_loc_inverse_iterative*/
    tuple<vector<double>,vector<double>,vector<double>> direct_loc_inverse_iterative(
        vector<double> row,
        vector<double> col,
        vector<double> alt,
        int nb_iter_max=10,
        bool fill_nan=false);

    /**get_alt_min_max*/
    array<double, 2> get_alt_min_max();

    /**los_extrema*/
    tuple<vector<double>,vector<double>,vector<double>> los_extrema(
        double row,
        double col,
        double alt_min,
        double alt_max,
        bool fill_nan=false,
        int epsg=4326);
    

    //-- getter --//

    /**get_num_col*/
    array<double, 20> get_num_col();
    /**get_den_col*/
    array<double, 20> get_den_col();
    /**get_num_row*/
    array<double, 20> get_num_row();
    /**get_den_row*/
    array<double, 20> get_den_row();

    /**get_num_lon*/
    array<double, 20> get_num_lon();
    /**get_den_lon*/
    array<double, 20> get_den_lon();
    /**get_num_lat*/
    array<double, 20> get_num_lat();
    /**get_den_lat*/
    array<double, 20> get_den_lat();

    /**get_offset_row*/
    double get_offset_row();
    /**get_scale_row*/
    double get_scale_row();
    /**get_offset_col*/
    double get_offset_col();
    /**get_scale_col*/
    double get_scale_col();
    /**get_offset_alt*/
    double get_offset_alt();
    /**get_scale_alt*/
    double get_scale_alt();
    /**get_offset_lon*/
    double get_offset_lon();
    /**get_scale_lon*/
    double get_scale_lon();
    /**get_offset_lat*/
    double get_offset_lat();
    /**get_scale_lat*/
    double get_scale_lat();
};

// function 

/**Compute polynomial equation"*/
double polynomial_equation(
    double xnorm,
    double ynorm,
    double znorm,
    const array<double, 20>* coeff);

/** compute_rational_function_polynomial unitary*/
tuple<double,double,double> compute_rational_function_polynomial_unitary(
    double lon_col_norm,
    double lat_row_norm,
    double alt_norm,
    array<double, 20> num_col,
    array<double, 20> den_col,
    array<double, 20> num_lin,
    array<double, 20> den_lin,

    //input
    double scale_lon_col,
    double offset_lon_col,
    double scale_lat_row,
    double offset_lat_row,
    double scale_alt,
    double offset_alt,

    //output
    double scale_col,
    double offset_col,
    double scale_lin,
    double offset_lin
);


/**Compute rational function polynomial. Useful to compute direct and inverse localization
        "using direct or inverse RPC."*/
tuple<vector<double>,vector<double>,vector<double>> compute_rational_function_polynomial(
    vector<double> lon_col_norm,
    vector<double> lat_row_norm,
    vector<double> alt_norm,
    array<double, 20> num_col,
    array<double, 20> den_col,
    array<double, 20> num_lin,
    array<double, 20> den_lin,

    //input
    double scale_lon_col,
    double offset_lon_col,
    double scale_lat_row,
    double offset_lat_row,
    double scale_alt,
    double offset_alt,

    //output
    double scale_col,
    double offset_col,
    double scale_lin,
    double offset_lin
);

/**Check if arrays have the same size and cut it if needed*/
tuple<vector<double>,
    vector<double>,
    vector<double>> check_sizes(vector<double> lon_col,
    vector<double> lat_row,
    vector<double>alt);

/** Compute derivative_polynomial_latitude*/
double derivative_polynomial_latitude(
    double lon_norm,
    double lat_norm,
    double alt_norm,
    const array<double, 20>* coeff);

/**Compute derivative_polynomial_longitude*/
double derivative_polynomial_longitude(
    double lon_norm,
    double lat_norm,
    double alt_norm,
    const array<double, 20>* coeff);