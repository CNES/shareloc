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
    map<string, double> rpc_params;//map<string, double> is a simple dict -> maybe inappropiate here
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
    double offset_lon;
    double scale_lon;
    double offset_lat;
    double scale_lat;


public:
    using GeoModelTemplate::GeoModelTemplate;

    /**Constructor*/
    RPC(array<double, 20> num_col,
        array<double, 20> den_col,
        array<double, 20> num_row,
        array<double, 20> den_row,
        array<double, 10> norm_coeffs);

    /**direct_loc_h*/
    vector<vector<double>> direct_loc_h(
        vector<double> row,
        vector<double> col,
        double alt,
        bool fill_nan=false);

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
        string dtm);// dtm is a python class not a string

    /**inverse_loc*/
    tuple<vector<double>,vector<double>,vector<double>> inverse_loc(
        vector<double> lon,
        vector<double> lat,
        double alt);

    /**filter_coordinates*/
    vector<vector<double>> filter_coordinates(
        vector<double> first_coord,
        vector<double> second_coord,
        bool fill_nan=false,
        string direction="direct");

    /**compute_loc_inverse_derivates*/
    tuple<vector<double>,
    vector<double>,
    vector<double>,
    vector<double>> compute_loc_inverse_derivates(
        vector<double> lon,
        vector<double> lat,
        vector<double> alt);

    /**direct_loc_inverse_iterative*/
    vector<vector<double>> direct_loc_inverse_iterative(
        vector<double> row,
        vector<double> col,
        double alt,
        int nb_iter_max=10,
        bool fill_nan=false);

    /**get_alt_min_max*/
    vector<double> get_alt_min_max();

    /**los_extrema*/
    vector<vector<double>> los_extrema(
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


