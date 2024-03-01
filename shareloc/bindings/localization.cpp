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

#include <vector>
#include <array>
#include <cmath>
#include <iostream>

#include "GeoModelTemplate.hpp"

using namespace std;

tuple<vector<double>,vector<double>,vector<double>> coloc(GeoModelTemplate const& geom1,
                                        GeoModelTemplate const& geom2,
                                        vector<double> const& vector_row,
                                        vector<double> const& vector_col,
                                        DTMIntersection const& elevation){

    size_t nb_points = vector_row.size();
    vector<double> lon(nb_points);
    vector<double> lat(nb_points);
    vector<double> alt(nb_points);

    tie(lon, lat, alt) = geom1.direct_loc_dtm(vector_row,vector_col, elevation);
    return geom2.inverse_loc(lon, lat, alt);
}

tuple<double,double,double> coloc(GeoModelTemplate const& geom1,
                                    GeoModelTemplate const& geom2,
                                    double row,
                                    double col,
                                    DTMIntersection const& elevation){

    double lon;
    double lat;
    double alt;

    tie(lon, lat, alt) = geom1.direct_loc_dtm(row,col, elevation);
    return geom2.inverse_loc(lon, lat, alt);
}

tuple<vector<double>,vector<double>,vector<double>> coloc(GeoModelTemplate const& geom1,
                                        GeoModelTemplate const& geom2,
                                        vector<double> const& vector_row,
                                        vector<double> const& vector_col,
                                        vector<double> const& elevation){

    size_t nb_points = vector_row.size();
    vector<double> lon(nb_points);
    vector<double> lat(nb_points);
    vector<double> alt(nb_points);

    tie(lon, lat, alt) = geom1.direct_loc_h(vector_row,vector_col,elevation);
    return geom2.inverse_loc(lon, lat, alt);
}

tuple<double,double,double> coloc(GeoModelTemplate const& geom1,
                                    GeoModelTemplate const& geom2,
                                    double row,
                                    double col,
                                    double elevation){

    double lon;
    double lat;
    double alt;

    tie(lon, lat, alt) = geom1.direct_loc_h(row,col,elevation);
    return geom2.inverse_loc(lon, lat, alt);
}