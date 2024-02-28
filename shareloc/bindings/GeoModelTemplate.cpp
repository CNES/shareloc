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

#include "GeoModelTemplate.hpp"

GeoModelTemplate::GeoModelTemplate() {

}
GeoModelTemplate::~GeoModelTemplate() {

}

std::tuple<double,double,double>\
 GeoModelTemplate::direct_loc_h(
    double row,
    double col,
    double alt, 
    bool fill_nan,
    bool using_direct_coef)const{
    std::tuple<double,double,double> vect;
    return vect;
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>\
 GeoModelTemplate::direct_loc_h(
    std::vector<double> const& row,
    std::vector<double> const& col,
    std::vector<double> const& alt, 
    bool fill_nan,
    bool using_direct_coef)const{
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> vect;
    return vect;
}

std::tuple<double,double,double>\
 GeoModelTemplate::direct_loc_dtm(
    double row,
    double col,
    DTMIntersection const& dtm) const{
    std::tuple<double,double,double> vect;
    return vect;
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>\
 GeoModelTemplate::direct_loc_dtm(
    std::vector<double> const& row,
    std::vector<double> const& col,
    DTMIntersection const& dtm) const{
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> vect;
    return vect;
}

std::tuple<double,double,double>\
 GeoModelTemplate::inverse_loc(
    double lon,
    double lat,
    double alt)const{
    std::tuple<double,double,double> res;
    return res;
    }

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>\
 GeoModelTemplate::inverse_loc(
    std::vector<double> const& lon,
    std::vector<double> const& lat,
    std::vector<double> const& alt)const{
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> res;
    return res;
}