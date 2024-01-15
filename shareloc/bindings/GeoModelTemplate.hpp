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
Cpp copy of GeoModelTemplate.py
Abstract class GeoModelTemplate.
    Child class: RPC
*/

#include <string>
#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace std;

/**
Abstract class GeoModelTemplate:
    Child: RPC
*/
class GeoModelTemplate {

private:
    string type;
    int epsg;
public:
    /**Constructor*/
    GeoModelTemplate();
    /**Destructor*/
    ~GeoModelTemplate();
    /**direct_loc_h*/
    virtual vector<vector<double>> direct_loc_h(
        vector<double> row,
        vector<double> col,
        double alt,
        bool fill_nan=false);
    /**direct_loc_dtm*/
    virtual vector<vector<double>> direct_loc_dtm(
        vector<double> row,
        vector<double> col,
        string dtm);
    /**inverse_loc*/
    virtual tuple<vector<double>,vector<double>,vector<double>> inverse_loc(
        vector<double> lon,
        vector<double> lat,
        vector<double> alt);
};
