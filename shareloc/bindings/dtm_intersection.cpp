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
Cpp copy of dtm_intersection.py
*/
#include "dtm_intersection.hpp"

//---- DTMIntersection methodes ----//


DTMIntersection::DTMIntersection(array<double, 20> dtm_image){//determiner comment passer les arg
    cout<<"Constructor DTMIntersection"<<endl;
}


double DTMIntersection::eq_plan(int i, array<double, 3> position){
    double res;
    return res;
}

array<double, 3> DTMIntersection::ter_to_index(array<double, 3> vect_ter){
    array<double, 3> res;
    return res;
}

vector<double> DTMIntersection::ter_to_indexs(vector<double> vect_ter){
    vector<double> res;
    return res;
}

array<double, 3> DTMIntersection::index_to_ter(array<double, 3> vect_ter){
    array<double, 3> res;
    return res;
}

array<double, 2> DTMIntersection::get_alt_offset(int epsg){//maybe unecessary
    array<double, 2> res;
    return res;
}

double DTMIntersection::interpolate(double pos_row, double pos_col){
    double res;
    return res;
}

tuple<bool, 
bool,
vector<double>,
bool,
vector<double>> DTMIntersection::intersect_dtm_cube(vector<double> los){
    tuple<bool,bool,vector<double>,bool,vector<double>> res;
    return res;
}

tuple<bool,
bool,
vector<double>> DTMIntersection::intersection(
    vector<double> los_index,
    vector<double> point_b, 
    double h_intersect){
    tuple<bool,bool,vector<double>> res;
    return res;
}




int main(){return 0;}// to delete