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

#include <string>
#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <array>
#include <algorithm>
#include <cmath>

using namespace std;

/**
Class DTMIntersection
Framework of the DTMIntersection python class.
*/

class DTMIntersection {

private:

        string dtm_file;
        vector<double> alt_data;
        double alt_min;
        double alt_max;
        double origin_x;
        double origin_y;
        double pixel_size_x;
        double pixel_size_y;
        vector<double> plane_coef_a;
        vector<double> plane_coef_b;
        vector<double> plane_coef_c;
        vector<double> plane_coef_d;
        double alt_min_cell;
        double alt_max_cell;
        double tol_z;// = 0.0001

        int epsg;

        vector<double> grid_row;
        vector<double> grid_col;

        vector<double> plans;

        vector<double> trans_inv; //affine.affine en python
        vector<double> transform;
        int nb_rows;
        int nb_columns;



public:


    /**Constructor*/
    DTMIntersection(array<double, 20> dtm_image);//determiner comment passer les arg

    /**eq_plan*/
    double eq_plan(int i, array<double, 3> position);

    /**ter_to_index*/
    array<double, 3> ter_to_index(array<double, 3> vect_ter);

    /**ter_to_indexs*/
    vector<double> ter_to_indexs(vector<double> vect_ter);

    /**index_to_ter*/
    array<double, 3> index_to_ter(array<double, 3> vect_ter);

    /**get_alt_offset*/
    array<double, 2> get_alt_offset(int epsg);//maybe unecessary

    /**interpolate*/
    double interpolate(double pos_row, double pos_col);

    /**intersect_dtm_cube*/
    tuple<bool,bool,vector<double>,bool,vector<double>> intersect_dtm_cube(vector<double> los);

    /**intersection*/
    tuple<bool,bool,vector<double>> intersection(
        vector<double> los_index,
        vector<double> point_b, 
        double h_intersect);

    //-- getter --//


    /**get_dtm_file*/
    string get_dtm_file();
    /**get_alt_data*/
    vector<double> get_alt_data();
    /**get_alt_min*/
    double get_alt_min();
    /**get_alt_max*/
    double get_alt_max();
    /**get_origin_x*/
    double get_origin_x();
    /**get_origin_y*/
    double get_origin_y();
    /**get_pixel_size_x*/
    double get_pixel_size_x();
    /**get_pixel_size_y*/
    double get_pixel_size_y();
    /**get_plane_coef_a*/
    vector<double> get_plane_coef_a();
    /**get_plane_coef_b*/
    vector<double> get_plane_coef_b();
    /**get_plane_coef_c*/
    vector<double> get_plane_coef_c();
    /**get_plane_coef_d*/
    vector<double> get_plane_coef_d();
    /**get_alt_min_cell*/
    double get_alt_min_cell();
    /**get_alt_max_cell*/
    double get_alt_max_cell();
    /**get_tol_z*/
    double get_tol_z();// = 0.0001
    /**get_epsg*/
    int get_epsg();
    /**get_grid_row*/
    vector<double> get_grid_row();
    /**get_grid_col*/
    vector<double> get_grid_col();
    /**get_plans*/
    vector<double> get_plans();
    /**get_trans_inv*/
    vector<double> get_trans_inv(); //affine.affine en python
    /**get_transform*/
    vector<double> get_transform();
    /**get_nb_rows*/
    int get_nb_rows();
    /**get_nb_columns*/
    int get_nb_columns();
};