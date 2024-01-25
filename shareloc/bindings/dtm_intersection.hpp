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
        array<double,6> plane_coef_a;
        array<double,6> plane_coef_b;
        array<double,6> plane_coef_c;
        array<double,6> plane_coef_d;
        vector<double> alt_min_cell;
        vector<double> alt_max_cell;
        double tol_z;// = 0.0001

        int epsg;

        vector<double> plans;

        array<double,6> trans_inv; //affine.affine en python
        array<double,6> transform;
        int nb_rows;
        int nb_columns;



public:


    /**Constructor*/
    DTMIntersection(
        int dtm_image_epsg,
        vector<double> dtm_image_alt_data,
        int dtm_image_nb_rows,
        int dtm_image_nb_columns,
        tuple<double,double,double,double,double,double> dtm_image_transform
    );//determiner comment passer les arg

    /**eq_plan*/
    double eq_plan(int i, array<double, 3> const& position);

    /**ter_to_index*/
    array<double, 3> ter_to_index(array<double, 3> const& vect_ter);

    /**ter_to_indexs*/
    vector<double> ter_to_indexs(vector<double> const& vect_ter);

    /**index_to_ter*/
    array<double, 3> index_to_ter(array<double, 3> const& vect_ter);

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
    /**get_plane_coef_a*/
    array<double,6> get_plane_coef_a();
    /**get_plane_coef_b*/
    array<double,6> get_plane_coef_b();
    /**get_plane_coef_c*/
    array<double,6> get_plane_coef_c();
    /**get_plane_coef_d*/
    array<double,6> get_plane_coef_d();
    /**get_alt_min_cell*/
    vector<double> get_alt_min_cell();
    /**get_alt_max_cell*/
    vector<double> get_alt_max_cell();
    /**get_tol_z*/
    double get_tol_z();// = 0.0001
    /**get_epsg*/
    int get_epsg();
    /**get_plans*/
    vector<double> get_plans();
    /**get_trans_inv*/
    array<double,6> get_trans_inv(); //affine.affine en python
    /**get_transform*/
    array<double,6> get_transform();
    /**get_nb_rows*/
    int get_nb_rows();
    /**get_nb_columns*/
    int get_nb_columns();
};



//-- Function --//

/**init_min_max*/
tuple<vector<double>,vector<double>> init_min_max(vector<double> const& alt_data,
                                                    int nb_rows,
                                                    int nb_columns);