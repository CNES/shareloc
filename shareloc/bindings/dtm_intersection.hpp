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

#ifndef DTM_INTERSECTION_H  
#define DTM_INTERSECTION_H


/**
Class DTMIntersection
Framework of the DTMIntersection python class.
*/

class DTMIntersection {

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
    string const& get_dtm_file() const noexcept {return m_dtm_file;};
    /**get_alt_data*/
    vector<double>  const&  get_alt_data() const noexcept {return m_alt_data;};
    /**get_alt_min*/
    double get_alt_min() const noexcept {return m_alt_min;};
    /**get_alt_max*/
    double get_alt_max() const noexcept {return m_alt_max;};
    /**get_origin_x*/
    double get_origin_x() const noexcept {return m_origin_x;};
    /**get_origin_y*/
    double get_origin_y() const noexcept {return m_origin_y;};
    /**get_pixel_size_x*/
    double get_pixel_size_x() const noexcept {return m_pixel_size_x;};
    /**get_pixel_size_y*/
    double get_pixel_size_y() const noexcept {return m_pixel_size_y;};
    /**get_plane_coef_a*/
    vector<double> const& get_plane_coef_a() const noexcept {return m_plane_coef_a;};
    /**get_plane_coef_b*/
    vector<double> const& get_plane_coef_b() const noexcept {return m_plane_coef_b;};
    /**get_plane_coef_c*/
    vector<double> const& get_plane_coef_c() const noexcept {return m_plane_coef_c;};
    /**get_plane_coef_d*/
    vector<double> const& get_plane_coef_d() const noexcept {return m_plane_coef_d;};
    /**get_alt_min_cell*/
    double get_alt_min_cell() const noexcept {return m_alt_min_cell;};
    /**get_alt_max_cell*/
    double get_alt_max_cell() const noexcept {return m_alt_max_cell;};
    /**get_tol_z*/
    double get_tol_z() const noexcept {return m_tol_z;};// = 0.0001
    /**get_epsg*/
    int get_epsg() const noexcept {return m_epsg;};
    /**get_grid_row*/
    vector<double> const& get_grid_row() const noexcept {return m_grid_row;};
    /**get_grid_col*/
    vector<double> const& get_grid_col() const noexcept {return m_grid_col;};
    /**get_plans*/
    vector<double> const& get_plans() const noexcept {return m_plans;};
    /**get_trans_inv*/
    vector<double> const& get_trans_inv() const noexcept {return m_trans_inv;}; //affine.affine en python
    /**get_transform*/
    vector<double> const& get_transform() const noexcept {return m_transform;};
    /**get_nb_rows*/
    int get_nb_rows() const noexcept {return m_nb_rows;};
    /**get_nb_columns*/
    int get_nb_columns() const noexcept {return m_nb_columns;};

private:

        string m_dtm_file;
        vector<double> m_alt_data;
        double m_alt_min;
        double m_alt_max;
        double m_origin_x;
        double m_origin_y;
        double m_pixel_size_x;
        double m_pixel_size_y;
        vector<double> m_plane_coef_a;
        vector<double> m_plane_coef_b;
        vector<double> m_plane_coef_c;
        vector<double> m_plane_coef_d;
        double m_alt_min_cell;
        double m_alt_max_cell;
        double m_tol_z;// = 0.0001

        int m_epsg;

        vector<double> m_grid_row;
        vector<double> m_grid_col;

        vector<double> m_plans;

        vector<double> m_trans_inv; //affine.affine en python
        vector<double> m_transform;
        int m_nb_rows;
        int m_nb_columns;
};


#endif