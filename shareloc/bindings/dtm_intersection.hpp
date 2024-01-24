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

#include <vector>
#include <string>
#include <array>
#include <tuple>

/**
  Class DTMIntersection
  Framework of the DTMIntersection python class.
 */

class DTMIntersection
{

public:

    /**Constructor*/
    DTMIntersection(std::array<double, 20> const& dtm_image);//determiner comment passer les arg

    /**eq_plan*/
    double eq_plan(int i, std::array<double, 3> position) const;

    /**ter_to_index*/
    std::array<double, 3> ter_to_index(std::array<double, 3> vect_ter) const;

    /**ter_to_indexs*/
    std::vector<double> ter_to_indexs(std::vector<double> const& vect_ter) const;

    /**index_to_ter*/
    std::array<double, 3> index_to_ter(std::array<double, 3> vect_ter) const;

    /**get_alt_offset*/
    std::array<double, 2> get_alt_offset(int epsg) const;//maybe unecessary

    /**interpolate*/
    double interpolate(double pos_row, double pos_col) const;

    /**intersect_dtm_cube*/
    std::tuple<bool,bool,std::vector<double>,bool,std::vector<double>> intersect_dtm_cube(std::vector<double> const& los) const;

    /**intersection*/
    std::tuple<bool,bool,std::vector<double>> intersection(
        std::vector<double> const& los_index,
        std::vector<double> const& point_b,
        double h_intersect) const;

    //-- getter --//


    /**get_dtm_file*/
    std::string const& get_dtm_file() const noexcept {return m_dtm_file;};
    /**get_alt_data*/
    std::vector<double>  const&  get_alt_data() const noexcept {return m_alt_data;};
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
    std::vector<double> const& get_plane_coef_a() const noexcept {return m_plane_coef_a;};
    /**get_plane_coef_b*/
    std::vector<double> const& get_plane_coef_b() const noexcept {return m_plane_coef_b;};
    /**get_plane_coef_c*/
    std::vector<double> const& get_plane_coef_c() const noexcept {return m_plane_coef_c;};
    /**get_plane_coef_d*/
    std::vector<double> const& get_plane_coef_d() const noexcept {return m_plane_coef_d;};
    /**get_alt_min_cell*/
    double get_alt_min_cell() const noexcept {return m_alt_min_cell;};
    /**get_alt_max_cell*/
    double get_alt_max_cell() const noexcept {return m_alt_max_cell;};
    /**get_tol_z*/
    double get_tol_z() const noexcept {return m_tol_z;};// = 0.0001
    /**get_epsg*/
    int get_epsg() const noexcept {return m_epsg;};
    /**get_grid_row*/
    std::vector<double> const& get_grid_row() const noexcept {return m_grid_row;};
    /**get_grid_col*/
    std::vector<double> const& get_grid_col() const noexcept {return m_grid_col;};
    /**get_plans*/
    std::vector<double> const& get_plans() const noexcept {return m_plans;};
    /**get_trans_inv*/
    std::vector<double> const& get_trans_inv() const noexcept {return m_trans_inv;}; //affine.affine en python
    /**get_transform*/
    std::vector<double> const& get_transform() const noexcept {return m_transform;};
    /**get_nb_rows*/
    int get_nb_rows() const noexcept {return m_nb_rows;};
    /**get_nb_columns*/
    int get_nb_columns() const noexcept {return m_nb_columns;};

private:
    std::string m_dtm_file;
    std::vector<double> m_alt_data;
    double m_alt_min;
    double m_alt_max;
    double m_origin_x;
    double m_origin_y;
    double m_pixel_size_x;
    double m_pixel_size_y;
    std::vector<double> m_plane_coef_a;
    std::vector<double> m_plane_coef_b;
    std::vector<double> m_plane_coef_c;
    std::vector<double> m_plane_coef_d;
    double m_alt_min_cell;
    double m_alt_max_cell;
    double m_tol_z;// = 0.0001

    int m_epsg;

    std::vector<double> m_grid_row;
    std::vector<double> m_grid_col;

    std::vector<double> m_plans;

    std::vector<double> m_trans_inv; //affine.affine en python
    std::vector<double> m_transform;
    int m_nb_rows;
    int m_nb_columns;
};


#endif
