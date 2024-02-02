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
#include <array>
#include <algorithm>
#include <cmath>

#include <pybind11/pybind11.h>
#include "pybind11/numpy.h"


/**
Class DTMIntersection
Framework of the DTMIntersection python class.
*/

class DTMIntersection
{

public:

    /**Constructor*/
    DTMIntersection(
        int dtm_image_epsg,
        pybind11::array_t<double, \
                    pybind11::array::c_style | pybind11::array::forcecast> dtm_image_alt_data,
        int dtm_image_nb_rows,
        int dtm_image_nb_columns,
        std::tuple<double,double,double,double,double,double> dtm_image_transform
    );//determiner comment passer les arg

    /**eq_plan*/
    double eq_plan(int i, std::array<double, 3> const& position)const;

    /**ter_to_index*/
    std::array<double, 3> ter_to_index(std::array<double, 3> const& vect_ter)const;

    /**ter_to_indexs*/
    std::vector<double> ter_to_indexs(std::vector<double> const& vect_ter);//maybe unecessary

    /**index_to_ter*/
    std::array<double, 3> index_to_ter(std::array<double, 3> const& vect_ter)const;

    /**interpolate*/
    double interpolate(double delta_shift_row, double delta_shift_col) const;

    /**intersect_dtm_cube*/
    std::tuple<bool,
    bool,
    std::vector<double>,
    bool,
    std::vector<double>> intersect_dtm_cube(std::vector<double> const& los) const;

    /**intersection*/
    std::tuple<bool,bool,std::vector<double>> intersection(
        std::vector<double> const& los_index,
        std::vector<double> const& point_b,
        double h_intersect) const;

    //-- getter --//

    /**get_alt_data*/
    std::vector<double>  const&  get_alt_data() const noexcept {return m_alt_data;};
    /**get_alt_min*/
    double get_alt_min() const noexcept {return m_alt_min;};
    /**get_alt_max*/
    double get_alt_max() const noexcept {return m_alt_max;};
    /**get_plane_coef_a*/
    std::array<double,6> const& get_plane_coef_a() const noexcept {return m_plane_coef_a;};
    /**get_plane_coef_b*/
    std::array<double,6> const& get_plane_coef_b() const noexcept {return m_plane_coef_b;};
    /**get_plane_coef_c*/
    std::array<double,6> const& get_plane_coef_c() const noexcept {return m_plane_coef_c;};
    /**get_plane_coef_d*/
    std::array<double,6> const& get_plane_coef_d() const noexcept {return m_plane_coef_d;};
    /**get_alt_min_cell*/
    std::vector<double> get_alt_min_cell() const noexcept {return m_alt_min_cell;};
    /**get_alt_max_cell*/
    std::vector<double> get_alt_max_cell() const noexcept {return m_alt_max_cell;};
    /**get_tol_z*/
    double get_tol_z() const noexcept {return m_tol_z;};// = 0.0001
    /**get_epsg*/
    int get_epsg() const noexcept {return m_epsg;};
    /**get_plans*/
    std::vector<double> const& get_plans() const noexcept {return m_plans;};
    /**get trans_inv*/
    std::array<double,6> const& get_trans_inv() const noexcept {return m_trans_inv;};
    /**get_transform*/
    std::array<double,6> const& get_transform() const noexcept {return m_transform;};
    /**get_nb_rows*/
    int get_nb_rows() const noexcept {return m_nb_rows;};
    /**get_nb_columns*/
    int get_nb_columns() const noexcept {return m_nb_columns;};




private:

        /**alt_data attribut*/
        std::vector<double> m_alt_data;
        /**alt_min attribut*/
        double m_alt_min;
        /**alt_max attribut*/
        double m_alt_max;
        /**plane_coef_a attribut*/
        std::array<double,6> m_plane_coef_a;
        /**plane_coef_b attribut*/
        std::array<double,6> m_plane_coef_b;
        /**plane_coef_c attribut*/
        std::array<double,6> m_plane_coef_c;
        /**plane_coef_d attribut*/
        std::array<double,6> m_plane_coef_d;
        /**alt_min_cell attribut*/
        std::vector<double> m_alt_min_cell;
        /**alt_max_cell attribut*/
        std::vector<double> m_alt_max_cell;
        /**tol_z attribut*/
        double m_tol_z;// = 0.0001

        /**epsg attribut*/
        int m_epsg;

        /**plans attribut*/
        std::vector<double> m_plans;

        /**trans_inv attribut*/
        std::array<double,6> m_trans_inv; //affine.affine en python
        /**transform attribut*/
        std::array<double,6> m_transform;
        /**nb_rows attribut*/
        int m_nb_rows;
        /**nb_columns attribut*/
        int m_nb_columns;

};



//-- Function --//

/**init_min_max*/
std::tuple<std::vector<double>,
std::vector<double>> init_min_max(std::vector<double> const& alt_data,
                                                    int nb_rows,
                                                    int nb_columns);

#endif