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

#ifndef RPC_H
#define RPC_H

#include "dtm_intersection.hpp"

#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <array>
#include <algorithm>
#include <numeric>

/**
  Class RPC
  Framework of the RPC python class.
 */

class RPC
{
public:

    /**Constructor*/
    RPC(bool inverse_coefficient,
        bool direct_coefficient,
        std::array<double, 20> const& num_col,
        std::array<double, 20> const& den_col,
        std::array<double, 20> const& num_row,
        std::array<double, 20> const& den_row,
        std::array<double, 20> const& num_lon,
        std::array<double, 20> const& den_lon,
        std::array<double, 20> const& num_lat,
        std::array<double, 20> const& den_lat,
        std::array<double, 10> const& norm_coeffs);

    /**direct_loc_h*/
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> direct_loc_h(
        std::vector<double> const& row,
        std::vector<double> const& col,
        std::vector<double> const& alt,
        bool fill_nan=false) const;

    /**direct_loc_grid_h*/
    std::tuple<std::vector<std::vector<double>>,std::vector<std::vector<double>>> direct_loc_grid_h(
        int row0,
        int col0,
        int steprow,
        int stepcol,
        int nbrow,
        int nbcol,
        double alt) const;

    /**direct_loc_dtm*/
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> direct_loc_dtm(
        std::vector<double> const& row,
        std::vector<double> const& col,
        DTMIntersection dtm) const;//override + dtm is a python class not a std::string

    /**inverse_loc unitary*/
    std::tuple<double,double,double> inverse_loc(
        double lon,
        double lat,
        double alt)const;

    /**inverse_loc*/
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>> inverse_loc(
        std::vector<double> const& lon,
        std::vector<double> const& lat,
        std::vector<double> const& alt)const;// override;

    /**filter_coordinates*/
    std::tuple<std::vector<bool>,std::vector<double>,std::vector<double>> filter_coordinates(
        std::vector<double> const& first_coord,
        std::vector<double> const& second_coord,
        bool fill_nan=false,
        std::string direction="direct")const;


    /**compute_loc_inverse_derivates unitary*/
    std::tuple<double, double, double, double> compute_loc_inverse_derivates(
        double lon,
        double lat,
        double alt) const;

    /**direct_loc_inverse_iterative*/
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
    direct_loc_inverse_iterative(
        std::vector<double> const& row,
        std::vector<double> const& col,
        std::vector<double> const& alt,
        int nb_iter_max=10,
        bool fill_nan=false)const;

    /**get_alt_min_max*/
    std::array<double, 2> get_alt_min_max()const;

    /**los_extrema*/
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>
    los_extrema(
        double row,
        double col,
        double alt_min,
        double alt_max,
        bool fill_nan=false) const;


    //-- getter --//

    /**get_num_col*/
    std::array<double, 20> const& get_num_col() const noexcept {return m_num_col;};
    /**get_den_col*/
    std::array<double, 20> const& get_den_col() const noexcept {return m_den_col;};
    /**get_num_row*/
    std::array<double, 20> const& get_num_row() const noexcept {return m_num_row;};
    /**get_den_row*/
    std::array<double, 20> const& get_den_row() const noexcept {return m_den_row;};

    /**get_num_lon*/
    std::array<double, 20> const& get_num_lon() const noexcept {return m_num_lon;};
    /**get_den_lon*/
    std::array<double, 20> const& get_den_lon() const noexcept {return m_den_lon;};
    /**get_num_lat*/
    std::array<double, 20> const& get_num_lat() const noexcept {return m_num_lat;};
    /**get_den_lat*/
    std::array<double, 20> const& get_den_lat() const noexcept {return m_den_lat;};

    /**get_offset_row*/
    double get_offset_row() const noexcept {return m_offset_row;};
    /**get_scale_row*/
    double get_scale_row() const noexcept {return m_scale_row;};
    /**get_offset_col*/
    double get_offset_col() const noexcept {return m_offset_col;};
    /**get_scale_col*/
    double get_scale_col() const noexcept {return m_scale_col;};
    /**get_offset_alt*/
    double get_offset_alt() const noexcept {return m_offset_alt;};
    /**get_scale_alt*/
    double get_scale_alt() const noexcept {return m_scale_alt;};
    /**get_offset_lon*/
    double get_offset_lon() const noexcept {return m_offset_lon;};
    /**get_scale_lon*/
    double get_scale_lon() const noexcept {return m_scale_lon;};
    /**get_offset_lat*/
    double get_offset_lat() const noexcept {return m_offset_lat;};
    /**get_scale_lat*/
    double get_scale_lat() const noexcept {return m_scale_lat;};

private:

    std::string m_datum;
    std::map<std::string, double> m_rpc_params;//it's a simple dict -> maybe inappropiate
    double m_lim_extrapol;

    std::vector<std::vector<double>> m_monomes;// to convert to array
    std::vector<std::vector<double>> m_monomes_deriv_1;
    std::vector<std::vector<double>> m_monomes_deriv_2;

    bool m_inverse_coefficient;
    bool m_direct_coefficient;

    alignas(64) std::array<double, 20> m_num_col;
    alignas(64) std::array<double, 20> m_den_col;
    alignas(64) std::array<double, 20> m_num_row;
    alignas(64) std::array<double, 20> m_den_row;

    alignas(64) std::array<double, 20> m_num_lon;
    alignas(64) std::array<double, 20> m_den_lon;
    alignas(64) std::array<double, 20> m_num_lat;
    alignas(64) std::array<double, 20> m_den_lat;

    std::array<double, 2> m_alt_minmax;

    double m_col0;
    double m_colmax;
    double m_row0;
    double m_rowmax;

    double m_offset_row;
    double m_scale_row;
    double m_offset_col;
    double m_scale_col;
    double m_offset_alt;
    double m_scale_alt;
    double m_offset_lon;//py: lon=x
    double m_scale_lon;//py: lon=x
    double m_offset_lat;//py: lat=y
    double m_scale_lat;//py: lat=y

};

// function

/**Compute pre_polynomial_equation*/
std::array<double, 20> pre_polynomial_equation(
        double xnorm,
        double ynorm,
        double znorm);

/**Compute polynomial equation"*/
double polynomial_equation(
        std::array<double, 20> const& norms,
        std::array<double, 20> const& coeffs);

/** compute_rational_function_polynomial unitary*/
std::tuple<double,double,double> compute_rational_function_polynomial_unitary(
    double lon_col_norm,
    double lat_row_norm,
    double alt_norm,
    std::array<double, 20> const& num_col,
    std::array<double, 20> const& den_col,
    std::array<double, 20> const& num_lin,
    std::array<double, 20> const& den_lin,

    //input
    double scale_lon_col,
    double offset_lon_col,
    double scale_lat_row,
    double offset_lat_row,
    double scale_alt,
    double offset_alt,

    //output
    double scale_col,
    double offset_col,
    double scale_lin,
    double offset_lin
);


/**Compute rational function polynomial. Useful to compute direct and inverse localization
        "using direct or inverse RPC."*/
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
compute_rational_function_polynomial(
    std::vector<double> const& lon_col_norm,
    std::vector<double> const& lat_row_norm,
    std::vector<double> const& alt_norm,
    std::array<double, 20> const& num_col,
    std::array<double, 20> const& den_col,
    std::array<double, 20> const& num_lin,
    std::array<double, 20> const& den_lin,

    //input
    double scale_lon_col,
    double offset_lon_col,
    double scale_lat_row,
    double offset_lat_row,
    double scale_alt,
    double offset_alt,

    //output
    double scale_col,
    double offset_col,
    double scale_lin,
    double offset_lin
);


/** Compute derivative_polynomial_latitude*/
double derivative_polynomial_latitude(
    double xnorm,
    double ynorm,
    double znorm,
    std::array<double, 20> const& coeff);


/**Compute derivative_polynomial_longitude*/
double derivative_polynomial_longitude(
    double xnorm,
    double ynorm,
    double znorm,
    std::array<double, 20> const& coeff);


/**Check if arrays have the same size and cut it if needed*/
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
check_sizes(
        std::vector<double> const& lon_col,
        std::vector<double> const& lat_row,
        std::vector<double> const& alt);

#endif
