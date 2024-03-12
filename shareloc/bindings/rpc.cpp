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

#include "rpc.hpp"

#include <stdexcept>
#include <iostream>
#include <cmath>

using namespace std;

//---- RPC methodes ----//

RPC::RPC(bool inverse_coefficient_input,
        bool direct_coefficient_input,
        array<double, 20> const& num_col_input,
        array<double, 20> const& den_col_input,
        array<double, 20> const& num_row_input,
        array<double, 20> const& den_row_input,
        array<double, 20> const& num_lon_input,
        array<double, 20> const& den_lon_input,
        array<double, 20> const& num_lat_input,
        array<double, 20> const& den_lat_input,
        array<double, 10> const& norm_coeffs){

    if(!inverse_coefficient_input && !direct_coefficient_input){
        throw runtime_error("C++: RPC: constructor: no RPC coeff");
    }

    m_inverse_coefficient = inverse_coefficient_input;
    m_direct_coefficient = direct_coefficient_input;

    if (inverse_coefficient_input){
        copy(num_col_input.begin(), num_col_input.end(), m_num_col.begin());
        copy(den_col_input.begin(), den_col_input.end(), m_den_col.begin());
        copy(num_row_input.begin(), num_row_input.end(), m_num_row.begin());
        copy(den_row_input.begin(), den_row_input.end(), m_den_row.begin());
    }
    if (direct_coefficient_input){
        copy(num_lon_input.begin(), num_lon_input.end(), m_num_lon.begin());
        copy(den_lon_input.begin(), den_lon_input.end(), m_den_lon.begin());
        copy(num_lat_input.begin(), num_lat_input.end(), m_num_lat.begin());
        copy(den_lat_input.begin(), den_lat_input.end(), m_den_lat.begin());
    }


    m_offset_lon = norm_coeffs[0];//offset_x
    m_scale_lon  = norm_coeffs[1];
    m_offset_lat = norm_coeffs[2];//offset_y
    m_scale_lat  = norm_coeffs[3];
    m_offset_alt = norm_coeffs[4];
    m_scale_alt  = norm_coeffs[5];
    m_offset_col = norm_coeffs[6];
    m_scale_col  = norm_coeffs[7];
    m_offset_row = norm_coeffs[8];
    m_scale_row  = norm_coeffs[9];

    m_lim_extrapol = 1.0001;

    m_alt_minmax = {m_offset_alt - m_scale_alt, m_offset_alt + m_scale_alt};
}

tuple<double,double,double> RPC::direct_loc_h(
    double row,
    double col,
    double alt,
    bool fill_nan,
    bool using_direct_coef) const
{
    if(!using_direct_coef && m_inverse_coefficient){
        return direct_loc_inverse_iterative(row, col, alt, 10, fill_nan);
    }else if(using_direct_coef && m_direct_coefficient){
        return compute_rational_function_polynomial_unitary(
            col,
            row,
            alt,
            m_num_lon,
            m_den_lon,
            m_num_lat,
            m_den_lat,
            fill_nan,
            "direct",
            m_scale_col,
            m_offset_col,
            m_scale_row,
            m_offset_row,
            m_scale_alt,
            m_offset_alt,
            m_scale_lon,
            m_offset_lon,
            m_scale_lat,
            m_offset_lat
        );
    }
    else{
        throw runtime_error("C++ : direct_loc_h: using_direct_coef doesn't\
         match with available coefficients");
    }
}

tuple<vector<double>,vector<double>,vector<double>> RPC::direct_loc_h(
    vector<double> const& row,
    vector<double> const& col,
    vector<double> const& alt,
    bool fill_nan,
    bool using_direct_coef) const
{

    if (!using_direct_coef && m_inverse_coefficient){
        return direct_loc_inverse_iterative(row, col, alt, 10, fill_nan);
    }
    else if(using_direct_coef && m_direct_coefficient){
        return compute_rational_function_polynomial(
            col,
            row,
            alt,
            m_num_lon,
            m_den_lon,
            m_num_lat,
            m_den_lat,
            fill_nan,
            "direct",
            m_scale_col,
            m_offset_col,
            m_scale_row,
            m_offset_row,
            m_scale_alt,
            m_offset_alt,
            m_scale_lon,
            m_offset_lon,
            m_scale_lat,
            m_offset_lat
        );
    }
    else{
        throw runtime_error("C++ : direct_loc_h: using_direct_coef doesn't\
         match with available coefficients");
    }
}

tuple<double,double,double> RPC::direct_loc_dtm(
    double row,
    double col,
    DTMIntersection const& dtm) const
{
    if(dtm.get_epsg() != 4326){
        throw runtime_error("C++ : direct_loc_dtm : epsg!=4326 -> Exiting");
        //cout<<"C++ : direct_loc_dtm : epsg!=4326 -> get_diff_alti_min/max initialized ?"<<endl;
    }

    double min_dtm = dtm.get_alt_min() - 1.0;
    double max_dtm = dtm.get_alt_max() + 1.0;
    
    vector<double> lon (2);
    vector<double> lat (2);
    vector<double> alt (2);
    bool var;//garbage variable
    bool solution;
    array<double,3> position_cube;
    double alti;
    vector<double> los_index_x (2);
    vector<double> los_index_y (2);
    vector<double> los_index_z (2);
    double position_x;
    double position_y;
    double position_z;

    tie(lon, lat, alt) = los_extrema(row, col, min_dtm, max_dtm);
    tie(solution, position_cube, alti, los_index_x, los_index_y, los_index_z) =\
    dtm.intersect_dtm_cube(lon, lat, alt);
    
    if(solution){
        tie(var, position_x, position_y, position_z) =\
        dtm.intersection(los_index_x, los_index_y, los_index_z, position_cube, alti);
    }
    else{
        position_x = numeric_limits<double>::quiet_NaN(); 
        position_y = numeric_limits<double>::quiet_NaN(); 
        position_z = numeric_limits<double>::quiet_NaN(); 
    }

    return {position_x,position_y,position_z};
}

tuple<vector<double>,vector<double>,vector<double>> RPC::direct_loc_dtm(
    vector<double> const& row,
    vector<double> const& col,
    DTMIntersection const& dtm) const
{
    if(dtm.get_epsg() != 4326){
        throw runtime_error("C++ : direct_loc_dtm : epsg!=4326 -> Exiting");
        //cout<<"C++ : direct_loc_dtm : epsg!=4326 -> get_diff_alti_min/max initialized ?"<<endl;
    }

    double min_dtm = dtm.get_alt_min() - 1.0;
    double max_dtm = dtm.get_alt_max() + 1.0;

    size_t nb_points = row.size();

    vector<double> lon (2);
    vector<double> lat (2);
    vector<double> alt (2);
    bool var;//garbage variable
    bool solution;
    array<double,3> position_cube;
    double alti;
    vector<double> los_index_x (2);
    vector<double> los_index_y (2);
    vector<double> los_index_z (2);
    double position_x;
    double position_y;
    double position_z;
    vector<double> res_lon (nb_points);
    vector<double> res_lat (nb_points);
    vector<double> res_alt (nb_points);

    for(size_t i = 0;i<nb_points;++i){

        tie(lon, lat, alt) = los_extrema(row[i], col[i], min_dtm, max_dtm);
        tie(solution, position_cube, alti, los_index_x, los_index_y, los_index_z) =\
        dtm.intersect_dtm_cube(lon, lat, alt);
        
        if(solution){
            tie(var, position_x, position_y, position_z) =\
            dtm.intersection(los_index_x, los_index_y, los_index_z, position_cube, alti);
        }
        else{
            position_x = numeric_limits<double>::quiet_NaN(); 
            position_y = numeric_limits<double>::quiet_NaN(); 
            position_z = numeric_limits<double>::quiet_NaN(); 
        }
        res_lon[i] = position_x;
        res_lat[i] = position_y;
        res_alt[i] = position_z;
    }
    
    return {res_lon,res_lat,res_alt};
}

tuple<double,double,double> RPC::inverse_loc(
    double lon,
    double lat,
    double alt)const
{

    auto [col_out, row_out, alt_out] = compute_rational_function_polynomial_unitary(
        lon,
        lat,
        alt,
        m_num_col,
        m_den_col,
        m_num_row,
        m_den_row,
        false,
        "inverse",
        m_scale_lon,
        m_offset_lon,
        m_scale_lat,
        m_offset_lat,
        m_scale_alt,
        m_offset_alt,
        m_scale_col,
        m_offset_col,
        m_scale_row,
        m_offset_row
    );

    return {row_out, col_out, alt_out};
}

tuple<vector<double>,vector<double>,vector<double>> RPC::inverse_loc(
    vector<double> const& lon,
    vector<double> const& lat,
    vector<double> const& alt) const
{

    auto [col_out, row_out, alt_res] = compute_rational_function_polynomial(
        lon,
        lat,
        alt,
        m_num_col,
        m_den_col,
        m_num_row,
        m_den_row,
        false,
        "inverse",
        m_scale_lon,
        m_offset_lon,
        m_scale_lat,
        m_offset_lat,
        m_scale_alt,
        m_offset_alt,
        m_scale_col,
        m_offset_col,
        m_scale_row,
        m_offset_row
    );
    return {row_out, col_out, alt_res};
}

tuple<double, double, double, double>
RPC::compute_loc_inverse_derivates(
    double lon,
    double lat,
    double alt)const
{
    //Normalisation
    double lon_norm = (lon - m_offset_lon)/m_scale_lon;
    double lat_norm = (lat - m_offset_lat)/m_scale_lat;
    double alt_norm = (alt - m_offset_alt)/m_scale_alt;

    alignas(64) array<double, 20> norms = pre_polynomial_equation(lon_norm, lat_norm, alt_norm);

    double num_dcol = polynomial_equation(norms, m_num_col);
    double den_dcol = polynomial_equation(norms, m_den_col);
    double num_drow = polynomial_equation(norms, m_num_row);
    double den_drow = polynomial_equation(norms, m_den_row);

    double num_dcol_dlon = derivative_polynomial_longitude(lon_norm, lat_norm, alt_norm, m_num_col);
    double den_dcol_dlon = derivative_polynomial_longitude(lon_norm, lat_norm, alt_norm, m_den_col);
    double num_drow_dlon = derivative_polynomial_longitude(lon_norm, lat_norm, alt_norm, m_num_row);
    double den_drow_dlon = derivative_polynomial_longitude(lon_norm, lat_norm, alt_norm, m_den_row);
                                         
    double num_dcol_dlat = derivative_polynomial_latitude(lon_norm, lat_norm, alt_norm, m_num_col);
    double den_dcol_dlat = derivative_polynomial_latitude(lon_norm, lat_norm, alt_norm, m_den_col);
    double num_drow_dlat = derivative_polynomial_latitude(lon_norm, lat_norm, alt_norm, m_num_row);
    double den_drow_dlat = derivative_polynomial_latitude(lon_norm, lat_norm, alt_norm, m_den_row);

    double dcol_dlon = m_scale_col / m_scale_lon * (num_dcol_dlon * den_dcol - den_dcol_dlon * num_dcol) / (den_dcol*den_dcol);
    double dcol_dlat = m_scale_col / m_scale_lat * (num_dcol_dlat * den_dcol - den_dcol_dlat * num_dcol) / (den_dcol*den_dcol);
    double drow_dlon = m_scale_row / m_scale_lon * (num_drow_dlon * den_drow - den_drow_dlon * num_drow) / (den_drow*den_drow);
    double drow_dlat = m_scale_row / m_scale_lat * (num_drow_dlat * den_drow - den_drow_dlat * num_drow) / (den_drow*den_drow);

    return {dcol_dlon, dcol_dlat, drow_dlon, drow_dlat};
}

tuple<double,double,double> RPC::direct_loc_inverse_iterative(
    double row,
    double col,
    double alt,
    int nb_iter_max,
    bool fill_nan)const
{

    double lon_out;
    double lat_out;

    // desired precision in pixels
    constexpr double eps = 1e-6;

    // Nan Filtering : if input nan -> output nan
    if(isnan(row) || isnan(col)){
        if(fill_nan){
            lon_out = m_offset_lon;
            lat_out = m_offset_lat;              
        }else{
            lon_out = numeric_limits<double>::quiet_NaN();
            lat_out = numeric_limits<double>::quiet_NaN();
        }
        return {lon_out, lat_out, alt};
    }
    else{
        lon_out = m_offset_lon;
        lat_out = m_offset_lat;
    }

    auto const [row_start, col_start, alt_start] = inverse_loc(lon_out, lat_out, alt);
    (void) alt_start; // disable "unused variable"

    // computing the residue between the sensor positions and those estimated
    //by the inverse localization
    double delta_col = col - col_start;
    double delta_row = row - row_start;


    // while the required precision is not achieved
    int iteration = 0;
    while ((abs(delta_col) > eps || abs(delta_row) > eps) && iteration < nb_iter_max){

        // partial derivatives
        auto const [dcol_dlon, dcol_dlat, drow_dlon, drow_dlat] = compute_loc_inverse_derivates(
            lon_out, lat_out, alt
        );


        double const det = dcol_dlon * drow_dlat - drow_dlon * dcol_dlat;

        double const delta_lon = (drow_dlat * delta_col - dcol_dlat * delta_row) / det;
        double const delta_lat = (-drow_dlon * delta_col + dcol_dlon * delta_row) / det;

        // update ground coordinates
        lon_out = lon_out + delta_lon;
        lat_out = lat_out + delta_lat;



        auto const [row_estim,col_estim,alt_estim] = inverse_loc(lon_out, lat_out, alt);
        (void) alt_estim; // disable "unused variable"

        // updating the residue between the sensor positions
        // and those estimated by the inverse localization
        delta_col = col - col_estim;
        delta_row = row - row_estim;

        ++iteration;
    }

    return {lon_out, lat_out, alt};
}

tuple<vector<double>,vector<double>,vector<double>> RPC::direct_loc_inverse_iterative(
    vector<double> const& row,
    vector<double> const& col,
    vector<double> const& alt,
    int nb_iter_max,
    bool fill_nan)const
{
    auto [col_norm,row_norm,alt_norm] = check_sizes(col,row,alt);

    // ** no suffix => _norm ** (same size) //

    size_t nb_points = col_norm.size();

    vector<double> lon_out(nb_points);
    vector<double> lat_out(nb_points);
    vector<bool> is_nan(nb_points);

    // desired precision in pixels
    constexpr double eps = 1e-6;

    // For all input point
    for (size_t i = 0;i<nb_points;++i){

        // Nan Filtering
        if(isnan(row_norm[i]) || isnan(col_norm[i])){
            if(fill_nan){
                lon_out[i] = m_offset_lon;
                lat_out[i] = m_offset_lat;              
            }else{
                lon_out[i] = numeric_limits<double>::quiet_NaN();
                lat_out[i] = numeric_limits<double>::quiet_NaN();
            }
            is_nan[i] = true;
            continue;
        }
        else{
            lon_out[i] = m_offset_lon;
            lat_out[i] = m_offset_lat;
            is_nan[i] = false;
        }


        // Initialisation
        auto const [row_start, col_start, alt_start] = inverse_loc(lon_out[i], lat_out[i], alt_norm[i]);
        (void) alt_start; // disable "unused variable"

        // computing the residue between the sensor positions and those estimated
        //by the inverse localization
        double delta_col = col_norm[i] - col_start;
        double delta_row = row_norm[i] - row_start;


        // while the required precision is not achieved
        int iteration = 0;
        while ((abs(delta_col) > eps || abs(delta_row) > eps) && iteration < nb_iter_max){

            // partial derivatives
            auto const [dcol_dlon, dcol_dlat, drow_dlon, drow_dlat] = compute_loc_inverse_derivates(
                lon_out[i], lat_out[i], alt_norm[i]
            );


            double const det = dcol_dlon * drow_dlat - drow_dlon * dcol_dlat;

            double const delta_lon = (drow_dlat * delta_col - dcol_dlat * delta_row) / det;
            double const delta_lat = (-drow_dlon * delta_col + dcol_dlon * delta_row) / det;

            // update ground coordinates
            lon_out[i] = lon_out[i]+delta_lon;
            lat_out[i] = lat_out[i]+delta_lat;



            auto const [row_estim,col_estim,alt_estim] = inverse_loc(lon_out[i], lat_out[i], alt_norm[i]);
            (void) alt_estim; // disable "unused variable"

            // updating the residue between the sensor positions
            // and those estimated by the inverse localization
            delta_col = col_norm[i] - col_estim;
            delta_row = row_norm[i] - row_estim;

            ++iteration;
        }
    }

    return {lon_out, lat_out, alt_norm};
}

array<double, 2> RPC::get_alt_min_max()const{
    array<double, 2> res = {m_offset_alt - m_scale_alt / 2.0, m_offset_alt + m_scale_alt / 2.0};
    return res;
}

tuple<vector<double>,vector<double>,vector<double>> RPC::los_extrema(
    double row,
    double col,
    double alt_min,
    double alt_max,
    bool fill_nan)const
{

    bool extrapolate = false;
    array<double, 2> los_alt_min_max;
    double los_alt_min = alt_min;
    double los_alt_max = alt_max;

    if (isnan(alt_min) || isnan(alt_max)){
        los_alt_min_max = get_alt_min_max();
        los_alt_min = los_alt_min_max[0];
        los_alt_max = los_alt_min_max[1];
    }else if(alt_min >= m_alt_minmax[0] && alt_max <= m_alt_minmax[1]){
        los_alt_min = alt_min;
        los_alt_max = alt_max;
    }else{
        extrapolate = true;
        los_alt_min_max = get_alt_min_max();
        los_alt_min = los_alt_min_max[0];
        los_alt_max = los_alt_min_max[1];
    }

    vector<double> row_array = {row,row};
    vector<double> col_array = {col,col};
    vector<double> alt_array = {los_alt_max,los_alt_min};


    auto [lon,lat,alt] = direct_loc_h(row_array, col_array, alt_array, fill_nan);

    if(extrapolate){
        double diff_lon = lon[0] - lon[1];
        double diff_lat = lat[0] - lat[1];
        double diff_alt = alt[0] - alt[1];

        double coeff_alt_max = (alt_max - alt[1]) / diff_alt;
        double coeff_alt_min = (alt_min - alt[1]) / diff_alt;

        lon[0] = lon[1] + diff_lon * coeff_alt_max;
        lat[0] = lat[1] + diff_lat * coeff_alt_max;
        alt[0] = alt[1] + diff_alt * coeff_alt_max;

        lon[1] = lon[1] + diff_lon * coeff_alt_min;
        lat[1] = lat[1] + diff_lat * coeff_alt_min;
        alt[1] = alt[1] + diff_alt * coeff_alt_min;
    }

    return {lon, lat, alt};
}

tuple<double,double,double>
RPC::compute_rational_function_polynomial_unitary(
    double lon_col,
    double lat_row,
    double alt,
    array<double, 20> const& num_col,
    array<double, 20> const& den_col,
    array<double, 20> const& num_lin,
    array<double, 20> const& den_lin,
    bool fill_nan,
    string direction,

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
) const
{
    double row_lat_out;
    double col_lon_out;
    double alt_out;
    if(isnan(lon_col) || isnan(lat_row)){

        if(fill_nan){
            if(direction=="direct"){
                col_lon_out = m_offset_lon;
                row_lat_out = m_offset_lat;
                alt_out = alt;
            }else{
                col_lon_out = m_offset_col;
                row_lat_out = m_offset_row;
                alt_out = alt;
            }
        }else{
            col_lon_out = numeric_limits<double>::quiet_NaN();
            row_lat_out = numeric_limits<double>::quiet_NaN();
            alt_out = alt;
        }
    }else{
        alt_out = alt;

        double lon_col_norm = (lon_col - offset_lon_col)/scale_lon_col;
        double lat_row_norm = (lat_row - offset_lat_row)/scale_lat_row;
        double alt_norm = (alt - offset_alt)/scale_alt;

        alignas(64) array<double, 20> norms = pre_polynomial_equation(lon_col_norm, lat_row_norm, alt_norm);
        double poly_num_col = polynomial_equation(norms, num_col);
        double poly_den_col = polynomial_equation(norms, den_col);
        double poly_num_lin = polynomial_equation(norms, num_lin);
        double poly_den_lin = polynomial_equation(norms, den_lin);


        col_lon_out = poly_num_col / poly_den_col * scale_col + offset_col;
        row_lat_out = poly_num_lin / poly_den_lin * scale_lin + offset_lin;
    };

    return {col_lon_out, row_lat_out, alt_out};
}

tuple<vector<double>,vector<double>,vector<double>> RPC::compute_rational_function_polynomial(
    vector<double> const& lon_col,
    vector<double> const& lat_row,
    vector<double> const& alt,
    array<double, 20> const& num_col,
    array<double, 20> const& den_col,
    array<double, 20> const& num_lin,
    array<double, 20> const& den_lin,
    bool fill_nan,
    string direction,

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
)const{
    auto [lon_col_norm,lat_row_norm,alt_norm] = check_sizes(lon_col,lat_row,alt);


    vector<double> col_lon_out(lon_col_norm.size());
    vector<double> row_lat_out(lon_col_norm.size());
    vector<double> alt_out(lon_col_norm.size());

    for(size_t i = 0;i<lon_col_norm.size();++i) {
        //--- Nan filtering
        if(isnan(lon_col_norm[i]) || isnan(lat_row_norm[i])){
            if(fill_nan){
                if(direction=="direct"){
                    col_lon_out[i] = m_offset_lon;
                    row_lat_out[i] = m_offset_lat;
                    alt_out[i] = alt[i];
                }else{
                    col_lon_out[i] = m_offset_col;
                    row_lat_out[i] = m_offset_row;
                    alt_out[i] = alt[i];
                }
            }else{
                col_lon_out[i] = numeric_limits<double>::quiet_NaN();
                row_lat_out[i] = numeric_limits<double>::quiet_NaN();
                alt_out[i] = alt[i];
            }
            continue;
        }

        alt_out[i] = alt_norm[i];

        //--- Normalisation
        lon_col_norm[i] = (lon_col_norm[i] - offset_lon_col)/scale_lon_col;
        lat_row_norm[i] = (lat_row_norm[i] - offset_lat_row)/scale_lat_row;
        alt_norm[i]     = (alt_norm[i] - offset_alt)/scale_alt;

        //-- Computation

        alignas(64) array<double, 20> norms = pre_polynomial_equation(lon_col_norm[i], lat_row_norm[i],alt_norm[i]);
        double poly_num_col = polynomial_equation(norms, num_col);
        double poly_den_col = polynomial_equation(norms, den_col);
        double poly_num_lin = polynomial_equation(norms, num_lin);
        double poly_den_lin = polynomial_equation(norms, den_lin);

        if (poly_den_col!=0 and poly_den_lin!=0){
            col_lon_out[i] = poly_num_col / poly_den_col * scale_col + offset_col;
            row_lat_out[i] = poly_num_lin / poly_den_lin * scale_lin + offset_lin;
        }
        else{
            throw runtime_error("C++ : compute_rational_function_polynomial: 0 divison");
        }
    }
    return {col_lon_out, row_lat_out, alt_out};
}

//
//
//
//
//
//---- Functions ----//
//
//
//
//
//


array<double, 20> pre_polynomial_equation(
        double xnorm,
        double ynorm,
        double znorm)
{
    return {
        1.0,
        xnorm,
        ynorm,
        znorm,
        xnorm * ynorm,
        xnorm * znorm,
        ynorm * znorm,
        xnorm * xnorm,
        ynorm * ynorm,
        znorm * znorm,
        xnorm * ynorm * znorm,
        xnorm * xnorm * xnorm,
        xnorm * ynorm * ynorm,
        xnorm * znorm * znorm,
        xnorm * xnorm * ynorm,
        ynorm * ynorm * ynorm,
        ynorm * znorm * znorm,
        xnorm * xnorm * znorm,
        ynorm * ynorm * znorm,
        znorm * znorm * znorm,
    };
}

double polynomial_equation(
        array<double, 20> const& norms,
        array<double, 20> const& coeffs)
{
    return
           norms[0]  * coeffs[0]
        +  norms[1]  * coeffs[1]
        +  norms[2]  * coeffs[2]
        +  norms[3]  * coeffs[3]
        +  norms[4]  * coeffs[4]
        +  norms[5]  * coeffs[5]
        +  norms[6]  * coeffs[6]
        +  norms[7]  * coeffs[7]
        +  norms[8]  * coeffs[8]
        +  norms[9]  * coeffs[9]
        +  norms[10] * coeffs[10]
        +  norms[11] * coeffs[11]
        +  norms[12] * coeffs[12]
        +  norms[13] * coeffs[13]
        +  norms[14] * coeffs[14]
        +  norms[15] * coeffs[15]
        +  norms[16] * coeffs[16]
        +  norms[17] * coeffs[17]
        +  norms[18] * coeffs[18]
        +  norms[19] * coeffs[19];
}

double derivative_polynomial_latitude(
    double lon_norm,
    double lat_norm,
    double alt_norm,
    const array<double, 20>& coeff)
{
    return
        coeff[2]
        + lon_norm * coeff[4]
        + alt_norm * coeff[6]
        + 2.0 * lat_norm * coeff[8]
        + lon_norm * alt_norm * coeff[10]
        + 2.0 * lat_norm * lon_norm * coeff[12]
        + lon_norm * lon_norm * coeff[14]
        + lat_norm * lat_norm * 3.0 * coeff[15]
        + alt_norm * alt_norm * coeff[16]
        + 2.0 * lat_norm * alt_norm * coeff[18];
}

double derivative_polynomial_longitude(
    double lon_norm,
    double lat_norm,
    double alt_norm,
    const array<double, 20>& coeff)
{
    return
        coeff[1]
        + lat_norm * coeff[4]
        + alt_norm * coeff[5]
        + 2.0 * lon_norm * coeff[7]
        + lat_norm * alt_norm * coeff[10]
        + lon_norm * lon_norm * 3.0 * coeff[11]
        + lat_norm * lat_norm * coeff[12]
        + alt_norm * alt_norm * coeff[13]
        + 2.0 * lon_norm * lat_norm * coeff[14]
        + 2.0 * lon_norm * alt_norm * coeff[17];
}




tuple<vector<double>, vector<double>, vector<double>>
check_sizes(
    vector<double> const& lon_col,
    vector<double> const& lat_row,
    vector<double> const& alt)
{

    vector<double> lon_col_norm;
    vector<double> lat_row_norm;
    vector<double> alt_norm;

    if(lat_row.size()<lon_col.size()){
        copy(lon_col.begin(), lon_col.begin()+lat_row.size(), back_inserter(lon_col_norm));
        lat_row_norm = lat_row;
    }else if (lat_row.size()>lon_col.size()){
        lon_col_norm = lon_col;
        copy(lat_row.begin(), lat_row.begin()+lon_col.size(), back_inserter(lat_row_norm));
    }else{
        lon_col_norm = lon_col;
        lat_row_norm = lat_row;
    };

    if (alt.size()!=lon_col_norm.size()){
        alt_norm.resize(lon_col_norm.size(),alt[0]);
    }else{
        alt_norm =alt;
    };

    return {lon_col_norm, lat_row_norm, alt_norm};
}
