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
#include "rpc.hpp"

//---- RPC methodes ----//

RPC::RPC(bool inverse_coefficient_input,
        bool direct_coefficient_input,
        array<double, 20> num_col_input,
        array<double, 20> den_col_input,
        array<double, 20> num_row_input,
        array<double, 20> den_row_input,
        array<double, 20> num_lon_input,
        array<double, 20> den_lon_input,
        array<double, 20> num_lat_input,
        array<double, 20> den_lat_input,
        array<double, 10> norm_coeffs):GeoModelTemplate(){

    if(!inverse_coefficient && !direct_coefficient){
        cerr<<"C++: RPC: constructor: no RPC coeff"<<endl;
    }

    inverse_coefficient = inverse_coefficient_input;
    direct_coefficient = direct_coefficient_input;

    if (inverse_coefficient){
        std::copy(num_col_input.begin(), num_col_input.end(), num_col.begin());
        std::copy(den_col_input.begin(), den_col_input.end(), den_col.begin());
        std::copy(num_row_input.begin(), num_row_input.end(), num_row.begin());
        std::copy(den_row_input.begin(), den_row_input.end(), den_row.begin());
    }
    if (direct_coefficient){
        std::copy(num_lon_input.begin(), num_lon_input.end(), num_lon.begin());
        std::copy(den_lon_input.begin(), den_lon_input.end(), den_lon.begin());
        std::copy(num_lat_input.begin(), num_lat_input.end(), num_lat.begin());
        std::copy(den_lat_input.begin(), den_lat_input.end(), den_lat.begin());
    }


    offset_lon = norm_coeffs[0];//offset_x
    scale_lon = norm_coeffs[1];
    offset_lat = norm_coeffs[2];//offset_y
    scale_lat = norm_coeffs[3];
    offset_alt = norm_coeffs[4];
    scale_alt = norm_coeffs[5];
    offset_col = norm_coeffs[6];
    scale_col = norm_coeffs[7];
    offset_row = norm_coeffs[8];
    scale_row = norm_coeffs[9];

    lim_extrapol = 1.0001;

    alt_minmax = {offset_alt - scale_alt, offset_alt + scale_alt};
}

tuple<vector<double>,vector<double>,vector<double>> RPC::direct_loc_h(
    vector<double> row,
    vector<double> col,
    vector<double> alt,
    bool fill_nan){

    vector<double> lon_out;
    vector<double> lat_out;
    vector<double> alt_out;
    if (direct_coefficient){
        // if(abs(col_norm[i])>lim_extrapol || 
        //     abs(row_norm[i])>lim_extrapol ||
        //     abs(alt_norm[i])>lim_extrapol){
        //     cout<<"Warning : normalisation values exceed lim_extrapol"<<endl;}

        tie(lon_out, lat_out, alt_out) = compute_rational_function_polynomial(
            col,
            row,
            alt,
            num_lon,
            den_lon,
            num_lat,
            den_lat,
            scale_col,
            offset_col,
            scale_row,
            offset_row,
            scale_alt,
            offset_alt,
            scale_lon,
            offset_lon,
            scale_lat,
            offset_lat
        );

    }else{
        tie(lon_out, lat_out, alt_out)  = direct_loc_inverse_iterative(row, col, alt, 10, fill_nan);

    }

    return make_tuple(lon_out, lat_out, alt_out);
}

tuple<vector<vector<double>>,vector<vector<double>>> RPC::direct_loc_grid_h(
    int row0,
    int col0,
    int steprow,
    int stepcol,
    int nbrow,
    int nbcol,
    double alt){
    tuple<vector<vector<double>>,vector<vector<double>>> res;
    return res;
}

vector<vector<double>> RPC::direct_loc_dtm(
    double row,
    double col,
    string dtm){// dtm intersection model ->python class
    vector<vector<double>> vect;
    return vect;
}

tuple<double,double,double> RPC::inverse_loc(
    double lon,
    double lat,
    double alt)
    /**
    Inverse localization

    :param lon: longitude position
    :type lon: double
    :param lat: latitude position
    :type lat: double
    :param alt: altitude
    :type alt: double
    :return: sensor position (row, col, alt)
    :rtype: tuple(double row position, double col position, double alt)
    */
    //for further info check inverse loc vectorized
{
    // if(abs(lon_norm[i])>lim_extrapol || 
    //     abs(lat_norm[i])>lim_extrapol ||
    //     abs(alt_norm[i])>lim_extrapol){
    //     //cout<<"Warning : normalisation values exceed lim_extrapol"<<endl;

    
    double row_out;
    double col_out;
    double alt_out;
    tie(col_out, row_out, alt_out) = compute_rational_function_polynomial_unitary(
        lon,
        lat,
        alt,
        num_col,
        den_col,
        num_row,
        den_row,
        scale_lon,
        offset_lon,
        scale_lat,
        offset_lat,
        scale_alt,
        offset_alt,
        scale_col,
        offset_col,
        scale_row,
        offset_row
    );
    return make_tuple(row_out, col_out, alt_out);
}

tuple<vector<double>,vector<double>,vector<double>> RPC::inverse_loc(
    vector<double> lon,
    vector<double> lat,
    vector<double> alt)
    /**
    Inverse localization

    :param lon: longitude position
    :type lon: float or 1D numpy.ndarray dtype=float64
    :param lat: latitude position
    :type lat: float or 1D numpy.ndarray dtype=float64
    :param alt: altitude
    :type alt: float
    :return: sensor position (row, col, alt)
    :rtype: tuple(1D np.array row position, 1D np.array col position, 1D np.array alt)
    */
    
    // TODO : Nan filtering

{
    //cout<<"C++ : Computing inverse loc"<<endl;

    //-- Useless c++ compiler return error if arg not good type (kept just in case)
    // if (!is_same<decltype(lon), vector<double>>::value){//if lon not a vector<double>
    //     try{
    //         lon = vector<double>{lon};
    //         lat = vector<double>{lat};
    //     }catch(...){
    //      cerr<<"Error : longitude/latitude not vector<double> and can't be converted"<<endl;
    //     }
    //     };
    // if (!is_same<decltype(alt), vector<double>>::value){//if alt not a vector<double>
    //     try{
    //         alt = vector<double>{alt};
    //     }catch(...){
    //         cerr<<"Error : altitude not vector<double> and can't be converted"<<endl;
    //     }
    //     };


    // if(abs(lon_norm[i])>lim_extrapol || 
    //     abs(lat_norm[i])>lim_extrapol ||
    //     abs(alt_norm[i])>lim_extrapol){
    //     //cout<<"Warning : normalisation values exceed lim_extrapol"<<endl;
    // }
  
    vector<double> row_out;
    vector<double> col_out;
    vector<double> alt_res;
    tie(col_out, row_out, alt_res) = compute_rational_function_polynomial(
        lon,
        lat,
        alt,
        num_col,
        den_col,
        num_row,
        den_row,
        scale_lon,
        offset_lon,
        scale_lat,
        offset_lat,
        scale_alt,
        offset_alt,
        scale_col,
        offset_col,
        scale_row,
        offset_row
    );
    return make_tuple(row_out, col_out, alt_res);
}


tuple<vector<bool>,vector<double>,vector<double>> RPC::filter_coordinates(// May be useless
    vector<double> first_coord,
    vector<double> second_coord,
    bool fill_nan,
    string direction){

        // assert none problematic inputs
        if (first_coord.size()!=second_coord.size()){
            cerr<<"Error : C++ : filter_coordinates : ";
            cerr<<"first_coord.size()!=second_coord.size()"<<endl;
            exit(EXIT_FAILURE);
        }

        // filter_nan computation 
        size_t size = second_coord.size();
        vector<bool> filter_nan(size);
        for (int i = 0; i < int(size); ++i) {
            filter_nan[i] = !(isnan(first_coord[i]) || isnan(second_coord[i]));
        }

        // Output computation
        vector<double>x_out(size,numeric_limits<double>::quiet_NaN());
        vector<double>y_out(size,numeric_limits<double>::quiet_NaN());   
        if (fill_nan){
            double out_x_nan_value;
            double out_y_nan_value;
            if (direction=="direct"){
                out_x_nan_value = offset_lon;
                out_y_nan_value = offset_lat;
            }else{
                out_x_nan_value = offset_col;
                out_y_nan_value = offset_row;              
            }
            fill(x_out.begin(), x_out.end(), out_x_nan_value);
            fill(y_out.begin(), y_out.end(), out_y_nan_value);
        }
    return make_tuple(filter_nan, x_out, y_out);
}

tuple<double,
double,
double,
double> RPC::compute_loc_inverse_derivates(
    double lon,
    double lat,
    double alt){


    //Normalisation

    double lon_norm = (lon - offset_lon)/scale_lon;
    double lat_norm = (lat - offset_lat)/scale_lat;
    double alt_norm = (alt - offset_alt)/scale_alt;

    double num_dcol = polynomial_equation(lon_norm, lat_norm, alt_norm, &num_col);
    double den_dcol = polynomial_equation(lon_norm, lat_norm, alt_norm, &den_col);
    double num_drow = polynomial_equation(lon_norm, lat_norm, alt_norm, &num_row);
    double den_drow = polynomial_equation(lon_norm, lat_norm, alt_norm, &den_row);

    double num_dcol_dlon = derivative_polynomial_longitude(lon_norm, lat_norm, alt_norm, &num_col);
    double den_dcol_dlon = derivative_polynomial_longitude(lon_norm, lat_norm, alt_norm, &den_col);
    double num_drow_dlon = derivative_polynomial_longitude(lon_norm, lat_norm, alt_norm, &num_row);
    double den_drow_dlon = derivative_polynomial_longitude(lon_norm, lat_norm, alt_norm, &den_row);

    double num_dcol_dlat = derivative_polynomial_latitude(lon_norm, lat_norm, alt_norm, &num_col);
    double den_dcol_dlat = derivative_polynomial_latitude(lon_norm, lat_norm, alt_norm, &den_col);
    double num_drow_dlat = derivative_polynomial_latitude(lon_norm, lat_norm, alt_norm, &num_row);
    double den_drow_dlat = derivative_polynomial_latitude(lon_norm, lat_norm, alt_norm, &den_row);

    double dcol_dlon = scale_col / scale_lon * (num_dcol_dlon * den_dcol - den_dcol_dlon * num_dcol) / (den_dcol*den_dcol);
    double dcol_dlat = scale_col / scale_lat * (num_dcol_dlat * den_dcol - den_dcol_dlat * num_dcol) / (den_dcol*den_dcol);
    double drow_dlon = scale_row / scale_lon * (num_drow_dlon * den_drow - den_drow_dlon * num_drow) / (den_drow*den_drow);
    double drow_dlat = scale_row / scale_lat * (num_drow_dlat * den_drow - den_drow_dlat * num_drow) / (den_drow*den_drow);


    return make_tuple(dcol_dlon, dcol_dlat, drow_dlon, drow_dlat);
    }



tuple<vector<double>,vector<double>,vector<double>> RPC::direct_loc_inverse_iterative(
    vector<double> row,
    vector<double> col,
    vector<double> alt,
    int nb_iter_max,
    bool fill_nan){

    vector<double> col_norm;
    vector<double> row_norm;
    vector<double> alt_norm;

    tie(col_norm,row_norm,alt_norm) = check_sizes(col,row,alt);

    // ** no suffix => _norm ** (same size) //

    size_t nb_point = col_norm.size();

    vector<double> lon_out(nb_point);//array<double, nb_point> => ‘size_t nb_point’ is not const
    vector<double> lat_out(nb_point);
    vector<bool> is_nan(nb_point);


    // desired precision in pixels
    double eps = 1e-6;

    double row_start;
    double col_start;
    double alt_start;

    int iteration;

    double delta_col;
    double delta_row;

    double drow_dlat;
    double drow_dlon;
    double dcol_dlat;
    double dcol_dlon;

    double det;
    double delta_lon;
    double delta_lat;

    double row_estim;
    double col_estim;
    double alt_estim;

    // For all input point
    for (int i = 0;i<int(nb_point);++i){


        // Nan Filtering : if input nan -> output nan
        if(isnan(row_norm[i]) || isnan(col_norm[i])){
            lon_out[i] = numeric_limits<double>::quiet_NaN();
            lat_out[i] = numeric_limits<double>::quiet_NaN();
            is_nan[i] = true;
            continue;
        }
        else{
            lon_out[i] = offset_lon;
            lat_out[i] = offset_lat;
            is_nan[i] = false;
        }


        // Initialisation /!\ CAN BE OUTSIDE THE LOOP but  alt_norm[i]
        tie(row_start, col_start, alt_start) = inverse_loc(lon_out[i], lat_out[i], alt_norm[i]);


        iteration = 0;
        // computing the residue between the sensor positions and those estimated
        //by the inverse localization
        delta_col = col_norm[i] - col_start;
        delta_row = row_norm[i] - row_start;


        // while the required precision is not achieved
        while ((fabs(delta_col) > eps || fabs(delta_row) > eps) && iteration < nb_iter_max){


            // partial derivatives
            tie(dcol_dlon, dcol_dlat, drow_dlon, drow_dlat) = compute_loc_inverse_derivates(
                lon_out[i], lat_out[i], alt_norm[i]
            );

            
            det = dcol_dlon * drow_dlat - drow_dlon * dcol_dlat;

            delta_lon = (drow_dlat * delta_col - dcol_dlat * delta_row) / det;
            delta_lat = (-drow_dlon * delta_col + dcol_dlon * delta_row) / det;

            // update ground coordinates 
            lon_out[i] = lon_out[i]+delta_lon;
            lat_out[i] = lat_out[i]+delta_lat;



            tie(row_estim,col_estim,alt_estim) = inverse_loc(lon_out[i], lat_out[i], alt_norm[i]);
            //5e-11 = maxerror w/r to python and alt_estim useless

            // updating the residue between the sensor positions
            // and those estimated by the inverse localization
            delta_col = col_norm[i] - col_estim;
            delta_row = row_norm[i] - row_estim;

            ++iteration;
        }
    }
    
    return make_tuple(lon_out, lat_out, alt_norm);
}

array<double, 2> RPC::get_alt_min_max(){
    array<double, 2> res = {offset_alt - scale_alt / 2.0, offset_alt + scale_alt / 2.0};
    return res;
}

tuple<vector<double>,vector<double>,vector<double>> RPC::los_extrema(
    double row,
    double col,
    double alt_min,
    double alt_max,
    bool fill_nan,
    int epsg){

    if(epsg!=4326){
         throw runtime_error("C++ : los_extrema : epsg!=4326 -> Exiting");
    }

    bool extrapolate = false;
    array<double, 2> los_alt_min_max;
    double los_alt_min = alt_min;
    double los_alt_max = alt_max;

    if (isnan(alt_min) || isnan(alt_max)){
        los_alt_min_max = get_alt_min_max();
        los_alt_min = los_alt_min_max[0];
        los_alt_max = los_alt_min_max[1];
    }else if(alt_min >= alt_minmax[0] && alt_max <= alt_minmax[1]){
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


    vector<double> lon;
    vector<double> lat;
    vector<double> alt;

    tie(lon,lat,alt) = direct_loc_h(row_array, col_array, alt_array, fill_nan);


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

    return make_tuple(lon, lat, alt);
}

array<double, 20> RPC::get_num_col(){return num_col;}
array<double, 20> RPC::get_den_col(){return den_col;}
array<double, 20> RPC::get_num_row(){return num_row;}
array<double, 20> RPC::get_den_row(){return den_row;}

array<double, 20> RPC::get_num_lon(){return num_lon;}
array<double, 20> RPC::get_den_lon(){return den_lon;}
array<double, 20> RPC::get_num_lat(){return num_lat;}
array<double, 20> RPC::get_den_lat(){return den_lat;}

double RPC::get_offset_row(){return offset_row;}
double RPC::get_scale_row(){return scale_row;}
double RPC::get_offset_col(){return offset_col;}
double RPC::get_scale_col(){return scale_col;}
double RPC::get_offset_alt(){return offset_alt;}
double RPC::get_scale_alt(){return scale_alt;}
double RPC::get_offset_lon(){return offset_lon;}
double RPC::get_scale_lon(){return scale_lon;}
double RPC::get_offset_lat(){return offset_lat;}
double RPC::get_scale_lat(){return scale_lat;}



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




double polynomial_equation(
    double xnorm,
    double ynorm,
    double znorm,
    const array<double, 20>* coeff){
    
    return
    (*coeff)[0]
    + (*coeff)[1] * xnorm
    + (*coeff)[2] * ynorm
    + (*coeff)[3] * znorm
    + (*coeff)[4] * xnorm * ynorm
    + (*coeff)[5] * xnorm * znorm
    + (*coeff)[6] * ynorm * znorm
    + (*coeff)[7] * xnorm*xnorm
    + (*coeff)[8] * ynorm*ynorm
    + (*coeff)[9] * znorm*znorm
    + (*coeff)[10] * xnorm * ynorm * znorm
    + (*coeff)[11] * xnorm*xnorm*xnorm
    + (*coeff)[12] * xnorm * ynorm*ynorm
    + (*coeff)[13] * xnorm * znorm*znorm
    + (*coeff)[14] * xnorm*xnorm * ynorm
    + (*coeff)[15] * ynorm*ynorm*ynorm
    + (*coeff)[16] * ynorm * znorm*znorm
    + (*coeff)[17] * xnorm*xnorm * znorm
    + (*coeff)[18] * ynorm*ynorm * znorm
    + (*coeff)[19] * znorm*znorm*znorm;
}

tuple<double,double,double> compute_rational_function_polynomial_unitary(
    double lon_col,
    double lat_row,
    double alt,
    array<double, 20> num_col,
    array<double, 20> den_col,
    array<double, 20> num_lin,
    array<double, 20> den_lin,

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
){
    double row_lat_out;
    double col_lon_out;
    double alt_out;
    if(isnan(lon_col) || isnan(lat_row)){
        col_lon_out = numeric_limits<double>::quiet_NaN();
        row_lat_out = numeric_limits<double>::quiet_NaN();
        alt_out = numeric_limits<double>::quiet_NaN();
    }else{
        alt_out = alt;

        double lon_col_norm = (lon_col - offset_lon_col)/scale_lon_col;
        double lat_row_norm = (lat_row - offset_lat_row)/scale_lat_row;
        double alt_norm = (alt - offset_alt)/scale_alt;

        double poly_num_col = polynomial_equation(lon_col_norm, lat_row_norm, alt_norm, &num_col);
        double poly_den_col = polynomial_equation(lon_col_norm, lat_row_norm, alt_norm, &den_col);
        double poly_num_lin = polynomial_equation(lon_col_norm, lat_row_norm, alt_norm, &num_lin);
        double poly_den_lin = polynomial_equation(lon_col_norm, lat_row_norm, alt_norm, &den_lin);
        col_lon_out = poly_num_col / poly_den_col * scale_col + offset_col;
        row_lat_out = poly_num_lin / poly_den_lin * scale_lin + offset_lin;
    };
    return make_tuple(col_lon_out, row_lat_out, alt_out);
}

tuple<vector<double>,vector<double>,vector<double>> compute_rational_function_polynomial(
    vector<double> lon_col,
    vector<double> lat_row,
    vector<double> alt,
    array<double, 20> num_col,
    array<double, 20> den_col,
    array<double, 20> num_lin,
    array<double, 20> den_lin,

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
){

    vector<double> lon_col_norm;
    vector<double> lat_row_norm;
    vector<double> alt_norm;

    tie(lon_col_norm,lat_row_norm,alt_norm) = check_sizes(lon_col,lat_row,alt);


    vector<double> col_lon_out(lon_col_norm.size());
    vector<double> row_lat_out(lon_col_norm.size());
    vector<double> alt_out(lon_col_norm.size());

    double poly_num_col;
    double poly_den_col;
    double poly_num_lin;
    double poly_den_lin;
    for(int i = 0;i<(int)lon_col_norm.size();++i){


        //--- Nan filtering
        if(isnan(lon_col_norm[i]) || isnan(lat_row_norm[i])){
            col_lon_out[i] = numeric_limits<double>::quiet_NaN();
            row_lat_out[i] = numeric_limits<double>::quiet_NaN();
            alt_out[i] = numeric_limits<double>::quiet_NaN();
            continue;
        }

        alt_out[i] = alt_norm[i];

        //--- Normalisation
        lon_col_norm[i] = (lon_col_norm[i] - offset_lon_col)/scale_lon_col;
        lat_row_norm[i] = (lat_row_norm[i] - offset_lat_row)/scale_lat_row;
        alt_norm[i] = (alt_norm[i] - offset_alt)/scale_alt;




        //-- Computation
        poly_num_col = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], &num_col);
        poly_den_col = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], &den_col);
        poly_num_lin = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], &num_lin);
        poly_den_lin = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], &den_lin);

        if (poly_den_col!=0 and poly_den_lin!=0){
            col_lon_out[i] = poly_num_col / poly_den_col * scale_col + offset_col;
            row_lat_out[i] = poly_num_lin / poly_den_lin * scale_lin + offset_lin;
        }
        else{
            cerr<<"C++ : compute_rational_function_polynomial: 0 divison"<<endl;
        }
    }
    tuple<vector<double>,vector<double>,vector<double>> res = make_tuple(col_lon_out, row_lat_out, alt_out);
    return res;
}

double derivative_polynomial_latitude(
    double lon_norm,
    double lat_norm,
    double alt_norm,
    const array<double, 20>* coeff){//array<double, 20> coeff){

    return
        (*coeff)[2]
        + (*coeff)[4] * lon_norm
        + (*coeff)[6] * alt_norm
        + 2 * (*coeff)[8] * lat_norm
        + (*coeff)[10] * lon_norm * alt_norm
        + 2 * (*coeff)[12] * lon_norm * lat_norm
        + (*coeff)[14] * lon_norm*lon_norm
        + 3 * (*coeff)[15] * lat_norm*lat_norm
        + (*coeff)[16] * alt_norm*alt_norm
        + 2 * (*coeff)[18] * lat_norm * alt_norm;
}

double derivative_polynomial_longitude(
    double lon_norm,
    double lat_norm,
    double alt_norm,
    const array<double, 20>* coeff){//array<double, 20> coeff){
    return
        (*coeff)[1]
        + (*coeff)[4] * lat_norm
        + (*coeff)[5] * alt_norm
        + 2 * (*coeff)[7] * lon_norm
        + (*coeff)[10] * lat_norm * alt_norm
        + 3 * (*coeff)[11] * lon_norm*lon_norm
        + (*coeff)[12] * lat_norm*lat_norm
        + (*coeff)[13] * alt_norm*alt_norm
        + 2 * (*coeff)[14] * lat_norm * lon_norm
        + 2 * (*coeff)[17] * lon_norm * alt_norm;
}


tuple<vector<double>,
    vector<double>,
    vector<double>> check_sizes(vector<double> lon_col,
    vector<double> lat_row,
    vector<double>alt){

    vector<double> lon_col_norm;
    vector<double> lat_row_norm;
    vector<double> alt_norm;

    if(lat_row.size()<lon_col.size()){
        cout<<" lat_row.size()!=lon_col.size() -> truncate lon_col"<<endl;
        copy(lon_col.begin(), lon_col.begin()+lat_row.size(), back_inserter(lon_col_norm));
        lat_row_norm = lat_row;
    }else if (lat_row.size()>lon_col.size()){
        cout<<" lat_row.size()!=lon_col.size() -> truncate lat_row "<<endl;
        lon_col_norm = lon_col;
        copy(lat_row.begin(), lat_row.begin()+lon_col.size(), back_inserter(lat_row_norm));
    }else{
        lon_col_norm = lon_col;
        lat_row_norm = lat_row;
    };

    if (alt.size()!=lon_col_norm.size()){
        cout<<" alt.size()!=lon_col_norm.size() -> alt vect = alt[0]"<<endl;
        alt_norm.resize(lon_col_norm.size(),alt[0]);
    }else{
        alt_norm =alt;
    };

    return make_tuple(lon_col_norm, lat_row_norm, alt_norm);
}