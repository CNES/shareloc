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

RPC::RPC(array<double, 20> num_col_input,
        array<double, 20> den_col_input,
        array<double, 20> num_row_input,
        array<double, 20> den_row_input,
        array<double, 10> norm_coeffs):GeoModelTemplate(){

    lim_extrapol = 1.0001;

    std::copy(num_col_input.begin(), num_col_input.end(), num_col.begin());
    std::copy(den_col_input.begin(), den_col_input.end(), den_col.begin());
    std::copy(num_row_input.begin(), num_row_input.end(), num_row.begin());
    std::copy(den_row_input.begin(), den_row_input.end(), den_row.begin());

    offset_lon = norm_coeffs[0];//offset_x
    scale_lon = norm_coeffs[1];
    offset_lat = norm_coeffs[2];//offset_t
    scale_lat = norm_coeffs[3];
    offset_alt = norm_coeffs[4];
    scale_alt = norm_coeffs[5];
    offset_col = norm_coeffs[6];
    scale_col = norm_coeffs[7];
    offset_row = norm_coeffs[8];
    scale_row = norm_coeffs[9];

}

vector<vector<double>> RPC::direct_loc_h(
    vector<double> row,
    vector<double> col,
    double alt,
    bool fill_nan){
    vector<vector<double>> vect;
    return vect;
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
    
    {
    // vector<double> num_col_default = {0.}; // CHANGE rpc.hpp
    // if (num_col!=num_col_default){ //so self.inverse_coefficient = True
    //     if (!is_same<decltype(lon), vector<double>>::value){//if lon not a vector<double>
    //         try{
    //             lon = vector<double>{lon};
    //             lat = vector<double>{lat};
    //         }catch(...){
    //             cerr<<"Error : longitude/latitude not vector<double> and can't be converted"<<endl;
    //         }
    //         };
    //     if (!is_same<decltype(alt), vector<double>>::value){//if alt not a vector<double>
    //         try{
    //             alt = vector<double>{alt};
    //         }catch(...){
    //             cerr<<"Error : altitude not vector<double> and can't be converted"<<endl;
    //         }
    //         };

        
    //     vector<double> lon_norm;
    //     vector<double> lat_norm;
    //     vector<double> alt_norm;
        
    //     // Check longitudes and latitudes sizes -> lon_norm and lat_norm
    //     if(lat.size()<lon.size()){
    //         cout<<"Warning : Inverse loc : compute_rational_function_polynomial : lat.size()!=lon.size() -> truncate lon"<<endl;
    //         copy(lon.begin(), lon.begin()+lat.size(), back_inserter(lon_norm));
    //         copy(lat.begin(), lat.begin(), back_inserter(lat_norm));
    //     }else if (lat.size()>lon.size()){
    //         cout<<"Warning : Inverse loc : compute_rational_function_polynomial : lat.size()!=lon.size() -> truncate lat "<<endl;
    //         copy(lon.begin(), lon.begin(), back_inserter(lon_norm));
    //         copy(lat.begin(), lat.begin()+lon.size(), back_inserter(lat_norm));
    //     }else{
    //         copy(lon.begin(), lon.begin(), back_inserter(lon_norm));
    //         copy(lat.begin(), lat.begin(), back_inserter(lat_norm));
    //     };

    //     //check altitude size -> lat_norm
    //     if (alt.size()!=lon_norm.size()){
    //         cout<<"Warning : Inverse loc : compute_rational_function_polynomial : alt.size()!=lon_norm.size() -> alt vect = alt[0]"<<endl;
    //         alt_norm.resize(lon_norm.size(),alt[0]);
    //     }else{
    //         copy(alt.begin(), alt.begin(), back_inserter(alt_norm));
    //     }

    //     vector<double> alt_res = alt_norm;//save not normalised alt

    //     //Normalisation
    //     for (int i =0;i<int(lon_norm.size());++i){
    //         lon_norm[i] = (lon_norm[i] - offset_lon)/scale_lon;
    //         lat_norm[i] = (lat_norm[i] - offset_lat)/scale_lat;
    //         alt_norm[i] = (alt_norm[i] - offset_alt)/scale_alt;

    //         if(abs(lon_norm[i])>lim_extrapol || abs(lat_norm[i])>lim_extrapol || abs(alt_norm[i])>lim_extrapol){
    //             cerr<<"Error : normalisation values exceed lim_extrapol"<<endl;
    //         }
    //     }

    //     vector<double> row_out, col_out;
    //     tie(row_out, col_out) = compute_rational_function_polynomial(
    //         lon_norm,
    //         lat_norm,
    //         alt_norm,
    //         num_col,
    //         den_col,
    //         num_row,
    //         den_row,
    //         scale_col,
    //         offset_col,
    //         scale_row,
    //         offset_row
    //     );
    //     return make_tuple(row_out, col_out, alt_res);
    // }else{
    //     cout<<"inverse localisation can't be performed, inverse coefficients have not been defined"<<endl;
    //     tuple<vector<double>, vector<double>, vector<double>> res;
    //     return res;
    //  };
    tuple<vector<double>,vector<double>,vector<double>> res;
    return res;
}


vector<vector<double>> RPC::filter_coordinates(
    vector<double> first_coord,
    vector<double> second_coord,
    bool fill_nan,
    string direction){
    vector<vector<double>> vect;
    return vect;
}

tuple<vector<double>,
vector<double>,
vector<double>,
vector<double>> RPC::compute_loc_inverse_derivates(
    vector<double> lon,
    vector<double> lat,
    vector<double> alt){
    tuple<vector<double>,vector<double>,vector<double>,vector<double>> res;
    return res;
}

vector<vector<double>> RPC::direct_loc_inverse_iterative(
    vector<double> row,
    vector<double> col,
    double alt,
    int nb_iter_max,
    bool fill_nan){
    vector<vector<double>> vect;
    return vect;
}

vector<double> RPC::get_alt_min_max(){
    vector<double> vect;
    return vect;
}

vector<vector<double>> RPC::los_extrema(
    double row,
    double col,
    double alt_min,
    double alt_max,
    bool fill_nan,
    int epsg){
    vector<vector<double>> vect;
    return vect;
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


//---- Functions ----//

double polynomial_equation(
    double xnorm,
    double ynorm,
    double znorm,
    array<double, 20> coeff){//const array<double, 20>& coeff
    
    return
    coeff[0]
    + coeff[1] * xnorm
    + coeff[2] * ynorm
    + coeff[3] * znorm
    + coeff[4] * xnorm * ynorm
    + coeff[5] * xnorm * znorm
    + coeff[6] * ynorm * znorm
    + coeff[7] * xnorm*xnorm
    + coeff[8] * ynorm*ynorm
    + coeff[9] * znorm*znorm
    + coeff[10] * xnorm * ynorm * znorm
    + coeff[11] * xnorm*xnorm*xnorm
    + coeff[12] * xnorm * ynorm*ynorm
    + coeff[13] * xnorm * znorm*znorm
    + coeff[14] * xnorm*xnorm * ynorm
    + coeff[15] * ynorm*ynorm*ynorm
    + coeff[16] * ynorm * znorm*znorm
    + coeff[17] * xnorm*xnorm * znorm
    + coeff[18] * ynorm*ynorm * znorm
    + coeff[19] * znorm*znorm*znorm;
}



tuple<vector<double>,vector<double>> compute_rational_function_polynomial(
    vector<double> lon_col_norm,
    vector<double> lat_row_norm,
    vector<double> alt_norm,
    array<double, 20> num_col,
    array<double, 20> den_col,
    array<double, 20> num_lin,
    array<double, 20> den_lin,
    double scale_col,
    double offset_col,
    double scale_lin,
    double offset_lin
){
    if (lon_col_norm.size() != alt_norm.size()){
        cerr<<"Error : not same number of altitudes and longitudes"<<endl;
        exit(EXIT_FAILURE);
    }

    vector<double> col_lat_out(lon_col_norm.size());
    vector<double> row_lon_out(lon_col_norm.size());

    double poly_num_col;
    double poly_den_col;
    double poly_num_lin;
    double poly_den_lin;
    for(int i = 0;i<(int)lon_col_norm.size();++i){
        poly_num_col = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], num_col);
        poly_den_col = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], den_col);
        poly_num_lin = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], num_lin);
        poly_den_lin = polynomial_equation(lon_col_norm[i], lat_row_norm[i], alt_norm[i], den_lin);
        col_lat_out[i] = poly_num_col / poly_den_col * scale_col + offset_col;
        row_lon_out[i] = poly_num_lin / poly_den_lin * scale_lin + offset_lin;
    }

    tuple<vector<double>,vector<double>> res = make_tuple(row_lon_out, col_lat_out);
    return res;
}



double derivative_polynomial_latitude(
    double lon_norm,
    double lat_norm,
    double alt_norm,
    vector<double> coeff){
    double res;
    return res;
}


double derivative_polynomial_longitude(
    double lon_norm,
    double lat_norm,
    double alt_norm,
    vector<double> coeff){
    double res;
    return res;
}


tuple<vector<double>,
vector<double>,
vector<double>,
vector<double>> compute_loc_inverse_derivates_optimized(
    vector<double> lon_norm,
    vector<double> lat_norm,
    vector<double> alt_norm,
    vector<double> num_col,
    vector<double> den_col,
    vector<double> num_lin,
    vector<double> den_lin,
    double scale_col,
    double scale_lon,
    double scale_lin,
    double scale_lat
){
    tuple<vector<double>,vector<double>,vector<double>,vector<double>> res;
    return res;
}