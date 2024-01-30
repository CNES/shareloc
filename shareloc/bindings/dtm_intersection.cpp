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

using namespace std;
namespace py = pybind11;


//---- DTMIntersection methodes ----//

DTMIntersection::DTMIntersection(
        int dtm_image_epsg,
        py::array_t<double, py::array::c_style | py::array::forcecast> dtm_image_alt_data,
        int dtm_image_nb_rows,
        int dtm_image_nb_columns,
        tuple<double,double,double,double,double,double> dtm_image_transform
    ){

    epsg = dtm_image_epsg;
    tol_z = 0.0001;

    nb_rows = dtm_image_nb_rows;
    nb_columns = dtm_image_nb_columns;

    py::buffer_info buf_info = dtm_image_alt_data.request();
    double* ptr = static_cast<double*>(buf_info.ptr);
    alt_data.insert(alt_data.end(), ptr, ptr + buf_info.size);

    tie(alt_min_cell,alt_max_cell) = init_min_max(alt_data,
                                                    dtm_image_nb_rows,
                                                    dtm_image_nb_columns);

    alt_min = *min_element(alt_data.begin(), alt_data.end());
    alt_max = *max_element(alt_data.begin(), alt_data.end());

    plane_coef_a = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    plane_coef_b = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    plane_coef_c = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0};
    plane_coef_d = {0.0,
                    dtm_image_nb_rows - 1.0,
                    0.0,
                    dtm_image_nb_columns - 1.0,
                    alt_min,alt_max};


    plans = {1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, dtm_image_nb_rows - 1.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, dtm_image_nb_columns - 1.0,
    0.0, 0.0, 1.0, alt_min,
    0.0, 0.0, 1.0, alt_max};//2D -> 1D nb_columns = 4


    apply([&](auto... args) { transform = {args...}; }, dtm_image_transform);
    

    //det = a*e-b*d
    double idet = 1.0/(transform[1]*transform[5]-transform[2]*transform[4]);

    double ra = transform[5] * idet;
    double rb = -transform[2] * idet;
    double rd = -transform[4] * idet;
    double re = transform[1] * idet;

    // invert transform=[c,a,b,f,d,e]
    trans_inv[1] = ra;
    trans_inv[2] = rb;
    trans_inv[0] = -transform[0] * ra - transform[3] * rb;
    trans_inv[4] = rd;
    trans_inv[5] = re;
    trans_inv[3] = -transform[0] * rd - transform[3] * re;

}


double DTMIntersection::eq_plan(int i, array<double, 3> const& position)const{

return plane_coef_a[i] * position[0]
    + plane_coef_b[i] * position[1]
    + plane_coef_c[i] * position[2]
    - plane_coef_d[i];

}

array<double, 3> DTMIntersection::ter_to_index(array<double, 3> const& vect_ter)const{

    double vx= vect_ter[0];//row
    double vy= vect_ter[1];//col

    array<double, 3> ter;

    ter[1] = (vx * trans_inv[1] + vy * trans_inv[2] + trans_inv[0])-0.5;
    ter[0] = (vx * trans_inv[4] + vy * trans_inv[5] + trans_inv[3])-0.5;
    ter[2] = vect_ter[2];

    return ter;
}

vector<double> DTMIntersection::ter_to_indexs(vector<double> const& vect_ter){//maybe unecessary
    vector<double> res;
    return res;
}

array<double, 3> DTMIntersection::index_to_ter(array<double, 3> const& vect_ter)const{

    
    double vx= vect_ter[1]+0.5;//col
    double vy= vect_ter[0]+0.5;//row

    array<double, 3> ter;

    ter[0] = vx * transform[1] + vy * transform[2] + transform[0];
    ter[1] =  vx * transform[4] + vy * transform[5] + transform[3];
    ter[2] = vect_ter[2];

    return ter;
}


double DTMIntersection::interpolate(double delta_shift_row, double delta_shift_col)const{

    //-- Initialise rows
    double lower_shift_row;
    if (delta_shift_row < 0.0){
        lower_shift_row = 0.0;
        }
    else if (delta_shift_row >= nb_rows - 1.0){
        lower_shift_row = nb_rows - 2.0;
        }     
    else{
        lower_shift_row = static_cast<int>(floor(delta_shift_row));
    }
    double upper_shift_row = lower_shift_row + 1.0;

    //- Initialise cols
    double lower_shift_col;
    if (delta_shift_col < 0.0){
        lower_shift_col = 0.0;
        }
    else if (delta_shift_col >= nb_columns - 1.0){
        lower_shift_col = nb_columns - 2.0;
        }     
    else{
        lower_shift_col = int(floor(delta_shift_col));
    }
    double upper_shift_col = lower_shift_col + 1.0;

    // (col_shift, row_shift) are subpixel distance to interpolate along each axis
    double col_shift = delta_shift_col - lower_shift_col;
    double row_shift = delta_shift_row - lower_shift_row;

    // Altitude
    double mati = 
        (1 - col_shift)*(1 - row_shift) * alt_data[lower_shift_row * nb_columns + lower_shift_col]
        + col_shift * (1 - row_shift) * alt_data[lower_shift_row * nb_columns + upper_shift_col]
        + (1 - col_shift) * row_shift * alt_data[upper_shift_row * nb_columns + lower_shift_col]
        + col_shift * row_shift * alt_data[upper_shift_row * nb_columns + upper_shift_col];
    return mati;
}

tuple<bool,
bool,
vector<double>,
bool,
vector<double>> DTMIntersection::intersect_dtm_cube(vector<double> const& los) const{
    tuple<bool,bool,vector<double>,bool,vector<double>> res;
    return res;
}

tuple<bool,
bool,
vector<double>> DTMIntersection::intersection(
    vector<double> const& los_index,
    vector<double> const& point_b,
    double h_intersect) const
{
    tuple<bool,bool,vector<double>> res;
    return res;
}



//-- function --//

tuple<vector<double>,
vector<double>> init_min_max(vector<double> const& alt_data,int nb_rows,int nb_columns)
{

vector<double> alt_min_cell ((nb_rows-1)*(nb_columns-1));
vector<double> alt_max_cell ((nb_rows-1)*(nb_columns-1));

for(int i = 0; i < nb_rows-1; ++i){
    for(int j = 0; j < nb_columns-1; ++j){

        double subarray1_i = alt_data[nb_columns * i + j];
        double subarray2_i = alt_data[nb_columns * (i+1) + j];
        double subarray3_i = alt_data[nb_columns * i + (j+1)];
        double subarray4_i = alt_data[nb_columns * (i+1) + (j+1)];

        alt_min_cell[(nb_columns-1)*i+j] = floor(min({subarray1_i,
                                                        subarray2_i,
                                                        subarray3_i,
                                                        subarray4_i}));

        alt_max_cell[(nb_columns-1)*i+j] = ceil(max({subarray1_i,
                                                        subarray2_i,
                                                        subarray3_i,
                                                        subarray4_i}));
        
    }
}

return make_tuple(alt_min_cell,alt_max_cell);
}
