/*
coding: utf8

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

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

#include "GeoModelTemplate.hpp"
#include "localization.cpp"

// fix 'M_PI': identifier not found (Windows)
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std;
namespace py = pybind11;


double compute_epipolar_angle(array<double,3> const& end_line,
                              array<double,3> const& start_line
                            ){

    // Compute the equation of the epipolar line y = a*x + b and define the epipolar angle
    double alpha = 0.0;

    if(end_line[1]==start_line[1]){

        if(end_line[0]>start_line[0]){
            alpha = 0.5 * M_PI;
        }
        else{
            alpha = -0.5 * M_PI;
        }
    }
    else{

        if(end_line[1]>start_line[1]){
            double slope = (end_line[0]-start_line[0])/\
                (end_line[1]-start_line[1]);
            alpha = atan(slope);
        }
        else{
            double slope = (end_line[0]-start_line[0])/\
                (end_line[1]-start_line[1]);
            alpha = M_PI + atan(slope);
        }

    }
    return alpha;
}

template <typename T>
tuple<double,
        double,
        double,
        double,
        double,
        double> moving_along_axis(GeoModelTemplate const& geom_model_left,
                                                        GeoModelTemplate const&  geom_model_right,
                                                        double                   current_coords_0,
                                                        double                   current_coords_1,
                                                        double                   current_coords_2,
                                                        double                   spacing,
                                                        const T&                 elevation,
                                                        int                      epi_step,
                                                        double                   epi_angles,
                                                        int                      axis
){

    //Check axis
    if (axis != 0 && axis != 1){
        throw runtime_error("axis value is not available");
    }

    epi_angles = epi_angles + (1 - axis) * M_PI / 2. ;

    double unit_vector_along_epi_x = epi_step * spacing * cos(epi_angles);
    double unit_vector_along_epi_y = epi_step * spacing * sin(epi_angles);

    double next_left_0 = current_coords_0 + unit_vector_along_epi_y;
    double next_left_1 = current_coords_1 + unit_vector_along_epi_x;
    double next_left_2 = current_coords_2;

    // Find the corresponding next pixels in the right image
    double next_right_0;
    double next_right_1;
    double next_right_2;

    tie(next_right_0, next_right_1, next_right_2) = coloc(
        geom_model_left, geom_model_right, next_left_0, next_left_1, elevation
    );

    return {next_left_0, next_left_1, next_left_2, next_right_0, next_right_1, next_right_2};
  

}

template <typename T>
tuple<array<double,3>,array<double,3>> compute_local_epipolar_line(
                                                    GeoModelTemplate const& geom_model_left,
                                                    GeoModelTemplate const& geom_model_right,
                                                    double                  left_point_0,
                                                    double                  left_point_1,
                                                    const T&                elevation,
                                                    double                  elevation_offset
){
    // Right correspondent of the left coordinate
    double right_corr_0;
    double right_corr_1;
    double right_corr_2;
    tie(right_corr_0, right_corr_1,right_corr_2) = coloc(
        geom_model_left, geom_model_right, left_point_0, left_point_1, elevation
    );
    double ground_elev = right_corr_2;

    // Find the beginning of the epipolar line in the left image, using right
    // correspondent at lower elevation
    right_corr_2 = ground_elev - elevation_offset;
    double epi_line_start_0;
    double epi_line_start_1;
    double epi_line_start_2;
    tie(epi_line_start_0, epi_line_start_1, epi_line_start_2) = coloc(
        geom_model_right, geom_model_left, right_corr_0, right_corr_1,right_corr_2
    );

    // Find the ending of the epipolar line in the left image, using right
    // correspondent at higher elevation
    right_corr_2 = ground_elev + elevation_offset;
    double epi_line_end_0;
    double epi_line_end_1;
    double epi_line_end_2;
    tie(epi_line_end_0, epi_line_end_1, epi_line_end_2) = coloc(
        geom_model_right, geom_model_left, right_corr_0, right_corr_1,right_corr_2
    );

    array<double,3> epi_line_start = {epi_line_start_0, epi_line_start_1, epi_line_start_2};
    array<double,3> epi_line_end = {epi_line_end_0, epi_line_end_1, epi_line_end_2};
    return {epi_line_start, epi_line_end};
}

template <typename T>
tuple<py::array_t<double>,
py::array_t<double>,
py::array_t<double>,
double> compute_strip_of_epipolar_grid(
    GeoModelTemplate const&  geom_model_left,
    GeoModelTemplate const&  geom_model_right,
    py::array_t<double, py::array::c_style | py::array::forcecast> left_positions_point_in,
    py::array_t<double, py::array::c_style | py::array::forcecast> right_positions_point_in,
    double spacing,
    int axis,
    int strip_size,
    int epi_step = 1,
    const T& elevation = 0.0,
    double elevation_offset = 50.0,
    py::array_t<double, py::array::c_style | py::array::forcecast> epipolar_angles_in = py::none()
){

    // Convert np.array to 1D std::vector<double>
    py::buffer_info buf_left = left_positions_point_in.request();
    double* left_positions_point = static_cast<double*>(buf_left.ptr);

    py::buffer_info buf_right = right_positions_point_in.request();
    double* right_positions_point = static_cast<double*>(buf_right.ptr);

    py::buffer_info buf_epipolar_angles = epipolar_angles_in.request();
    double* ptr = static_cast<double*>(buf_epipolar_angles.ptr);
    vector<double> epipolar_angles;
    epipolar_angles.reserve(buf_epipolar_angles.size);
    epipolar_angles.insert(epipolar_angles.end(), ptr, ptr + buf_epipolar_angles.size);

    //get dim of arrays
    int nb_rows = buf_left.shape[0];
    int nb_cols = buf_left.shape[1];

    // Instantiate mean baseline ratio
    double mean_baseline_ratio = 0.;

    // compute the output size depending on axis
    array<int,3> size_shape;
    if(axis==0){
        size_shape = {strip_size, nb_cols, 3};
    }else{
        size_shape = {nb_rows, strip_size, 3};
    }

    // rename
    double* current_left_point = left_positions_point;
    double* current_right_point = right_positions_point;

    // copy first elements in output grids
    int grid_rows = size_shape[0];
    int grid_cols = size_shape[1];
    int nb_index = grid_rows * grid_cols * 3;

    vector<double> left_grid (nb_index);
    vector<double> right_grid (nb_index);


    // if epipolar angles is not given as input we first compute the epipolar direction, otherwise
    // we consider that baseline ratio has been already computed thus already_computed_ratio
    // index is updated
    double already_computed_ratio;
    array<double,3> local_epi_start_ar;
    array<double,3> local_epi_end_ar;
    vector<double> epipolar_angles_computed(grid_rows*grid_cols);
    vector<double> epipolar_angles_out(grid_rows*grid_cols);
    double local_baseline_ratio;

    int i_row_col;
    int i_row_col_0;
    int i_row_col_1;
    int i_row_col_2;

    int i_grid_cols_0;
    int i_grid_cols_1;
    int i_grid_cols_2;

    int i_strip_0;
    int i_strip_1;
    int i_strip_2;

    //check if we have to compute the first epipolar_angle
    bool compute_first_angle;
    if(isnan(epipolar_angles[0])&&epipolar_angles.size()==1){
        compute_first_angle = true;
        already_computed_ratio = 0.0;
    }else{
        compute_first_angle = false;
        already_computed_ratio = int(epipolar_angles.size());
    }

        
    for(int row=0;row<nb_rows;++row){
        for(int col=0;col<nb_cols;++col){

            // indexes
            i_row_col = row*grid_cols+col;
            i_row_col_0 = row*nb_cols*3 + col*3 + 0;
            i_row_col_1 = row*nb_cols*3 + col*3 + 1;
            i_row_col_2 = row*nb_cols*3 + col*3 + 2;

            i_grid_cols_0 = row*grid_cols*3 + col*3 + 0;
            i_grid_cols_1 = row*grid_cols*3 + col*3 + 1;
            i_grid_cols_2 = row*grid_cols*3 + col*3 + 2;

            i_strip_0 = row*strip_size*3 + col*3 + 0;
            i_strip_1 = row*strip_size*3 + col*3 + 1;
            i_strip_2 = row*strip_size*3 + col*3 + 2;

            //save input
            left_grid[i_grid_cols_0] = current_left_point[i_row_col_0];
            left_grid[i_grid_cols_1] = current_left_point[i_row_col_1];
            left_grid[i_grid_cols_2] = current_left_point[i_row_col_2];
            right_grid[i_grid_cols_0] = current_right_point[i_row_col_0];
            right_grid[i_grid_cols_1] = current_right_point[i_row_col_1];
            right_grid[i_grid_cols_2] = current_right_point[i_row_col_2];


            if(compute_first_angle){
                //compute epipolar angle
                tie(local_epi_start_ar, local_epi_end_ar) = compute_local_epipolar_line(
                    geom_model_left,
                    geom_model_right,
                    current_left_point[i_row_col_0],
                    current_left_point[i_row_col_1],
                    elevation,
                    elevation_offset);

                epipolar_angles_computed[i_row_col] = compute_epipolar_angle(local_epi_end_ar, local_epi_start_ar);

                local_baseline_ratio = sqrt(
                    (local_epi_end_ar[1] - local_epi_start_ar[1]) * (local_epi_end_ar[1] - local_epi_start_ar[1])
                    + (local_epi_end_ar[0] - local_epi_start_ar[0]) * (local_epi_end_ar[0] - local_epi_start_ar[0])
                ) / (2. * elevation_offset);

                mean_baseline_ratio += local_baseline_ratio;
            }else{
                epipolar_angles_computed[i_row_col] = epipolar_angles[row*nb_cols+col];
            }

            //save computed epipolar angle
            epipolar_angles_out[i_row_col] = epipolar_angles_computed[i_row_col];




            // Grid bode computation
            // 1/ move along axis
            // 2/ compute local epipolar line
            // 3/ compute locale epipolar angle
            // 4/ fill the output
            // 5/ compute and increment the baseline ratio

            // axis = 1 : nb_rows=1 and stripe_size=nb_cols
            //            iterator : row=0 and col=point
            // axis = 0 : stripe_size=nb_rows and nb_cols=1
            //            iterator : row=point and col=point
            for (int point=1;point<strip_size;++point){

                tie(current_left_point[i_row_col_0],
                current_left_point[i_row_col_1],
                current_left_point[i_row_col_2],
                current_right_point[i_row_col_0],
                current_right_point[i_row_col_1],
                current_right_point[i_row_col_2]) = moving_along_axis(
                    geom_model_left,
                    geom_model_right,
                    current_left_point[i_row_col_0],
                    current_left_point[i_row_col_1],
                    current_left_point[i_row_col_2],
                    spacing,
                    elevation,
                    epi_step,
                    epipolar_angles_computed[i_row_col], axis
                );

                tie(local_epi_start_ar, local_epi_end_ar) = compute_local_epipolar_line(
                    geom_model_left,
                    geom_model_right,
                    current_left_point[i_row_col_0],
                    current_left_point[i_row_col_1],
                    elevation,
                    elevation_offset
                );

                epipolar_angles_computed[i_row_col] = compute_epipolar_angle(local_epi_end_ar, local_epi_start_ar);

                // Stock values
                if (axis == 0){
                    //size_shape = {strip_size, nb_cols, 3};

                    epipolar_angles_out[i_row_col+point*grid_cols] = epipolar_angles_computed[i_row_col];
                    left_grid[i_row_col_0 + point*nb_cols*3] = current_left_point[i_row_col_0];
                    left_grid[i_row_col_1 + point*nb_cols*3] = current_left_point[i_row_col_1];
                    left_grid[i_row_col_2 + point*nb_cols*3] = current_left_point[i_row_col_2];
                    right_grid[i_row_col_0 + point*nb_cols*3] = current_right_point[i_row_col_0];
                    right_grid[i_row_col_1 + point*nb_cols*3] = current_right_point[i_row_col_1];
                    right_grid[i_row_col_2 + point*nb_cols*3] = current_right_point[i_row_col_2];

                }
                else{
                    //size_shape = {nb_rows, strip_size, 3};

                    epipolar_angles_out[i_row_col+point] = epipolar_angles_computed[row*grid_cols+col];
                    left_grid[i_strip_0+point*3] = current_left_point[i_row_col_0];
                    left_grid[i_strip_1+point*3] = current_left_point[i_row_col_1];
                    left_grid[i_strip_2+point*3] = current_left_point[i_row_col_2];
                    right_grid[i_strip_0+point*3] = current_right_point[i_row_col_0];
                    right_grid[i_strip_1+point*3] = current_right_point[i_row_col_1];
                    right_grid[i_strip_2+point*3] = current_right_point[i_row_col_2];

                }

                local_baseline_ratio = sqrt(
                    (local_epi_end_ar[1] - local_epi_start_ar[1]) * (local_epi_end_ar[1] - local_epi_start_ar[1])
                    + (local_epi_end_ar[0] - local_epi_start_ar[0]) * (local_epi_end_ar[0] - local_epi_start_ar[0])
                ) / (2. * elevation_offset);

                mean_baseline_ratio += local_baseline_ratio;
            }
        }
    }

    // Compute the mean baseline ratio
    mean_baseline_ratio /= grid_rows * grid_cols - already_computed_ratio;

    //to numpy array
    size_t grid_rows_out = (size_t)grid_rows;
    size_t grid_cols_out = (size_t)grid_cols;
    size_t nb_coord = 3;


    auto left_grid_out = py::array_t<double>(
        {grid_rows_out, grid_cols_out,nb_coord},
        left_grid.data()
    );

    auto right_grid_out = py::array_t<double>(
        {grid_rows_out, grid_cols_out,nb_coord},
        right_grid.data()
    );

    auto np_epipolar_angles_out = py::array_t<double>(
        {grid_rows_out, grid_cols_out},
        epipolar_angles_out.data()
    );

    return {left_grid_out, right_grid_out, np_epipolar_angles_out, mean_baseline_ratio};

}
