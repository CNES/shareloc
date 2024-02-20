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

/**
Rectification function
*/

#include <vector>
#include <array>
#include <cmath>
#include <iostream>

#include "GeoModelTemplate.hpp"
#include "localization.cpp"

using namespace std;


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

// vector<double> compute_epipolar_angle(vector<double> const& end_line_row,
//                                       vector<double> const& end_line_col,
//                                       vector<double> const& end_line_alt,
//                                       vector<double> const& start_line_row,
//                                       vector<double> const& start_line_col,
//                                       vector<double> const& start_line_alt
//                                       ){

//     // Compute the equation of the epipolar line y = a*x + b and define the epipolare angle
//     vector<double> alpha(start_line_row.size(),0.0);

//     bool same_col_positive;
//     bool same_col_negative;
//     bool diff_col_pos;
//     bool init_diff_col_pos = false;//has diff_col_pos already been true ?
//     size_t index_diff_col_pos_0;
//     double slope;
//     bool diff_col_neg;
//     bool init_diff_col_neg = false;//has diff_col_neg already been true ?
//     size_t index_diff_col_neg_0;
//     for(size_t i = 0;i<start_line_row.size();++i){

//         same_col_positive = (end_line_col[i]==start_line_col[i])&&\
//                             (end_line_row[i]>start_line_row[i]);
//         if(same_col_positive){
//             alpha[i] = 0.5 * M_PI;
//             }

//         same_col_negative = (end_line_col[i]==start_line_col[i])&&\
//                             (end_line_row[i]<=start_line_row[i]);
//         if(same_col_negative){
//             alpha[i] = -0.5 * M_PI;
//             }

//     diff_col_pos = (end_line_col[i] != start_line_col[i])&&(end_line_col[i]>start_line_col[i]);
//     diff_col_pos = (end_line_col[i] != start_line_col[i])&&(end_line_col[i]>start_line_col[i]);
//         if(diff_col_pos && !init_diff_col_pos){
//             init_diff_col_pos = true;
//             index_diff_col_pos_0 =i;
//         }
//         if(init_diff_col_pos){
//             slope = (end_line_row[index_diff_col_pos_0]-start_line_row[index_diff_col_pos_0])/\
//                     (end_line_col[index_diff_col_pos_0]-start_line_col[index_diff_col_pos_0]);
//             if(diff_col_pos){
//                 alpha[i] = atan(slope);
//                 }
//         }

//    diff_col_neg = (end_line_col[i] != start_line_col[i])&&(end_line_col[i]<=start_line_col[i]);
//    diff_col_neg = (end_line_col[i] != start_line_col[i])&&(end_line_col[i]<=start_line_col[i]);
//         if(diff_col_neg && !init_diff_col_neg){
//             init_diff_col_neg = true;
//             index_diff_col_neg_0 =i;
//         }
//         if(init_diff_col_neg){
//             slope = (end_line_row[index_diff_col_neg_0]-start_line_row[index_diff_col_neg_0])/\
//                     (end_line_col[index_diff_col_neg_0]-start_line_col[index_diff_col_neg_0]);
//             if(diff_col_neg){
//                 alpha[i] = M_PI + atan(slope);
//                 }
//         }

        
//     }

//     return alpha;
// }


template <typename T>
tuple<array<double,3>,array<double,3>> moving_along_axis(GeoModelTemplate const& geom_model_left,
                                                        GeoModelTemplate const&  geom_model_right,
                                                        array<double,3>          current_coords,
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

    double next_left_0 = current_coords[0] + unit_vector_along_epi_y;
    double next_left_1 = current_coords[1] + unit_vector_along_epi_x;
    double next_left_2 = current_coords[2];

    // Find the corresponding next pixels in the right image
    double next_right_0;
    double next_right_1;
    double next_right_2;

    tie(next_right_0, next_right_1, next_right_2) = coloc(
        geom_model_left, geom_model_right, next_left_0, next_left_1, elevation
    );

    array<double,3> next_left = {next_left_0, next_left_1, next_left_2};
    array<double,3> next_right = {next_right_0, next_right_1, next_right_2};
    return {next_left,next_right};

}

template <typename T>
tuple<array<double,3>,array<double,3>> compute_local_epipolar_line(
                                                    GeoModelTemplate const& geom_model_left,
                                                    GeoModelTemplate const& geom_model_right,
                                                    array<double,3>         left_point,
                                                    const T&                elevation,
                                                    double                  elevation_offset
){
    // Right correspondent of the left coordinate
    double right_corr_0;
    double right_corr_1;
    double right_corr_2;
    tie(right_corr_0, right_corr_1,right_corr_2) = coloc(
        geom_model_left, geom_model_right, left_point[0], left_point[1], elevation
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