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
Thes purpose of this module is only to "bind" cpp code.
It gives to the compiler the instructions to compile the usefull cpp code into an .so file
which is callable in a python code as a python module.
*/

#include <vector>
#include <array>
#include <cmath>
#include <iostream>

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
