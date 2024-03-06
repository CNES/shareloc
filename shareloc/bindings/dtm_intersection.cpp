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

#include "dtm_intersection.hpp"
#include <iostream>

using namespace std;
namespace py = pybind11;


//---- DTMIntersection methods ----//

DTMIntersection::DTMIntersection(){}

DTMIntersection::DTMIntersection(
        int dtm_image_epsg,
        py::array_t<double, py::array::c_style | py::array::forcecast> dtm_image_alt_data,
        int dtm_image_nb_rows,
        int dtm_image_nb_columns,
        tuple<double,double,double,double,double,double> dtm_image_transform
    ){

    m_epsg = dtm_image_epsg;
    m_tol_z = 0.0001;

    m_nb_rows = dtm_image_nb_rows;
    m_nb_columns = dtm_image_nb_columns;

    py::buffer_info buf_info = dtm_image_alt_data.request();
    double* ptr = static_cast<double*>(buf_info.ptr);
    m_alt_data.insert(m_alt_data.end(), ptr, ptr + buf_info.size);

    tie(m_alt_min_cell,m_alt_max_cell) = init_min_max(m_alt_data,
                                                    dtm_image_nb_rows,
                                                    dtm_image_nb_columns);

    m_alt_min = *min_element(m_alt_data.begin(), m_alt_data.end());
    m_alt_max = *max_element(m_alt_data.begin(), m_alt_data.end());

    m_plane_coef_a = {1.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    m_plane_coef_b = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    m_plane_coef_c = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0};
    m_plane_coef_d = {0.0,
                    dtm_image_nb_rows - 1.0,
                    0.0,
                    dtm_image_nb_columns - 1.0,
                    m_alt_min,m_alt_max};


    m_plans = {1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, dtm_image_nb_rows - 1.0,
    0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, dtm_image_nb_columns - 1.0,
    0.0, 0.0, 1.0, m_alt_min,
    0.0, 0.0, 1.0, m_alt_max};//2D -> 1D : nb_columns = 4


    apply([&](auto... args) { m_transform = {args...}; }, dtm_image_transform);
    

    //det = a*e-b*d
    double idet = 1.0/(m_transform[1]*m_transform[5]-m_transform[2]*m_transform[4]);

    double ra = m_transform[5] * idet;
    double rb = -m_transform[2] * idet;
    double rd = -m_transform[4] * idet;
    double re = m_transform[1] * idet;

    // invert transform=[c,a,b,f,d,e]
    m_trans_inv[1] = ra;
    m_trans_inv[2] = rb;
    m_trans_inv[0] = -m_transform[0] * ra - m_transform[3] * rb;
    m_trans_inv[4] = rd;
    m_trans_inv[5] = re;
    m_trans_inv[3] = -m_transform[0] * rd - m_transform[3] * re;

}


double DTMIntersection::eq_plan(int i, array<double, 3> const& position)const{

return m_plane_coef_a[i] * position[0]
    + m_plane_coef_b[i] * position[1]
    + m_plane_coef_c[i] * position[2]
    - m_plane_coef_d[i];

}

array<double, 3> DTMIntersection::ter_to_index(array<double, 3> const& vect_ter)const{

    double vx= vect_ter[0];//row
    double vy= vect_ter[1];//col

    array<double, 3> ter;

    ter[1] = (vx * m_trans_inv[1] + vy * m_trans_inv[2] + m_trans_inv[0])-0.5;
    ter[0] = (vx * m_trans_inv[4] + vy * m_trans_inv[5] + m_trans_inv[3])-0.5;
    ter[2] = vect_ter[2];

    return ter;
}

array<double, 3> DTMIntersection::index_to_ter(array<double, 3> const& vect_ter)const{

    
    double vx= vect_ter[1]+0.5;//col
    double vy= vect_ter[0]+0.5;//row

    array<double, 3> ter;

    ter[0] = vx * m_transform[1] + vy * m_transform[2] + m_transform[0];
    ter[1] =  vx * m_transform[4] + vy * m_transform[5] + m_transform[3];
    ter[2] = vect_ter[2];

    return ter;
}

py::array_t<double> DTMIntersection::get_footprint_corners()const{

    vector<double> res = {-0.5,-0.5,\
                        -0.5,this->get_nb_columns()-0.5,\
                        this->get_nb_rows()-0.5,this->get_nb_columns()-0.5,\
                        this->get_nb_rows()-0.5,-0.5};
    
    return py::array_t<double>({4, 2}, res.data());
}

double DTMIntersection::interpolate(double delta_shift_row, double delta_shift_col)const{

    //-- Initialise rows
    double lower_shift_row;
    if (delta_shift_row < 0.0){
        lower_shift_row = 0.0;
        }
    else if (delta_shift_row >= m_nb_rows - 1.0){
        lower_shift_row = m_nb_rows - 2.0;
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
    else if (delta_shift_col >= m_nb_columns - 1.0){
        lower_shift_col = m_nb_columns - 2.0;
        }     
    else{
        lower_shift_col = int(floor(delta_shift_col));
    }
    double upper_shift_col = lower_shift_col + 1.0;

    // (col_shift, row_shift) are subpixel distance to interpolate along each axis
    double col_shift = delta_shift_col - lower_shift_col;
    double row_shift = delta_shift_row - lower_shift_row;

    // Altitude
    double interp_value = 
        (1-col_shift)*(1-row_shift) * m_alt_data[lower_shift_row * m_nb_columns + lower_shift_col]
        + col_shift*(1 - row_shift) * m_alt_data[lower_shift_row * m_nb_columns + upper_shift_col]
        + (1 - col_shift)*row_shift * m_alt_data[upper_shift_row * m_nb_columns + lower_shift_col]
        + col_shift * row_shift * m_alt_data[upper_shift_row * m_nb_columns + upper_shift_col];
    return interp_value;
}




tuple<bool, 
array<double,3>,
double,
vector<double>,
vector<double>,
vector<double>> DTMIntersection::intersect_dtm_cube(vector<double> const& los_x,
                                                    vector<double> const& los_y,
                                                    vector<double> const& los_z)const{

    array<double,3> point_b;
    double h_intersect;

    vector<double> coord_col_i;
    vector<double> coord_row_i;
    vector<double> coord_alt_i;
    vector<double> alti_layer_i;

    size_t nbalt = los_x.size();

    vector<double> los_x_index(nbalt);
    vector<double> los_y_index(nbalt);
    vector<double> los_z_index(nbalt);

    for (size_t i =0;i<nbalt;++i ){//method ter_to_indexS
        array<double,3> los_index = ter_to_index({los_x[i],los_y[i],los_z[i]}); //to optimize arg
        los_x_index[i] = los_index[0];
        los_y_index[i] = los_index[1];
        los_z_index[i] = los_index[2];

    }

    // -----------------------------------------------------------------------
    // Number of valid intersections found
    int nbi = 0;

    double los_x_hat;
    double los_y_hat;
    double los_z_hat;

    double los_a;
    double los_hat_onplane;

    array<double,3> s_a;

    // -----------------------------------------------------------------------
    // We loop on the 6 planes of the DTM cube
    for(int plane_index=0;plane_index<6;++plane_index){

        // -----------------------------------------------------------------------
        // Init the vertex of the geometric line of sight
        los_x_hat = los_x_index[0];
        los_y_hat = los_y_index[0];
        los_z_hat = los_z_index[0];

        // -----------------------------------------------------------------------
        // Init position parallel to the plane
        // print self.plans[plane_index,:-1]
        // los_hat_onplane = (self.plans[plane_index,:-1]*los_hat).sum() - self.d[plane_index]
        // print los_hat_onplane
        los_hat_onplane = eq_plan(plane_index, {los_x_hat,los_y_hat,los_z_hat});
        // -----------------------------------------------------------------------
        // Loop on line of sight segments
        // and control if we cross or not the current face plane_index of DTM cube

        for(size_t alti_layer=0;alti_layer<nbalt;++alti_layer){
            // -----------------------------------------------------------------------
            // Transfer point B into point A
            los_a = los_hat_onplane;
            s_a = {los_x_hat,los_y_hat,los_z_hat};

            // -----------------------------------------------------------------------
            // Reinit point B
            los_x_hat = los_x_index[alti_layer];// we iterate on different los points
            los_y_hat = los_y_index[alti_layer];
            los_z_hat = los_z_index[alti_layer];

            // -----------------------------------------------------------------------
            // Init position parallel to the plane
            los_hat_onplane = eq_plan(plane_index,{los_x_hat,los_y_hat,los_z_hat});

            // -----------------------------------------------------------------------
            // Intersection test: los_a and los_hat_onplane with opposite sign
            if (los_a * los_hat_onplane <= 0){
                if(los_a == 0 && los_hat_onplane == 0){
                    throw runtime_error("C++ : too many solutions in DTM intersection");
                }
                // -----------------------------------------------------------------------
                if(los_a == 0 ){
                    // A is solution (it is on the same plane)
                    coord_col_i.push_back(s_a[0]);
                    coord_row_i.push_back(s_a[1]);
                    coord_alt_i.push_back(s_a[2]);
                    alti_layer_i.push_back(alti_layer - 1.0);
                }
                // -----------------------------------------------------------------------
                else if(los_hat_onplane == 0.0 ){
                    // A is solution (it is on the same plane)
                    coord_col_i.push_back(los_x_hat);
                    coord_row_i.push_back(los_y_hat);
                    coord_alt_i.push_back(los_z_hat);
                    alti_layer_i.push_back(alti_layer);
                }
                // -----------------------------------------------------------------------
                else{
                    // -----------------------------------------------------------------------
                    // A and B are on either side of the plane
                    // Intersection interpolation coefficients between A and B
                    double interp_coef_a = los_hat_onplane / (los_hat_onplane - los_a);
                    double interp_coef_b = -los_a / (los_hat_onplane - los_a);
                    // Assignment or interpolation
                    // NB : to avoid test problems
                    //      <BeOnCube> (see further)
                    // . coordinate <u> (line)
                    // -----------------------------------------------------------------------
                    if(plane_index < 2.0){
                        coord_col_i.push_back(m_plane_coef_d[plane_index]);
                        }
                    // -----------------------------------------------------------------------
                    else{
                        coord_col_i.push_back(interp_coef_a * s_a[0] + interp_coef_b * los_x_hat);
                        }
                    // -----------------------------------------------------------------------
                    // . coordinate <v> (column)
                    if(1.0<plane_index && plane_index<4.0){
                        coord_row_i.push_back(m_plane_coef_d[plane_index]);
                        }
                    else{
                        coord_row_i.push_back(interp_coef_a * s_a[1] + interp_coef_b * los_y_hat);
                        }
                    // -----------------------------------------------------------------------
                    // . coordinate <z> (altitude)
                    if(plane_index > 3.0){
                        coord_alt_i.push_back(m_plane_coef_d[plane_index]);
                        }
                    // -----------------------------------------------------------------------
                    else{
                        coord_alt_i.push_back(interp_coef_a * s_a[2] + interp_coef_b * los_z_hat);
                        }

                    alti_layer_i.push_back(alti_layer - interp_coef_a);
                }
                // -----------------------------------------------------------------------
                // Incrementing the number of intersections found
                nbi ++;
                // -----------------------------------------------------------------------
                // Switch to the next face of the cube
                break;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Sorting points along line of sight (there are at least two)
    // Arrange them in ascending and descending order
    double dtmp;
    for(int alti_layer = 0;alti_layer<nbi;++alti_layer){
        for (int next_alti_layer = alti_layer + 1;next_alti_layer<nbi;++next_alti_layer){
            if (alti_layer_i[next_alti_layer] < alti_layer_i[alti_layer]){

                dtmp = coord_col_i[alti_layer];
                coord_col_i[alti_layer] = coord_col_i[next_alti_layer];
                coord_col_i[next_alti_layer] = dtmp;

                dtmp = coord_row_i[alti_layer];
                coord_row_i[alti_layer] = coord_row_i[next_alti_layer];
                coord_row_i[next_alti_layer] = dtmp;

                dtmp = coord_alt_i[alti_layer];
                coord_alt_i[alti_layer] = coord_alt_i[next_alti_layer];
                coord_alt_i[next_alti_layer] = dtmp;

                dtmp = alti_layer_i[alti_layer];
                alti_layer_i[alti_layer] = alti_layer_i[next_alti_layer];
                alti_layer_i[next_alti_layer] = dtmp;

            }
        }
    }
    // -----------------------------------------------------------------------
    // Filtering points not located on the cube
    int alti_layer = 0;
    while(alti_layer < nbi){
        // test inside the cube
        bool test_on_cube = (
            (coord_col_i[alti_layer] >= m_plane_coef_d[0]) &&
            (coord_col_i[alti_layer] <= m_plane_coef_d[1]) &&
            (coord_row_i[alti_layer] >= m_plane_coef_d[2]) &&
            (coord_row_i[alti_layer] <= m_plane_coef_d[3]) &&
            (coord_alt_i[alti_layer] >= m_plane_coef_d[4]) &&
            (coord_alt_i[alti_layer] <= m_plane_coef_d[5])
        );

        if (!test_on_cube){
            // We translate all the following points (we overwrite this invalid point)
            for (int next_alti_layer = alti_layer + 1;next_alti_layer<nbi;++next_alti_layer){

                coord_col_i[next_alti_layer - 1] = coord_col_i[next_alti_layer];
                coord_row_i[next_alti_layer - 1] = coord_row_i[next_alti_layer];
                coord_alt_i[next_alti_layer - 1] = coord_alt_i[next_alti_layer];
                alti_layer_i[next_alti_layer - 1] = alti_layer_i[next_alti_layer];    

            }
            --nbi;
        }
        else{
            ++alti_layer;
        }
    }

    // -----------------------------------------------------------------------
    // No solution if 0 or 1 single point is found (we have tangent to the cube)
    if(nbi<2){
        return make_tuple(false, point_b, h_intersect, los_x_index, los_y_index, los_z_index);
    }

    // -----------------------------------------------------------------------
    // There are only 2 points left so we cross the cube
    // LAIG-FA-MAJA-2168-CNES: no more filtering on identical points. 
    // There may be a number of points > 2
    // Init the current point
    // DTM coordinates
    point_b[0] = coord_col_i[0];
    point_b[1] = coord_row_i[0];
    point_b[2] = coord_alt_i[0];
    // point_dtm is the first intersection with the cube (line, column)
    // -----------------------------------------------------------------------
    // h is gld 3D
    h_intersect = alti_layer_i[0];
    // h_intersect is the h interpolation index (not integer)
    // -----------------------------------------------------------------------
    // End, return

    return make_tuple(true, point_b, h_intersect, los_x_index, los_y_index, los_z_index);
}

tuple<bool,double,double,double> DTMIntersection::intersection(
    vector<double> const& los_x_index,//rpc -> size=2
    vector<double> const& los_y_index,
    vector<double> const& los_z_index,
    array<double, 3> const& point_b,
    double h_intersect) const
{
    size_t npl = los_x_index.size();
    array<double,3> point_r;

    vector<double> alti(npl+1);
    for(double index_alti=npl;index_alti>-1.0;index_alti-=1.0){alti[index_alti]=index_alti;}

    array<double,3> p_1 = point_b;
    double h_intersect_p1 = h_intersect;
    int n_row = m_nb_rows;
    int n_col = m_nb_columns;

    // 1 - Init and preliminary tests
    //   1.1 - Test if the vertex is above the DTM
    //       - Compute DTM altitude ? vertex position
    double alti_1 = interpolate(p_1[0], p_1[1]);
    //       - Compute the altitude difference to the DTM
    double d_alti_1 = p_1[2] - alti_1;

    //       - Test if the new top point is above the DTM
    if (d_alti_1<0.0){
            //       - The Point is below the DTM
            //          . means that the line of sight goes into the DTM by the side
            //          . then below, no solution.
        double point_r_x = point_r[0];
        double point_r_y = point_r[1];
        double point_r_z = point_r[2];
        return make_tuple(false,point_r_x,point_r_y,point_r_z);
    };

    //   1.2 - Init the rank of the first vertex of the line of sight
    size_t i_0 = size_t(floor(h_intersect_p1)); //size_t ?

    //   1.3 - Init the starting point (in p_2)
    array<double,3> p_2 = point_b;
    double h_intersect_p2 = h_intersect;

    // 2. - Loop on the grid planes
    size_t nb_planes = los_x_index.size();// TODO:use npl, in python as well
    while(i_0<(nb_planes-1)){// /! int vs size_t
        // 2.1 - Init current vertex of line of sight
        double col_0 = los_x_index[i_0];
        double row_0 = los_y_index[i_0];
        double z_0 = los_z_index[i_0];
        double z_1 = los_z_index[i_0+1];

        // 2.2 - Init line of sight DTM
        double los_dtm_x = los_x_index[i_0+1] - los_x_index[i_0];
        double los_dtm_y = los_y_index[i_0+1] - los_y_index[i_0];
        double los_dtm_z = los_z_index[i_0+1] - los_z_index[i_0];

        // 2.3 - Test if line of sight is vertical
        if(los_dtm_x == 0 && los_dtm_y==0){

            // 2.3.1 - LOS is  vertical:
            //    - Compute DTM altitude ? vertex position
            alti_1 = interpolate(col_0,row_0);

            //    Test if the next plane is above DTM
            if(los_z_index[i_0+1]<=alti_1){
                // Init exit point
                p_1[0] = col_0;
                p_1[1] = row_0;
                p_1[2] = alti_1;
                point_r = index_to_ter(p_1);

                double point_r_x = point_r[0];
                double point_r_y = point_r[1];
                double point_r_z = point_r[2];

                return make_tuple(true,point_r_x,point_r_y,point_r_z);
            }
            // Positioning on next vertex
            i_0 += 1;
        }
        else{
            // 2.3.2 - LOS is not vertical :
            //         it will fly over the DTM
            //         . we can continue
            //         . remains to demonstrate that the LOS will crash on DTM...
            //
            // Init starting point
            // Init its LOS x-axis
            double a_2 = h_intersect_p2 - i_0;

            // Fixed an FA DG 10 bug, when  resetting, a_2 value can be 1
            // Which has the consequence of jumping to the next segment of the LOS
            // Thus, one can lose points on a slice of altitude
            // To be sure to scan the aiming segment, we reset
            // the starting point of the segment, i.e. a_2 = 0
            if(a_2>=1.0){a_2 = 0.0;};

            // Init first intersected mesh
            //  - Init mesh index
            int col_c = int(floor(p_2[0]));
            int row_c = int(floor(p_2[1]));

            // NB :    caution before to start:
            //        . we put ourselves on the right side of the mesh
            //        . in principle, you should not leave the DTM
            // We enter from the bottom, the DTM mesh is the previous one
            if(p_2[0] == col_c && los_dtm_x<0.0){col_c=1;}
            // We enter from the left, the DTM mesh is the previous one
            if(p_2[1] == row_c && los_dtm_y<0.0){row_c=1;}

            // LDD - We're already out of bounds, we stop
            if(!(a_2<1 && -1<col_c && col_c<(n_row - 1) && -1<row_c && row_c<(n_col - 1))){

                double point_r_x = point_r[0];
                double point_r_y = point_r[1];
                double point_r_z = point_r[2];

                return make_tuple(false,point_r_x,point_r_y,point_r_z);
            }

            // Iterative search loop of the intersected cell
            while(a_2<1 && -1<col_c && col_c<(n_row - 1) && -1<row_c && row_c<(n_col - 1)){
                // - Min and max altitudes of the mesh
                double h_i = m_alt_min_cell[col_c*(m_nb_columns-1)+row_c];
                double h_s = m_alt_max_cell[col_c*(m_nb_columns-1)+row_c];

                // - Transfer: the low point becomes the high point
                // a1 = a_2;
                // p_1 becomes p_2
                p_1 = p_2;
                h_intersect_p1 = h_intersect_p2;

                // 4.2 - Determination of a new low point
                //      - LOS orientation test
                if(los_dtm_x==0.0){
                    // 4.2.1 - LOS is completely oriented  east-west
                    //   p_2[0] = p_1[0] ; // useless, is already init
                    //       - LOS orientation test
                    if(los_dtm_y<0.0){

                        // 4.2.1.1 - LOS goes due west
                        p_2[1] = row_c;
                        row_c -= 1.0;

                    }
                    else{

                        // 4.2.1.2 - LOS goes due east
                        row_c += 1.0;
                        p_2[1] = row_c;

                    }

                    a_2 = (p_2[1] - row_0)/los_dtm_y;
                    p_2[2] = z_0 + a_2 * los_dtm_z;
                }
                else if(los_dtm_y==0.0){
                    // 4.2.2 - LOS is oriented north-south
                    //  p_2[1] = p_1[1] ;
                    //       - LOS orientation test
                    if(los_dtm_x<0.0){
                        // 4.2.2.1 - LOS goes due north
                        p_2[0] = col_c;
                        col_c -= 1.0;
                    }
                    else{
                        // 4.2.2.2 - LOS goes due south
                        col_c += 1.0;
                        p_2[0] = col_c;
                    }

                    a_2 = (p_2[0] - col_0)/los_dtm_x;
                    p_2[2] = z_0 + a_2 * los_dtm_z;
                }
                else{
                    // 4.2.3 - Any other LOS here
                    //            - Determination of the exit side
                    if (los_dtm_x < 0.0 && los_dtm_x <= los_dtm_y && los_dtm_x<=-los_dtm_y){
                        // 4.2.3.1 - LOS is mainly oriented north
                        //             - Intersect with north side
                        a_2 = (col_c - col_0)/los_dtm_x;
                        p_2[1] = row_0 + a_2 * los_dtm_y;

                        if(p_2[1]>row_c && p_2[1]<(row_c+1)){
                            // LOS goes out by the north
                            p_2[0] = col_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            col_c -= 1.0;
                        }

                        else if (p_2[1]<row_c){
                            // LOS goes out by the west
                            a_2 = (row_c - row_0) / los_dtm_y;
                            p_2[0] = col_0 + a_2 * los_dtm_x;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            row_c -= 1.0;
                        }

                        else if(p_2[1] > (row_c + 1)){
                            // LOS goes out by the east
                            row_c += 1.0;
                            a_2 = (row_c - row_0) / los_dtm_y;
                            p_2[0] = col_0 + a_2 * los_dtm_x;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;      
                        }

                        else if(p_2[1]==row_c){
                            // LOS goes out by the north-west corner
                            p_2[0] = col_c;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            col_c -= 1.0;
                            row_c -= 1.0;
                        }

                        else if(p_2[1] == (row_c + 1)){
                            // LOS goes out by the north-east corner
                            p_2[0] = col_c;
                            row_c += 1.0;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            col_c -= 1.0;       
                        }
                    }
                    else if(los_dtm_y>0.0 && los_dtm_y>=los_dtm_x && los_dtm_y >=-los_dtm_x){
                        // 4.2.3.2 - LOS is mainly oriented east
                        //         - Intersect with east side
                        a_2 = (row_c + 1 - row_0) / los_dtm_y;
                        p_2[0] = col_0 + a_2 * los_dtm_x;

                        if(p_2[0] > col_c && p_2[0] < (col_c + 1)){
                            //  LOS goes out by the east
                            row_c += 1.0;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;

                        }
                        else if (p_2[0] < col_c){
                            //  LOS goes out by the north
                            p_2[0] = col_c;
                            a_2 = (col_c - col_0) / los_dtm_x;
                            p_2[1] = row_0 + a_2 * los_dtm_y;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            col_c -= 1.0;

                        }
                        else if(p_2[0] > (col_c + 1)){
                            // LOS goes out by the south
                            col_c += 1.0;
                            p_2[0] = col_c;
                            a_2 = (col_c - col_0) / los_dtm_x;
                            p_2[1] = row_0 + a_2 * los_dtm_y;
                            p_2[2] = z_0 + a_2 * los_dtm_z;

                        }

                        else if(p_2[0] == col_c){
                            // LOS goes out by the north-east corner
                            row_c += 1.0;
                            p_2[0] = col_c;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            col_c -= 1.0;

                        }

                        else if(p_2[0] == (col_c + 1)){
                            // LOS goes out by the south-east corner
                            col_c += 1.0;
                            row_c += 1.0;
                            p_2[0] = col_c;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;

                        }
                    }
                    else if(los_dtm_x > 0.0 && los_dtm_x >= los_dtm_y && los_dtm_x >= -los_dtm_y){
                        // 4.2.3.3 - LOS is mainly oriented south
                        //         - Intersect with south side
                        a_2 = (col_c + 1.0 - col_0) / los_dtm_x;
                        p_2[1] = row_0 + a_2 * los_dtm_y;

                        if(p_2[1] > row_c && p_2[1] < (row_c + 1.0)){
                            // LOS goes out by the south
                            col_c += 1.0;
                            p_2[0] = col_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;

                        }
                        else if(p_2[1] < row_c){
                            // LOS goes out by the west
                            a_2 = (row_c - row_0) / los_dtm_y;
                            p_2[0] = col_0 + a_2 * los_dtm_x;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            row_c -= 1.0;

                        }
                        else if(p_2[1] > row_c + 1.0){
                            // LOS goes out by the east
                            row_c += 1.0;
                            a_2 = (row_c - row_0) / los_dtm_y;
                            p_2[0] = col_0 + a_2 * los_dtm_x;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;

                        }
                        else if(p_2[1] == row_c){
                            // LOS goes out by the south-west corner
                            col_c += 1.0;
                            p_2[0] = col_c;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            row_c -= 1.0;
                        }

                        else if(p_2[1] == row_c + 1.0){
                            // LOS goes out by the south-east corner
                            col_c += 1.0;
                            row_c += 1.0;
                            p_2[0] = col_c;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;

                        }
                    }else if(los_dtm_y < 0.0 && los_dtm_y <= los_dtm_x && los_dtm_y <= -los_dtm_x){
                        //  4.2.3.4 - VLOS is mainly oriented west
                        //          - Intersect with west side
                        a_2 = (row_c - row_0) / los_dtm_y;
                        p_2[0] = col_0 + a_2 * los_dtm_x;

                        if(p_2[0] > col_c && p_2[0] < col_c + 1.0){
                            //  LOS goes out by the west
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            row_c -= 1.0;

                        }
                        else if(p_2[0] < col_c){
                            //  LOS goes out by the north
                            p_2[0] = col_c;
                            a_2 = (col_c - col_0) / los_dtm_x;
                            p_2[1] = row_0 + a_2 * los_dtm_y;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            col_c -= 1.0;

                        }
                        else if(p_2[0] > (col_c + 1)){
                            //  LOS goes out by the south
                            col_c += 1.0;
                            p_2[0] = col_c;
                            a_2 = (col_c - col_0) / los_dtm_x;
                            p_2[1] = row_0 + a_2 * los_dtm_y;
                            p_2[2] = z_0 + a_2 * los_dtm_z;                       

                        }
                        else if(p_2[0] == col_c){
                            //  LOS goes out by the north-west corner
                            p_2[0] = col_c;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            col_c -= 1.0;
                            row_c -= 1.0;

                        }
                        else if(p_2[0] == (col_c + 1.0)){
                            //  LOS goes out by the south-west corner
                            col_c += 1.0;
                            p_2[0] = col_c;
                            p_2[1] = row_c;
                            p_2[2] = z_0 + a_2 * los_dtm_z;
                            row_c -= 1.0;                         

                        }
                    }
                }

                // LDD - min and max bounds of the "layer" Checking
                bool b_intersect = false;
                if(p_2[2] > z_0){
                    // We've gone too high, and that's not good!!!
                    b_intersect = !(((p_1[2] > h_s) && (z_0 > h_s)) || \
                                    ((p_1[2] < h_i) && (z_0 < h_i)));
                }
                else if(p_2[2] < z_1){
                    // We went too low, and that's not good either!!!
                    // (even if it already makes more sense)
                    b_intersect = !(((p_1[2] > h_s) && (z_1 > h_s)) || \
                                    ((p_1[2] < h_i) && (z_1 < h_i)));
                }
                else{
                    b_intersect = !(((p_1[2] > h_s) && (p_2[2] > h_s)) || \
                                    ((p_1[2] < h_i) && (p_2[2] < h_i))); 
                }

                // 5. LOS intersection test with the cube
                if(b_intersect){

                    // There is intersection between LOS and the cube
                    // 5.1 - DTM Altitudes
                    alti_1 = interpolate(p_1[0], p_1[1]);
                    double h_2 = interpolate(p_2[0], p_2[1]);

                    // 5.2 - Altitude differences with DTM
                    double d_alti_1 = p_1[2] - alti_1;
                    double d_2 = p_2[2] - h_2;

                    // 5.3 - Intersection test with DTM
                    if(d_alti_1 * d_2 <= 0.0){

                        // There is intersection between los and the DTM
                        // 5.3.1 - Compute of approximate solution
                        d_2 = 2.0 * m_tol_z;
                        double col_a = p_2[0];
                        double row_a = p_2[1];
                        double z_a = h_2;

                        double c_h;
                        double z_v;
                        while(abs(d_2)>m_tol_z){

                            // 5.3.1.1 - Linear interpolation coefficient of H
                            c_h = (p_1[2] - alti_1) / ((h_2 - alti_1) - (p_2[2] - p_1[2]));

                            // 5.3.1.2 - position of the interpolated point
                            col_a = p_1[0] + c_h * (p_2[0] - p_1[0]);
                            row_a = p_1[1] + c_h * (p_2[1] - p_1[1]);
                            z_a = p_1[2] + c_h * (p_2[2] - p_1[2]);

                            // 5.3.1.3 - Altitude of the interpolated point
                            z_v = interpolate(col_a, row_a);

                            // 5.3.1.4 - Altitude difference of the interpolated point
                            d_2 = z_v - z_a;

                            // 5.3.1.5 - Update
                            if (d_2<0.0){
                                // Update of the top point
                                p_1[0] = col_a;
                                p_1[1] = row_a;
                                p_1[2] = z_a;
                                alti_1 = z_v;                            
                            }
                            else{
                                // Update of the low point
                                p_2[0] = col_a;
                                p_2[1] = row_a;
                                p_2[2] = z_a;
                                h_2 = z_v;                            
                            }
                        }
                        // End, return
                        p_1[0] = col_a;
                        p_1[1] = row_a;
                        p_1[2] = z_a;

                        point_r = index_to_ter(p_1);
                        
                        double point_r_x = point_r[0];
                        double point_r_y = point_r[1];
                        double point_r_z = point_r[2];

                        return make_tuple(true, point_r_x, point_r_y, point_r_z);
                    }
                }
            }

            // End loop on meshes
            // Test if we are still in the DTM cube si on est toujours dans le cube DTM
            if(a_2>=1.0){
                // Change of plane
                i_0 +=1;

                // Loading into p_2 of the new vertex
                p_2[0] = los_x_index[i_0];
                p_2[1] = los_y_index[i_0];
                p_2[2] = los_z_index[i_0];
                h_intersect_p2 = alti[i_0];

            }
            else{
                // LDD - We looped on the meshes, 
                // we found nothing and we did not reach the next plane
                // It means we're getting out of the grip, no need to continue
                double point_r_x = point_r[0];
                double point_r_y = point_r[1];
                double point_r_z = point_r[2];

                return make_tuple(false, point_r_x, point_r_y, point_r_z);           
            }
        }// End of general case (LOS not vertical)

    }

    // End loop on the vertices
    // End, return
    double point_r_x = point_r[0];
    double point_r_y = point_r[1];
    double point_r_z = point_r[2];

   return make_tuple(false, point_r_x, point_r_y, point_r_z);
}

//-- function --//

tuple<vector<double>,
vector<double>> init_min_max(vector<double> const& alt_data,int nb_rows,int nb_columns)
{

    vector<double> alt_min_cell ((nb_rows-1)*(nb_columns-1));
    vector<double> alt_max_cell ((nb_rows-1)*(nb_columns-1));

    for(int i = 0; i < nb_rows-1; ++i){
        for(int j = 0; j < nb_columns-1; ++j){

            double val_sub_array_flat_1 = alt_data[nb_columns * i + j];
            double val_sub_array_flat_2 = alt_data[nb_columns * (i+1) + j];
            double val_sub_array_flat_3 = alt_data[nb_columns * i + (j+1)];
            double val_sub_array_flat_4 = alt_data[nb_columns * (i+1) + (j+1)];

            alt_min_cell[(nb_columns-1)*i+j] = floor(min({val_sub_array_flat_1,
                                                            val_sub_array_flat_2,
                                                            val_sub_array_flat_3,
                                                            val_sub_array_flat_4}));

            alt_max_cell[(nb_columns-1)*i+j] = ceil(max({val_sub_array_flat_1,
                                                            val_sub_array_flat_2,
                                                            val_sub_array_flat_3,
                                                            val_sub_array_flat_4}));
            
        }
    }

return make_tuple(alt_min_cell,alt_max_cell);
}



py::array_t<double> DTMIntersection::intersection_n_los_dtm(
    py::array_t<double, py::array::c_style | py::array::forcecast> los_input
    ) const
{

    py::buffer_info buf = los_input.request();
    double* los = static_cast<double*>(buf.ptr);

    int nb_points = buf.shape[0];
    int nb_alt = buf.shape[1];

    vector<double> los_i_x(nb_alt);
    vector<double> los_i_y(nb_alt);
    vector<double> los_i_z(nb_alt);

    bool var;//garbadge variable
    bool solution;
    array<double,3> position_cube;
    double alti;
    vector<double> los_index_x (2);
    vector<double> los_index_y (2);
    vector<double> los_index_z (2);
    double position_x;
    double position_y;
    double position_z;
    vector<double> res_lon(nb_points);
    vector<double> res_lat(nb_points);
    vector<double> res_alt(nb_points);


    for(int i = 0;i<nb_points;++i){

        for(int alt =0;alt<nb_alt;++alt){
            los_i_x[alt] = los[i * buf.shape[1] * buf.shape[2] + alt * buf.shape[2] + 0];
            los_i_y[alt] = los[i * buf.shape[1] * buf.shape[2] + alt * buf.shape[2] + 1];
            los_i_z[alt] = los[i * buf.shape[1] * buf.shape[2] + alt * buf.shape[2] + 2];
        }

        tie(solution, position_cube, alti, los_index_x, los_index_y, los_index_z) =\
        intersect_dtm_cube(los_i_x, los_i_y, los_i_z);
        
        if(solution){
            tie(var, position_x, position_y, position_z) =\
            intersection(los_index_x, los_index_y, los_index_z, position_cube, alti);
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

    //Cast output into np.array

    vector<double> res;
    res.reserve(nb_points*3);

    res.insert(res.end(), res_lon.begin(), res_lon.end());
    res.insert(res.end(), res_lat.begin(), res_lat.end());
    res.insert(res.end(), res_alt.begin(), res_alt.end());


    const size_t rows = nb_points;
    const size_t cols = 3;


    auto result = py::array_t<double>(
        {rows, cols},
        {sizeof(double), rows * sizeof(double)}, // strides
        res.data()
    );
    
    return result;
}