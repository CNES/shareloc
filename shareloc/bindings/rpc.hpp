#include "GeoModelTemplate.hpp"
#include "GeoModelTemplate.cpp"
#include <map>

using std::map;

/*
Squelette de la classe RPC.py de l'issue 221

*/

class RPC : public GeoModelTemplate{
private:
    string type;
    int epsg;
    string datum;
    map<string, double> rpc_params; //map<string, double> is a simple dict -> maybe inappropiate here
    double lim_extrapol;

    vector<vector<double>> monomes;
    vector<vector<double>> monomes_deriv_1;
    vector<vector<double>> monomes_deriv_2;

    bool inverse_coefficient;
    bool direct_coefficient;

    vector<double> num_col;
    vector<double> den_col;
    vector<double> num_row;
    vector<double> den_row;

    vector<double> num_x;
    vector<double> den_x;
    vector<double> num_y;
    vector<double> den_y;

    vector<double> alt_minmax;

    double col0;
    double colmax;
    double row0;
    double rowmax;

    double offset_row;
    double scale_row;
    double offset_col;
    double scale_col;
    double offset_alt;
    double scale_alt;
    
public:
    using GeoModelTemplate::GeoModelTemplate;

    vector<vector<double>> direct_loc_h(vector<double> row, vector<double> col, double alt, bool fill_nan=false);
    tuple<vector<vector<double>>,vector<vector<double>>> direct_loc_grid_h(int row0, int col0, int steprow, int stepcol, int nbrow, int nbcol, double alt);
    vector<vector<double>> direct_loc_dtm(double row, double col, string dtm);// dtm is a python class not a string
    tuple<vector<double>,vector<double>,vector<double>> inverse_loc(vector<double> lon, vector<double> lat, double alt);
    vector<vector<double>> filter_coordinates(vector<double> first_coord, vector<double> second_coord, bool fill_nan=false, string direction="direct");
    tuple<vector<double>,vector<double>,vector<double>,vector<double>> compute_loc_inverse_derivates(vector<double> lon, vector<double> lat, vector<double> alt);
    vector<vector<double>> direct_loc_inverse_iterative(vector<double> row, vector<double> col, double alt, int nb_iter_max=10, bool fill_nan=false);
    vector<double> get_alt_min_max();
    vector<vector<double>> los_extrema(double row, double col, double alt_min, double alt_max, bool fill_nan=false, int epsg=4326);
};