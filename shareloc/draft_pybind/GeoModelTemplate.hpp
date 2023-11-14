#include <string>
#include <iostream>
#include <vector>
#include <tuple>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::tuple;

class GeoModelTemplate {
private:
    string geomodel_path;
    string type;
public:
    GeoModelTemplate(string geomodel_path_arg);
    ~GeoModelTemplate();
    vector<vector<double>> direct_loc_h(vector<double> row, vector<double> col, double alt, bool fill_nan=false);
    vector<vector<double>> direct_loc_dtm(vector<double> row, vector<double> col, string dtm);
    tuple<vector<double>,vector<double>,vector<double>> inverse_loc(vector<double> lon, vector<double> lat, double alt);
};
