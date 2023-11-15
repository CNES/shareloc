


GeoModelTemplate::GeoModelTemplate() {
    cout<<"GeoModelTemplate : constructor"<<endl;
}
GeoModelTemplate::~GeoModelTemplate() {
    cout<<"GeoModelTemplate : destructor"<<endl;
}

vector<vector<double>> GeoModelTemplate::direct_loc_h(vector<double> row, vector<double> col, double alt, bool fill_nan){
    cout<<"GeoModelTemplate : direct_loc_h"<<endl;
    vector<vector<double>> vect;
    return vect;
}

vector<vector<double>> GeoModelTemplate::direct_loc_dtm(vector<double> row, vector<double> col, string dtm){
    cout<<"GeoModelTemplate : direct_loc_dtm"<<endl;
    vector<vector<double>> vect;
    return vect;
}

tuple<vector<double>,vector<double>,vector<double>> GeoModelTemplate::inverse_loc(vector<double> lon, vector<double> lat, double alt){
    cout<<"GeoModelTemplate : inverse_loc"<<endl;
    tuple<vector<double>,vector<double>,vector<double>> res;
    return res;
}