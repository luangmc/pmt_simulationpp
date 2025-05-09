// pmt_hits.h
#ifndef PMT_HITS_H
#define PMT_HITS_H

#include <string>
#include <vector>
#include <map>
#include <nlohmann/json.hpp>

using namespace std;

struct PMTData {
    int hits;
    double arrival_time;
};

class PhotonPropagation {
public:
    PhotonPropagation(const vector<double>& x_0,
                      const vector<double>& y_0,
                      const vector<int>& n_photons,
                      const vector<double>& arr_times,
                      const map<string, string>& options);

    map<string, map<string, int>> sim_pmt_hits_with_eq(double x, double y, int n_photons);

    map<string, map<string, PMTData>> pmt_hits();

private:
    map<string, string> options_;
    vector<double> x_0_, y_0_;
    vector<int> n_photons_;
    vector<double> arr_times_;
    double z_0_ = 0.0;
};

#endif
