// pmt_hits.cxx
#include "pmt_hits.h"
#include <fstream>
#include <cmath>
#include <random>

using namespace std;

PhotonPropagation::PhotonPropagation(const vector<double>& x_0,
                                     const vector<double>& y_0,
                                     const vector<int>& n_photons,
                                     const vector<double>& arr_times,
                                     const map<string, string>& options)
    : x_0_(x_0), y_0_(y_0), n_photons_(n_photons), arr_times_(arr_times), options_(options) {
}

map<string, map<string, int>>
PhotonPropagation::sim_pmt_hits_with_eq(double x, double y, int n_photons) {
    map<string, map<string, int>> hits;
    random_device rd;
    mt19937 gen(rd());

    const auto& pmt_positions = options_["pmt_positions"];
    double r_pmt = stod(options_["pmt_radius"]);
    double z_pmt = stod(options_["dist_gem_pmt"]);
    double n = 3.6;

    for (int i = 1; i <= stoi(options_["pmt_number"]); ++i) {
        string pmt_name = "pmt_" + to_string(i);
        double x_pmt = stod(options_[pmt_name + "_x"]);
        double y_pmt = stod(options_[pmt_name + "_y"]);

        double R = sqrt(pow(x_pmt - x, 2) + pow(y_pmt - y, 2) + pow(z_pmt - z_0_, 2));
        double mean = n_photons * pow(r_pmt, 2) * pow(z_pmt, 2) / (4 * pow(R, n));
        poisson_distribution<> d(mean);

        hits["pmt_hits"][pmt_name] = d(gen);
    }
    return hits;
}

map<string, map<string, PMTData>> PhotonPropagation::pmt_hits() {
    map<string, map<string, PMTData>> all_hits;

    for (size_t i = 0; i < x_0_.size(); ++i) {
        string voxel_name = "voxel_" + to_string(i);
        auto hits_raw = sim_pmt_hits_with_eq(x_0_[i], y_0_[i], n_photons_[i]);

        const auto& pmt_hits = hits_raw["pmt_hits"];
        for (const auto& [pmt, hit_count] : pmt_hits) {
            PMTData data;
            data.hits = hit_count;
            data.arrival_time = arr_times_[i];
            all_hits[voxel_name][pmt] = data;
        }
    }

    return all_hits;
}
