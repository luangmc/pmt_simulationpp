// simulation.h
#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include <vector>
#include <map>
#include <variant>
#include "pmt_hits.h"

using namespace std;

class SignalSimulation {
public:
    SignalSimulation(const map<string, map<string, PMTData>>& hits_dict,
                     const map<string, string>& options);
    vector<double> signal_time();
    void set_t();
    vector<double> load_txt_array(const string& filename);
    vector<double> compute_noise(const vector<double>& psd, const string& digitizer);
    void gen_noise();
    double fwhm2std(double fwhm);
    double transit_time();
    vector<double> expgaussian(const vector<double>& x, double G, double cen, double sig, double lambda);
    void spe_signal(double arr_time, vector<double>& fast_wf_aux, vector<double>& slow_wf_aux);
    void gen_signal(int nr_photons, double arr_time, vector<double>& fast_wf_aux, vector<double>& slow_wf_aux);
    vector<double> quantization(const vector<double>& signal);
    void pmt_signal(const string& pmt, const vector<string>& voxel_keys, const vector<double>& arrival_time,
                    map<string, vector<double>>& fast_signal,
                    map<string, vector<double>>& slow_signal);
    void simulated_signals(map<string, vector<double>>& fast_signal,
                           map<string, vector<double>>& slow_signal);

private:
    map<string, map<string, PMTData>> ptc_hits_;
    map<string, string> options_;
    string digitizers_;
    int fast_window_len_;
    int slow_window_len_;
    double Fs_fast_;
    double Fs_slow_;
    vector<int> sample_fast_, sample_slow_;
    vector<double> t_fast_, t_slow_;
    map<string, vector<double>> fast_noise_, slow_noise_;
};

#endif
