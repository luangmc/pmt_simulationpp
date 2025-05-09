// simulation.cxx
#include "signal_generation.h"
#include <fstream>
#include <cmath>
#include <random>
#include <complex>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/core/utility.hpp>

using namespace cv;
using namespace std;
using namespace chrono;

SignalSimulation::SignalSimulation(const map<string, map<string, PMTData>>& hits_dict,
                                   const map<string, string>& options)
    : ptc_hits_(hits_dict), options_(options) {
    digitizers_ = options_["digitizers"];
    fast_window_len_ = stoi(options_["fast_window_len"]);
    slow_window_len_ = stoi(options_["slow_window_len"]);
    Fs_fast_ = stod(options_["fast_freq"]);
    Fs_slow_ = stod(options_["slow_freq"]);
    set_t();
    gen_noise();
}

vector<double> SignalSimulation::signal_time() {
    vector<double> times;
    for (const auto& [voxel, data] : ptc_hits_) {
        for (const auto& [pmt, pmt_data] : data) {
            times.push_back(pmt_data.arrival_time);
        }
    }
    return times;
}

void SignalSimulation::set_t() {
    t_fast_.resize(fast_window_len_);
    t_slow_.resize(slow_window_len_);
    for (int i = 0; i < fast_window_len_; ++i) t_fast_[i] = i / Fs_fast_;
    for (int i = 0; i < slow_window_len_; ++i) t_slow_[i] = i / Fs_slow_;
}

double SignalSimulation::fwhm2std(double fwhm) {
    return fwhm / (2.0 * sqrt(2.0 * log(2.0)));
}

double SignalSimulation::transit_time() {
    double mu = stod(options_["transit_time"]);
    double fwhm = stod(options_["transit_time_spread"]);
    double sigma = fwhm2std(fwhm);
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> d(mu, sigma);
    return d(gen);
}

vector<double> SignalSimulation::expgaussian(const vector<double>& x, double G, double cen, double sig, double lambda) {
    vector<double> wf(x.size());
    double Q = G * (-1.6e-19);
    double R = 50.0;

    for (size_t i = 0; i < x.size(); ++i) {
        double exp_arg = lambda * (cen - x[i] + (lambda * sig * sig / 2.0));
        double erfc_arg = (cen + lambda * sig * sig - x[i]) / (sig * sqrt(2.0));
    
        if (!isfinite(exp_arg) || !isfinite(erfc_arg) || exp_arg > 700.0 || abs(erfc_arg) > 30.0) {
            wf[i] = 0.0;
            continue;
        }

        double pulse = (lambda / 2.0) * exp(exp_arg) * erfc(erfc_arg);

        if (!isfinite(pulse)) {
            wf[i] = 0.0;
        } else {
            wf[i] = Q * R * pulse;
        }
    }
    return wf;
}

void SignalSimulation::spe_signal(double arr_time, vector<double>& fast_wf_aux, vector<double>& slow_wf_aux) {
    double gain = stod(options_["pmt_gain"]), std = stod(options_["pmt_sigma"]), lambda = stod(options_["pmt_lambda"]);
    double scale = stod(options_["exp_dispersion_scale"]);

    random_device rd;
    mt19937 gen(rd());
    exponential_distribution<double> exp_dist(scale);

    double disp = exp_dist(gen) * 1e9;
    double mean = (transit_time() + arr_time + disp) * 1e-9;
    
    float shift_fast = 200 / Fs_fast_;
    float shift_slow = 1466 / Fs_slow_;

    if (digitizers_ == "Both" || digitizers_ == "Fast") {
        auto wf = expgaussian(t_fast_, gain, mean + shift_fast, std, lambda);
        for (size_t i = 0; i < wf.size(); ++i) fast_wf_aux[i] += wf[i];
    }
    if (digitizers_ == "Both" || digitizers_ == "Slow") {
        auto wf = expgaussian(t_slow_, gain, mean + shift_slow, std, lambda);
        for (size_t i = 0; i < wf.size(); ++i) slow_wf_aux[i] += wf[i];
    }
}

void SignalSimulation::gen_signal(int nr_photons, double arr_time, vector<double>& fast_wf_aux, vector<double>& slow_wf_aux) {
    for (int i = 0; i < nr_photons; ++i)
        spe_signal(arr_time, fast_wf_aux, slow_wf_aux);
}

vector<double> SignalSimulation::quantization(const vector<double>& signal) {
    int n_bits = 12;
    int levels = 1 << n_bits;

    vector<double> scaled_signal(signal.size());
    for (size_t i = 0; i < signal.size(); ++i) {
        scaled_signal[i] = round(signal[i] * levels); 
    }

    double min_val = *min_element(scaled_signal.begin(), scaled_signal.end());
    double max_val = *max_element(scaled_signal.begin(), scaled_signal.end());
    double interval = (max_val - min_val) / (levels - 1);

    vector<double> result(signal.size());
    for (size_t i = 0; i < scaled_signal.size(); ++i) {
        double quantized = round((scaled_signal[i] - min_val) / interval);
        quantized = clamp(quantized, 0.0, static_cast<double>(levels - 1));
        result[i] = (quantized * interval + min_val) / (levels - 1);  
    }
    return result;
}

void SignalSimulation::pmt_signal(const string& pmt, const vector<string>& voxel_keys,
                                   const vector<double>& arrival_time,
                                   map<string, vector<double>>& fast_signal,
                                   map<string, vector<double>>& slow_signal) {
    vector<double> fast_wf_aux(fast_window_len_, 0.0);
    vector<double> slow_wf_aux(slow_window_len_, 0.0);
                                   
    for (size_t i = 0; i < voxel_keys.size(); ++i) {
        const auto& pmt_data = ptc_hits_[voxel_keys[i]].at(pmt);  // Acesso estruturado
        gen_signal(pmt_data.hits, arrival_time[i], fast_wf_aux, slow_wf_aux);
    }

    for (size_t i = 0; i < fast_wf_aux.size(); ++i) fast_wf_aux[i] += fast_noise_[pmt][i];
    for (size_t i = 0; i < slow_wf_aux.size(); ++i) slow_wf_aux[i] += slow_noise_[pmt][i];

    //fast_signal[pmt] = fast_wf_aux;
    //slow_signal[pmt] = slow_wf_aux;
    fast_signal[pmt] = quantization(fast_wf_aux);
    slow_signal[pmt] = quantization(slow_wf_aux);
}

void SignalSimulation::simulated_signals(map<string, vector<double>>& fast_signal,
                                         map<string, vector<double>>& slow_signal) {   
    auto start = high_resolution_clock::now();
                                                                             
    vector<string> pmts = {"pmt_1", "pmt_2", "pmt_3", "pmt_4"};
    vector<string> voxel_keys;
    for (const auto& [key, _] : ptc_hits_) voxel_keys.push_back(key);
    vector<double> arrival_time = signal_time();

    for (const auto& pmt : pmts)
        pmt_signal(pmt, voxel_keys, arrival_time, fast_signal, slow_signal);

    fast_signal["time"] = t_fast_;
    slow_signal["time"] = t_slow_;
    
    auto end = high_resolution_clock::now();
    cout << "Signal generation took "
         << duration_cast<milliseconds>(end - start).count()
         << " ms" << endl;
}

vector<double> SignalSimulation::load_txt_array(const string& filename) {
    ifstream file(filename);
    vector<double> data;
    double val;
    while (file >> val) {
        data.push_back(val);
    }
    return data;
}

void SignalSimulation::gen_noise() {
    auto start = high_resolution_clock::now();
    vector<string> pmts = {"pmt_1", "pmt_2", "pmt_3", "pmt_4"};
    for (const auto& pmt : pmts) {
        vector<double> fast_psd = load_txt_array(options_["fast_noise_path_" + pmt]);
        vector<double> slow_psd = load_txt_array(options_["slow_noise_path_" + pmt]);

        //fast_noise_[pmt] = vector<double>(fast_psd.size(), 0.0);
        //slow_noise_[pmt] = vector<double>(slow_psd.size(), 0.0);

        fast_noise_[pmt] = compute_noise(fast_psd, "Fast");
        slow_noise_[pmt] = compute_noise(slow_psd, "Slow");
    }

    auto end = high_resolution_clock::now();
    cout << "Gen noise took "
         << duration_cast<milliseconds>(end - start).count()
         << " ms" << endl;
}

vector<double> SignalSimulation::compute_noise(const vector<double>& psd, const string& digitizer) {
    int N = (digitizer == "Fast") ? fast_window_len_ : slow_window_len_;
    double fs = (digitizer == "Fast") ? Fs_fast_ : Fs_slow_;

    vector<double> PSD = psd;
    for (int i = psd.size() - 2; i > 0; --i) {
        PSD.push_back(psd[i]);
    }

    size_t len = PSD.size();
    vector<complex<double>> Nf(len);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> phase_dist(-M_PI, M_PI);

    // Positive spectrum
    for (size_t i = 0; i < psd.size(); ++i) {
        double phase = phase_dist(gen);
        double magnitude = sqrt(psd[i] * fs / N);
        Nf[i] = polar(magnitude, phase);
    }

    // Negative spectrum
    for (size_t i = psd.size(); i < len; ++i) {
        Nf[i] = conj(Nf[len - i]);
    }

    // Convert Nf to cv::Mat for OpenCV processing
    Mat complexInput(len, 1, CV_64FC2);
    for (size_t i = 0; i < len; ++i) {
        complexInput.at<Vec2d>(i)[0] = Nf[i].real();
        complexInput.at<Vec2d>(i)[1] = Nf[i].imag();
    }

    // Perform IFFT using OpenCV
    Mat complexOutput;
    dft(complexInput, complexOutput, DFT_INVERSE | DFT_REAL_OUTPUT);

    // Convert the IFFT result back to a vector of doubles
    vector<double> noise(N);
    for (int n = 0; n < N; ++n) {
        noise[n] = complexOutput.at<double>(n);
    }

    return noise;
}
