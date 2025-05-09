// main.cxx
#include "pmt_hits.h"
#include "signal_generation.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <chrono>

using namespace std;
using namespace chrono;

void save_signals_to_csv(const map<string, vector<double>>& signals, const string& filename);

vector<double> load_txt_array(const string& filename);

vector<int> load_txt_array_int(const string& filename);

void ReadConfig(string name, map<string,string>& options);

int main() {
    auto start_total = high_resolution_clock::now();

    vector<double> x0 = load_txt_array("../6keV_arrays/x_0.txt");
    vector<double> y0 = load_txt_array("../6keV_arrays/y_0.txt");
    vector<int> photons = load_txt_array_int("../6keV_arrays/n_photons.txt");
    vector<double> arrival_times = load_txt_array("../6keV_arrays/arr_times.txt");

    cout << "Start PMT simulation" << endl;

    string config_path = "../config/ConfigFile_new.txt";
    map<string,string> options;
    ReadConfig(config_path, options);	

    // PMT Hits
    cout << "Photon Propagation - start" << endl;
    auto start_propagation = high_resolution_clock::now();

    PhotonPropagation propag(x0, y0, photons, arrival_times, options);
    auto hits_dict = propag.pmt_hits();

    auto end_propagation = high_resolution_clock::now();
    cout << "Photon Propagation - end ("
         << duration_cast<milliseconds>(end_propagation - start_propagation).count()
         << " ms)" << endl;

    // PMT Signal generation
    cout << "Signal simulation - start" << endl;
    auto start_simulation = high_resolution_clock::now();

    SignalSimulation simulator(hits_dict, options);
    map<string, vector<double>> fast_signal, slow_signal;
    simulator.simulated_signals(fast_signal, slow_signal);

    auto end_simulation = high_resolution_clock::now();
    cout << "Signal simulation - end ("
         << duration_cast<milliseconds>(end_simulation - start_simulation).count()
         << " ms)" << endl;

    auto end_total = high_resolution_clock::now();
    cout << "End PMT Simulation ("
         << duration_cast<milliseconds>(end_total - start_total).count()
         << " ms)" << endl;

    save_signals_to_csv(fast_signal, "../output/fast_signal.csv");
    save_signals_to_csv(slow_signal, "../output/slow_signal.csv");     

    return 0;
}

void save_signals_to_csv(const map<string, vector<double>>& signals, const string& filename) {
    ofstream file(filename);
    if (file.is_open()) {
        file << "KEY,Index,Signal\n";
        for (const auto& pair : signals) {
            const string& pmt = pair.first;
            const vector<double>& signal = pair.second;
            for (size_t i = 0; i < signal.size(); ++i) {
                file << pmt << "," << i << "," << signal[i] << "\n";
            }
        }
        file.close();
        cout << "Signals saved in " << filename << endl;
    } else {
        cerr << "Error csv" << endl;
    }
}

vector<double> load_txt_array(const string& filename) {
    ifstream file(filename);
    vector<double> data;

    if (!file.is_open()) {
        cerr << "Error opening the file: " << filename << "\n";
        return data;
    }

    double val;
    while (file >> val) {
        data.push_back(val);
    }

    return data;
}

vector<int> load_txt_array_int(const string& filename) {
    ifstream infile(filename);
    vector<int> data;
    
    if (!infile.is_open()) {
        cerr << "Error opening the file: " << filename << endl;
        return data;
    }

    double temp;
    while (infile >> temp) {
        int val = static_cast<int>(round(temp)); 
        data.push_back(val);
    }

    return data;
}

void ReadConfig(string name, map<string,string>& options)
{
    ifstream config(name.c_str());
	
    string line;
    while(getline(config,line))
    {
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return isspace(c);}),line.end());
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return c=='\'';}),line.end());
        if(line[0] == '#' || line.empty() || line[0] == '{' || line[0] == '}') continue;
		
        auto delim1= line.find(":");
        auto delim2= line.find(",");
        if(delim2==string::npos) delim2=line.size();
        auto index= line.substr(0,delim1);
        auto val= line.substr(delim1+1,min(delim2,line.size())-delim1-1 );
        options[index]=val;
    }
	
}
