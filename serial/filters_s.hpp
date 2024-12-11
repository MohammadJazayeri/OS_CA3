#ifndef FILTERSS_H
#define FILTERSS_H

#include <vector>
using namespace std;


vector<float> apply_band_pass_filter(vector<float> audio_data, float delta_f);
vector<float> apply_notch_filter(vector<float> audio_data, float f0, int n);
vector<float> apply_FIR_filter(const vector<float> input, const vector<float>& h);
vector<float> apply_IIR_filter(const vector<float> input, const vector<float>& a, const vector<float>& b);

#endif