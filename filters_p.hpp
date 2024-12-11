#ifndef FILTERSP_H
#define FILTERSP_H

#include <vector>
using namespace std;

struct Filter_chunk_args {
    vector<float>* input_data;
    vector<float>* output_data;
    vector<float> a;
    vector<float> b;
    vector<float> h;
    float delta_f;
    float f0;
    int n;
    size_t start_idx;
    size_t end_idx;   
};


void* apply_notch_filter_chunk(void* args);
void* apply_band_pass_filter_chunk(void* args);
void* apply_FIR_filter_chunk(void* args);
void* apply_IIR_filter_chunk(void* args);
void divide_and_run_filter(void* (*filter_function)(void*),
 vector<float>& input_data, vector<float>& output_data, Filter_chunk_args& baseArgs, int num_threads);

#endif