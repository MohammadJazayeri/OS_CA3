#include "filters_p.hpp"
#include <cmath>
#include <pthread.h>

using namespace std;


void* apply_notch_filter_chunk(void* args) {
    Filter_chunk_args* params = (Filter_chunk_args*)args;
    for (size_t k = params->start_idx; k < params->end_idx; ++k) {
        float Hf = 1.0f / (pow((*params->input_data)[k] / params->f0, 2 * params->n) + 1.0f);
        (*params->output_data)[k] = Hf;
    }
    return nullptr;
}

void* apply_band_pass_filter_chunk(void* args) {
    Filter_chunk_args* params = (Filter_chunk_args*)args;
    for (size_t k = params->start_idx; k < params->end_idx; ++k) {
        float Hf = ((*params->input_data)[k] * (*params->input_data)[k]) / ((*params->input_data)[k] * (*params->input_data)[k] + params->delta_f * params->delta_f);
        (*params->output_data)[k] = Hf;
    }
    return nullptr;
}

void* apply_FIR_filter_chunk(void* args) {
    Filter_chunk_args* params = (Filter_chunk_args*)args;
    size_t M = params->h.size();
    for (size_t n = params->start_idx; n < params->end_idx; ++n) {
        (*params->output_data)[n] = 0;
        for (size_t k = 0; k < M; ++k) {
            if (n >= k) {
                (*params->output_data)[n] += params->h[k] * (*params->input_data)[n - k];
            }
        }
    }
    return nullptr;
}

void* apply_IIR_filter_feedforward_chunk(void* args) {
    Filter_chunk_args* params = (Filter_chunk_args*)args;
    size_t M = params->b.size();

    for (size_t n = params->start_idx; n < params->end_idx; ++n) {
        (*params->output_data)[n] = 0;
        for (size_t k = 0; k < M; ++k) {
            if (n >= k) {
                (*params->output_data)[n] += params->b[k] * (*params->input_data)[n - k];
            }
        }
    }

    return nullptr;
}

void sequential_feedback_IIR(vector<float>& output_data, const vector<float>& a) {
    size_t N = a.size();
    for (size_t n = 0; n < output_data.size(); ++n) {
        for (size_t j = 1; j < N; ++j) {
            if (n >= j) {
                output_data[n] -= a[j] * output_data[n - j];
            }
        }
    }
}


void divide_and_run_filter(void* (*filter_function)(void*), vector<float>& input_data, vector<float>& output_data, Filter_chunk_args& baseArgs, int num_threads) {
    pthread_t threads[num_threads];
    Filter_chunk_args args[num_threads];
    size_t chunk_size = input_data.size() / num_threads;

    for (int i = 0; i < num_threads; ++i) {
        args[i] = baseArgs;
        args[i].start_idx = i * chunk_size;
        args[i].end_idx = (i == num_threads - 1) ? input_data.size() : (i + 1) * chunk_size;
        args[i].input_data = &input_data;
        args[i].output_data = &output_data;
        pthread_create(&threads[i], nullptr, filter_function, &args[i]);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], nullptr);
    }
}

void divide_and_run_IIR_filter(vector<float>& input_data, vector<float>& output_data, Filter_chunk_args& baseArgs, int num_threads) {

    divide_and_run_filter(apply_IIR_filter_feedforward_chunk, input_data, output_data, baseArgs, num_threads);

    sequential_feedback_IIR(output_data, baseArgs.a);
}