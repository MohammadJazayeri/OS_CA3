#include "filters_s.hpp"
#include <cmath>

vector<float> apply_band_pass_filter(vector<float> audio_data, float delta_f) {
    size_t N = audio_data.size();
    vector<float> filtered_data(audio_data.size());

    for (size_t k = 0; k < N; ++k) {
        float Hf = (audio_data[k] * audio_data[k]) / (audio_data[k] * audio_data[k] + delta_f * delta_f);
        filtered_data[k] = Hf;
    }

    return filtered_data;
}

vector<float> apply_notch_filter(vector<float> audio_data, float f0, int n) {
    size_t N = audio_data.size();
    vector<float> filtered_data(audio_data.size());
    for (size_t k = 0; k < N; ++k) {
        float Hf = 1.0f / (pow(audio_data[k] / f0, 2 * n) + 1.0f);
        filtered_data[k] = Hf;
    }

    return filtered_data;
}

vector<float> apply_FIR_filter(const vector<float> input, const vector<float>& h) {
    size_t M = h.size();
    size_t N = input.size();
    std::vector<float> output(N);

    for (size_t n = 0; n < N; ++n) {
        output[n] = 0;
        for (size_t k = 0; k < M; ++k) {
            if (n >= k) {
                output[n] += h[k] * input[n - k];
            }
        }
    }
    return output;
}

vector<float> apply_IIR_filter(const vector<float> input, const vector<float>& a, const vector<float>& b) {
    size_t M = b.size();
    size_t N = a.size();
    size_t L = input.size();
    std::vector<float> output(L);

    for (size_t n = 0; n < L; ++n) {
        output[n] = 0;

        for (size_t k = 0; k < M; ++k) {
            if (n >= k) {
                output[n] += b[k] * input[n - k];
            }
        }

        for (size_t j = 1; j < N; ++j) {
            if (n >= j) {
                output[n] -= a[j] * output[n - j];
            }
        }
    }
    return output;
}