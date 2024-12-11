#include <iostream>
#include <sndfile.h>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <pthread.h>
#include <chrono>
#include <random>
#include "read_write.hpp"
#include "filters_p.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Invalid arguments\n");
        return 1;
    }

    string inputFile = argv[1];
    string outputFile_band = "band_parallel.wav";
    string outputFile_notch = "notch_parallel.wav";
    string outputFile_IIR = "IIR_parallel.wav";
    string outputFile_FIR = "FIR_parallel.wav";

    SF_INFO fileInfo;
    vector<float> audioData, a, b, h, output_data_notch, output_data_band, output_data_FIR, output_data_IIR;
    memset(&fileInfo, 0, sizeof(fileInfo));

    generate_random_floats(a, 10000, 0.0f, 4.0f);
    generate_random_floats(b, 10000, 0.0f, 4.0f);
    generate_random_floats(h, 10000, 0.0f, 4.0f);

    auto start = high_resolution_clock::now();

    readWavFile(inputFile, audioData, fileInfo);

    auto end = high_resolution_clock::now();

    auto read_duration = duration_cast<milliseconds>(end - start);

    output_data_notch.resize(audioData.size());
    output_data_band.resize(audioData.size());
    output_data_FIR.resize(audioData.size());
    output_data_IIR.resize(audioData.size());

    float delta_f = 0.1f;
    float f0 = 0.00001f;
    int n = 2;

    Filter_chunk_args notch_args = {&audioData, &output_data_notch, a, b, h, delta_f, f0, n};
    Filter_chunk_args band_args = {&audioData, &output_data_band, a, b, h, delta_f, f0, n};
    Filter_chunk_args FIR_args = {&audioData, &output_data_FIR, a, b, h, delta_f, f0, n};
    Filter_chunk_args IIR_args = {&audioData, &output_data_IIR, a, b, h, delta_f, f0, n};

    int num_threads = 6;

    auto start_notch = high_resolution_clock::now();
    divide_and_run_filter(apply_notch_filter_chunk, audioData, output_data_notch, notch_args, num_threads);
    auto end_notch = high_resolution_clock::now();
    auto notch_duration = duration_cast<milliseconds>(end_notch - start_notch);

    auto start_band = high_resolution_clock::now();
    divide_and_run_filter(apply_band_pass_filter_chunk, audioData, output_data_band, band_args, num_threads);
    auto end_band = high_resolution_clock::now();
    auto band_duration = duration_cast<milliseconds>(end_band - start_band);

    auto start_FIR = high_resolution_clock::now();
    divide_and_run_filter(apply_FIR_filter_chunk, audioData, output_data_FIR, FIR_args, num_threads);
    auto end_FIR = high_resolution_clock::now();
    auto FIR_duration = duration_cast<milliseconds>(end_FIR - start_FIR);

    auto start_IIR = high_resolution_clock::now();
    // divide_and_run_filter(apply_IIR_filter_chunk, audioData, output_data_IIR, IIR_args, num_threads);
    divide_and_run_IIR_filter(audioData, output_data_IIR, IIR_args, num_threads);
    auto end_IIR = high_resolution_clock::now();
    auto IIR_duration = duration_cast<milliseconds>(end_IIR - start_IIR);

    writeWavFile(outputFile_band, output_data_band, fileInfo);
    writeWavFile(outputFile_notch, output_data_notch, fileInfo);
    writeWavFile(outputFile_FIR, output_data_FIR, fileInfo);
    writeWavFile(outputFile_IIR, output_data_IIR, fileInfo);

    print_info(read_duration.count(), notch_duration.count(), band_duration.count(), FIR_duration.count(),IIR_duration.count());
    return 0;
}