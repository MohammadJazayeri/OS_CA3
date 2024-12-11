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
#include "filters_s.hpp"

using namespace std::chrono;
using namespace std;

int main(int argc, char* argv[]) {
    if(argc != 2){
        printf("Invalid arguments\n");
    }
    string inputFile = argv[1];
    string outputFile_band = "band_serial.wav";
    string outputFile_notch = "notch_serial.wav";
    string outputFile_IIR = "IIR_serial.wav";
    string outputFile_FIR = "FIR_serial.wav";

    SF_INFO fileInfo;
    vector<float> audio_data, a, b, h;
    memset(&fileInfo, 0, sizeof(fileInfo));

    generate_random_floats(a, 10000, 0.0f, 4.0f);
    generate_random_floats(b, 10000, 0.0f, 4.0f);
    generate_random_floats(h, 10000, 0.0f, 4.0f);

    auto start = high_resolution_clock::now();

    readWavFile(inputFile, audio_data, fileInfo);

    auto end = high_resolution_clock::now();

    auto read_duration = duration_cast<milliseconds>(end - start);

    vector<float> filtered_data_notch(audio_data.size());
    vector<float> filtered_data_band(audio_data.size());
    vector<float> filtered_data_FIR(audio_data.size());
    vector<float> filtered_data_IIR(audio_data.size());

    float delta_f = 0.1f;
    float f0 = 0.00001f;
    int n = 2;

    start = high_resolution_clock::now();

    filtered_data_notch = apply_notch_filter(audio_data, f0, n);

    end = high_resolution_clock::now();
    auto notch_duration = duration_cast<milliseconds>(end - start);

    start = high_resolution_clock::now();

    filtered_data_band = apply_band_pass_filter(audio_data, delta_f);

    end = high_resolution_clock::now();
    auto band_duration = duration_cast<milliseconds>(end - start);

    start = high_resolution_clock::now();

    filtered_data_FIR = apply_FIR_filter(audio_data, h);

    end = high_resolution_clock::now();
    auto FIR_duration = duration_cast<milliseconds>(end - start);

    start = high_resolution_clock::now();

    filtered_data_IIR = apply_IIR_filter(audio_data, a, b);

    end = high_resolution_clock::now();
    auto IIR_duration = duration_cast<milliseconds>(end - start);
    
    writeWavFile(outputFile_band, filtered_data_band, fileInfo);
    writeWavFile(outputFile_notch, filtered_data_notch, fileInfo);
    writeWavFile(outputFile_FIR, filtered_data_FIR, fileInfo);
    writeWavFile(outputFile_IIR, filtered_data_IIR, fileInfo);

    print_info(read_duration.count(), notch_duration.count(), band_duration.count(), FIR_duration.count(),IIR_duration.count());
    return 0;
}