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


// void generate_random_floats(vector <float>& input, int size, float low, float high) {

//     random_device rd;
//     mt19937 generator(rd());
//     uniform_real_distribution<float> distribution(low, high);

//     for (int i = 0; i < size; ++i) {
//         float random_val = distribution(generator);
//        input.push_back(random_val);
//     }

//     return;
// }

// vector<float> apply_band_pass_filter(vector<float> audio_data, float delta_f) {
//     size_t N = audio_data.size();
//     vector<float> filtered_data(audio_data.size());

//     for (size_t k = 0; k < N; ++k) {
//         float Hf = (audio_data[k] * audio_data[k]) / (audio_data[k] * audio_data[k] + delta_f * delta_f);
//         filtered_data[k] = Hf;
//     }

//     return filtered_data;
// }

// vector<float> apply_notch_filter(vector<float> audio_data, float f0, int n) {
//     size_t N = audio_data.size();
//     vector<float> filtered_data(audio_data.size());
//     for (size_t k = 0; k < N; ++k) {
//         float Hf = 1.0f / (pow(audio_data[k] / f0, 2 * n) + 1.0f);
//         filtered_data[k] = Hf;
//     }

//     return filtered_data;
// }

// vector<float> apply_FIR_filter(const vector<float> input, const vector<float>& h) {
//     size_t M = h.size();
//     size_t N = input.size();
//     std::vector<float> output(N);

//     for (size_t n = 0; n < N; ++n) {
//         output[n] = 0;
//         for (size_t k = 0; k < M; ++k) {
//             if (n >= k) {
//                 output[n] += h[k] * input[n - k];
//             }
//         }
//     }
//     return output;
// }

// vector<float> apply_IIR_filter(const vector<float> input, const vector<float>& a, const vector<float>& b) {
//     size_t M = b.size();
//     size_t N = a.size();
//     size_t L = input.size();
//     std::vector<float> output(L);

//     for (size_t n = 0; n < L; ++n) {
//         output[n] = 0;

//         for (size_t k = 0; k < M; ++k) {
//             if (n >= k) {
//                 output[n] += b[k] * input[n - k];
//             }
//         }

//         for (size_t j = 1; j < N; ++j) {
//             if (n >= j) {
//                 output[n] -= a[j] * output[n - j];
//             }
//         }
//     }
//     return output;
// }

// void readWavFile(const std::string& inputFile, std::vector<float>& data, SF_INFO& fileInfo) {
//     SNDFILE* inFile = sf_open(inputFile.c_str(), SFM_READ, &fileInfo);
//     if (!inFile) {
//         std::cerr << "Error opening input file: " << sf_strerror(NULL) << std::endl;
//         exit(1);
//     }

//     data.resize(fileInfo.frames * fileInfo.channels);
//     sf_count_t numFrames = sf_readf_float(inFile, data.data(), fileInfo.frames);
//     if (numFrames != fileInfo.frames) {
//         std::cerr << "Error reading frames from file." << std::endl;
//         sf_close(inFile);
//         exit(1);
//     }

//     sf_close(inFile);
//     std::cout << "Successfully read " << numFrames << " frames from " << inputFile << std::endl;
// }

// void writeWavFile(const std::string& outputFile, const std::vector<float>& data,  SF_INFO fileInfo) {
//     sf_count_t originalFrames = fileInfo.frames;
//     SNDFILE* outFile = sf_open(outputFile.c_str(), SFM_WRITE, &fileInfo);
//     if (!outFile) {
//         std::cerr << "Error opening output file: " << sf_strerror(NULL) << std::endl;
//         exit(1);
//     }

//     sf_count_t numFrames = sf_writef_float(outFile, data.data(), originalFrames);
//     if (numFrames != originalFrames) {
//         std::cerr << "Error writing frames to file." << std::endl;
//         sf_close(outFile);
//         exit(1);
//     }

//     sf_close(outFile);
//     std::cout << "Successfully wrote " << numFrames << " frames to " << outputFile << std::endl;
// }

// void print_info(int read_duration, int notch_duration, int band_duration, int FIR_duration, int IIR_duration){
//     cout << "Read: " << read_duration << " ms" << endl;
//     cout << "Notch Filter Time: " << notch_duration << " ms" << endl;
//     cout << "Band-Pass Filter Time: " << band_duration << " ms" << endl;
//     cout << "FIR Filter Time: " << FIR_duration << " ms" << endl;
//     cout << "IIR Filter Time: " << IIR_duration << " ms" << endl;

//     auto total_duration = notch_duration + band_duration + FIR_duration + IIR_duration;
//     cout << "Total Execution Time: " << total_duration << " ms" << endl;
// }

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