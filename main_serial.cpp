#include <iostream>
#include <sndfile.h>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <pthread.h>
#include <chrono>
#include <random>
using namespace std::chrono;
using namespace std;


void generate_random_floats(vector <float>& input, int size, float low, float high) {

    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<float> distribution(low, high);

    for (int i = 0; i < size; ++i) {
        float random_val = distribution(generator);
       input.push_back(random_val);
    }

    return;
}

vector<float> applyBandPassFilter(vector<float> audioData, float deltaF) {
    size_t N = audioData.size();
    vector<float> filteredData(audioData.size());

    for (size_t k = 0; k < N; ++k) {
        float Hf = (audioData[k] * audioData[k]) / (audioData[k] * audioData[k] + deltaF * deltaF);
        filteredData[k] = Hf;
    }

    return filteredData;
}

vector<float> applyNotchFilter(vector<float> audioData, float f0, int n) {
    size_t N = audioData.size();
    vector<float> filteredData(audioData.size());
    for (size_t k = 0; k < N; ++k) {
        float Hf = 1.0f / (pow(audioData[k] / f0, 2 * n) + 1.0f);
        filteredData[k] = Hf;
    }

    return filteredData;
}

vector<float> applyFIRFilter(const vector<float> input, const vector<float>& h) {
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

vector<float> applyIIRFilter(const vector<float> input, const vector<float>& a, const vector<float>& b) {
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

void readWavFile(const std::string& inputFile, std::vector<float>& data, SF_INFO& fileInfo) {
    SNDFILE* inFile = sf_open(inputFile.c_str(), SFM_READ, &fileInfo);
    if (!inFile) {
        std::cerr << "Error opening input file: " << sf_strerror(NULL) << std::endl;
        exit(1);
    }

    data.resize(fileInfo.frames * fileInfo.channels);
    sf_count_t numFrames = sf_readf_float(inFile, data.data(), fileInfo.frames);
    if (numFrames != fileInfo.frames) {
        std::cerr << "Error reading frames from file." << std::endl;
        sf_close(inFile);
        exit(1);
    }

    sf_close(inFile);
    std::cout << "Successfully read " << numFrames << " frames from " << inputFile << std::endl;
}

void writeWavFile(const std::string& outputFile, const std::vector<float>& data,  SF_INFO fileInfo) {
    sf_count_t originalFrames = fileInfo.frames;
    SNDFILE* outFile = sf_open(outputFile.c_str(), SFM_WRITE, &fileInfo);
    if (!outFile) {
        std::cerr << "Error opening output file: " << sf_strerror(NULL) << std::endl;
        exit(1);
    }

    sf_count_t numFrames = sf_writef_float(outFile, data.data(), originalFrames);
    if (numFrames != originalFrames) {
        std::cerr << "Error writing frames to file." << std::endl;
        sf_close(outFile);
        exit(1);
    }

    sf_close(outFile);
    std::cout << "Successfully wrote " << numFrames << " frames to " << outputFile << std::endl;
}

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
    vector<float> audioData, a, b, h;
    memset(&fileInfo, 0, sizeof(fileInfo));

    generate_random_floats(a, 10000, 0.0f, 4.0f);
    generate_random_floats(b, 10000, 0.0f, 4.0f);
    generate_random_floats(h, 10000, 0.0f, 4.0f);

    auto start = high_resolution_clock::now();

    readWavFile(inputFile, audioData, fileInfo);

    auto end = high_resolution_clock::now();

    auto read_duration = duration_cast<milliseconds>(end - start);

    vector<float> filteredData_notch(audioData.size());
    vector<float> filteredData_band(audioData.size());
    vector<float> filteredData_FIR(audioData.size());
    vector<float> filteredData_IIR(audioData.size());

    float deltaF = 0.1f;
    float f0 = 0.00001f;
    int n = 2;

    start = high_resolution_clock::now();

    filteredData_notch = applyNotchFilter(audioData, f0, n);

    end = high_resolution_clock::now();
    auto notch_duration = duration_cast<milliseconds>(end - start);

    start = high_resolution_clock::now();

    filteredData_band = applyBandPassFilter(audioData, deltaF);

    end = high_resolution_clock::now();
    auto band_duration = duration_cast<milliseconds>(end - start);

    start = high_resolution_clock::now();

    filteredData_FIR = applyFIRFilter(audioData, h);

    end = high_resolution_clock::now();
    auto FIR_duration = duration_cast<milliseconds>(end - start);

    start = high_resolution_clock::now();

    filteredData_IIR = applyIIRFilter(audioData, a, b);

    end = high_resolution_clock::now();
    auto IIR_duration = duration_cast<milliseconds>(end - start);
    
    writeWavFile(outputFile_band, filteredData_band, fileInfo);
    writeWavFile(outputFile_notch, filteredData_notch, fileInfo);
    writeWavFile(outputFile_FIR, filteredData_FIR, fileInfo);
    writeWavFile(outputFile_IIR, filteredData_IIR, fileInfo);

    cout << "Read: " << read_duration.count() << " ms" << endl;
    cout << "Band-pass Filter: " << band_duration.count() << " ms" << endl;
    cout << "Notch Filter: " << notch_duration.count() << " ms" << endl;
    cout << "IRR Filter: " << IIR_duration.count() << " ms" << endl;
    cout << "FIR Filter: " << FIR_duration.count() << " ms" << endl;
    cout << "Execution: " << read_duration.count() +  band_duration.count() + notch_duration.count() + IIR_duration.count() + FIR_duration.count()<< " ms" << endl;
    return 0;
}