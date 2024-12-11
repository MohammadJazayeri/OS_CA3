#include "read_write.hpp"
#include <iostream>
#include <random>


void print_info(int read_duration, int notch_duration, int band_duration, int FIR_duration, int IIR_duration) {
    cout << "Read: " << read_duration << " ms" << endl;
    cout << "Notch Filter Time: " << notch_duration << " ms" << endl;
    cout << "Band-Pass Filter Time: " << band_duration << " ms" << endl;
    cout << "FIR Filter Time: " << FIR_duration << " ms" << endl;
    cout << "IIR Filter Time: " << IIR_duration << " ms" << endl;

    auto total_duration = notch_duration + band_duration + FIR_duration + IIR_duration;
    cout << "Total Execution Time: " << total_duration << " ms" << endl;
}

void writeWavFile(const string& outputFile, const vector<float>& data,  SF_INFO fileInfo) {
    sf_count_t originalFrames = fileInfo.frames;
    SNDFILE* outFile = sf_open(outputFile.c_str(), SFM_WRITE, &fileInfo);
    if (!outFile) {
        cerr << "Error opening output file: " << sf_strerror(NULL) << endl;
        exit(1);
    }

    sf_count_t numFrames = sf_writef_float(outFile, data.data(), originalFrames);
    if (numFrames != originalFrames) {
        cerr << "Error writing frames to file." << endl;
        sf_close(outFile);
        exit(1);
    }

    sf_close(outFile);
    cout << "Successfully wrote " << numFrames << " frames to " << outputFile << endl;
}

void readWavFile(const string& inputFile, vector<float>& data, SF_INFO& fileInfo) {
    SNDFILE* inFile = sf_open(inputFile.c_str(), SFM_READ, &fileInfo);
    if (!inFile) {
        cerr << "Error opening input file: " << sf_strerror(NULL) << endl;
        exit(1);
    }

    data.resize(fileInfo.frames * fileInfo.channels);
    sf_count_t numFrames = sf_readf_float(inFile, data.data(), fileInfo.frames);
    if (numFrames != fileInfo.frames) {
        cerr << "Error reading frames from file." << endl;
        sf_close(inFile);
        exit(1);
    }

    sf_close(inFile);
    cout << "Successfully read " << numFrames << " frames from " << inputFile << endl;
}

void generate_random_floats(vector<float>& input, int size, float low, float high) {
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<float> distribution(low, high);
    for (int i = 0; i < size; ++i) {
        input.push_back(distribution(generator));
    }
}