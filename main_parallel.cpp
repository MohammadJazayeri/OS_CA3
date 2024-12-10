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

// Function prototypes
void* applyNotchFilterThread(void* args);
void* applyBandPassFilterThread(void* args);
void* applyFIRFilterThread(void* args);
void* applyIIRFilterThread(void* args);

// Shared data structure to pass arguments to threads
struct FilterArgs {
    vector<float>* inputData;
    vector<float> outputData;
    vector<float> coefficientsA;
    vector<float> coefficientsB;
    vector<float> coefficientsH;
    float deltaF, sampleRate, f0;
    int n;
};

// Notch Filter
void* applyNotchFilterThread(void* args) {
    FilterArgs* data = (FilterArgs*)args;
    size_t N = data->inputData->size();
    for (size_t k = 0; k < N; ++k) {
        float f = (k < N / 2 ? k : k - N) * data->sampleRate / N;
        float Hf = 1.0f / (pow(f / data->f0, 2 * data->n) + 1.0f);
        (*data->inputData)[k] *= Hf;
    }
    data->outputData = *(data->inputData);
    return nullptr;
}

// Band-pass Filter
void* applyBandPassFilterThread(void* args) {
    FilterArgs* data = (FilterArgs*)args;
    size_t N = data->inputData->size();
    for (size_t k = 0; k < N; ++k) {
        float freq = (float)k * data->sampleRate / N;
        float Hf = (freq * freq) / (freq * freq + data->deltaF * data->deltaF);
        (*data->inputData)[k] *= Hf;
    }
    data->outputData = *(data->inputData);
    return nullptr;
}

// FIR Filter
void* applyFIRFilterThread(void* args) {
    FilterArgs* data = (FilterArgs*)args;
    size_t M = data->coefficientsH.size();
    size_t N = data->inputData->size();
    vector<float> output(N, 0);
    for (size_t n = 0; n < N; ++n) {
        for (size_t k = 0; k < M; ++k) {
            if (n >= k) output[n] += data->coefficientsH[k] * (*data->inputData)[n - k];
        }
    }
    data->outputData = output;
    return nullptr;
}

// IIR Filter
void* applyIIRFilterThread(void* args) {
    FilterArgs* data = (FilterArgs*)args;
    size_t M = data->coefficientsB.size();
    size_t N = data->coefficientsA.size();
    size_t L = data->inputData->size();
    vector<float> output(L, 0);
    for (size_t n = 0; n < L; ++n) {
        for (size_t k = 0; k < M; ++k) {
            if (n >= k) output[n] += data->coefficientsB[k] * (*data->inputData)[n - k];
        }
        for (size_t j = 1; j < N; ++j) {
            if (n >= j) output[n] -= data->coefficientsA[j] * output[n - j];
        }
    }
    data->outputData = output;
    return nullptr;
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

void writeWavFile(const std::string& outputFile, const std::vector<float>& data,  SF_INFO& fileInfo) {
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
    if (argc != 2) {
        printf("Invalid arguments\n");
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = "filtered_output.wav";

    SF_INFO fileInfo;
    vector<float> audioData, a, b, h;
    memset(&fileInfo, 0, sizeof(fileInfo));

    generate_random_floats(a, 2000, -1.0f, 1.0f);
    generate_random_floats(b, 2000, -1.0f, 1.0f);
    generate_random_floats(h, 2000, -1.0f, 1.0f);

    readWavFile(inputFile, audioData, fileInfo);

    // Define thread data
    pthread_t notchThread, bandThread, FIRThread, IIRThread;
    FilterArgs notchArgs, bandArgs, FIRArgs, IIRArgs;

    float deltaF = 30000.0f;
    float sampleRate = fileInfo.samplerate;
    float f0 = 945;
    int n = 3;

    // Initialize arguments for each thread
    notchArgs = {&audioData, {}, {}, {}, {}, 0, sampleRate, f0, n};
    bandArgs = {&audioData, {}, {}, {}, {}, deltaF, sampleRate, 0, 0};
    FIRArgs = {&audioData, {}, {}, {}, h, 0, 0, 0, 0};
    IIRArgs = {&audioData, {}, a, b, {}, 0, 0, 0, 0};

    auto start = high_resolution_clock::now();
    auto start_notch = high_resolution_clock::now();
    auto start_band = high_resolution_clock::now();
    auto start_IIR = high_resolution_clock::now();
    auto start_FIR = high_resolution_clock::now();

    pthread_create(&notchThread, nullptr, applyNotchFilterThread, &notchArgs);
    pthread_create(&bandThread, nullptr, applyBandPassFilterThread, &bandArgs);
    pthread_create(&FIRThread, nullptr, applyFIRFilterThread, &FIRArgs);
    pthread_create(&IIRThread, nullptr, applyIIRFilterThread, &IIRArgs);

    // Wait for threads to finish
    pthread_join(notchThread, nullptr);
    auto end_notch = high_resolution_clock::now();
    pthread_join(bandThread, nullptr);
    auto end_band = high_resolution_clock::now();
    pthread_join(FIRThread, nullptr);
    auto end_IIR = high_resolution_clock::now();
    pthread_join(IIRThread, nullptr);
    auto end_FIR = high_resolution_clock::now();

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start);
    auto notch_duration = duration_cast<milliseconds>(end_notch - start_notch);
    auto band_duration = duration_cast<milliseconds>(end_band - start_band);
    auto FIR_duration = duration_cast<milliseconds>(end_FIR - start_FIR);
    auto IIR_duration = duration_cast<milliseconds>(end_IIR - start_IIR);
    writeWavFile(outputFile, bandArgs.outputData, fileInfo);

    cout << "Band-pass Filter: " << band_duration.count() << " ms" << endl;
    cout << "Notch Filter: " << notch_duration.count() << " ms" << endl;
    cout << "IRR Filter: " << IIR_duration.count() << " ms" << endl;
    cout << "FIR Filter: " << FIR_duration.count() << " ms" << endl;
    cout << "All filters applied in parallel: " << duration.count() << " ms" << endl;

    return 0;
}