// #include <iostream>
// #include <sndfile.h>
// #include <vector>
// #include <string>
// #include <cstring>
// #include <cmath>
// #include <pthread.h>
// #include <chrono>
// #include <random>

// using namespace std::chrono;
// using namespace std;

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

// // Function prototypes
// void* applyNotchFilterThread(void* args);
// void* applyBandPassFilterThread(void* args);
// void* applyFIRFilterThread(void* args);
// void* applyIIRFilterThread(void* args);

// // Shared data structure to pass arguments to threads
// struct FilterArgs {
//     vector<float>* inputData;
//     vector<float> outputData;
//     vector<float> coefficientsA;
//     vector<float> coefficientsB;
//     vector<float> coefficientsH;
//     float deltaF, sampleRate, f0;
//     int n;
// };

// // Notch Filter
// void* applyNotchFilterThread(void* args) {
//     FilterArgs* data = (FilterArgs*)args;
//     size_t N = data->inputData->size();
//     for (size_t k = 0; k < N; ++k) {
//         float f = (k < N / 2 ? k : k - N) * data->sampleRate / N;
//         float Hf = 1.0f / (pow(f / data->f0, 2 * data->n) + 1.0f);
//         (*data->inputData)[k] *= Hf;
//     }
//     data->outputData = *(data->inputData);
//     return nullptr;
// }

// // Band-pass Filter
// void* applyBandPassFilterThread(void* args) {
//     FilterArgs* data = (FilterArgs*)args;
//     size_t N = data->inputData->size();
//     for (size_t k = 0; k < N; ++k) {
//         float freq = (float)k * data->sampleRate / N;
//         float Hf = (freq * freq) / (freq * freq + data->deltaF * data->deltaF);
//         (*data->inputData)[k] *= Hf;
//     }
//     data->outputData = *(data->inputData);
//     return nullptr;
// }

// // FIR Filter
// void* applyFIRFilterThread(void* args) {
//     FilterArgs* data = (FilterArgs*)args;
//     size_t M = data->coefficientsH.size();
//     size_t N = data->inputData->size();
//     vector<float> output(N, 0);
//     for (size_t n = 0; n < N; ++n) {
//         for (size_t k = 0; k < M; ++k) {
//             if (n >= k) output[n] += data->coefficientsH[k] * (*data->inputData)[n - k];
//         }
//     }
//     data->outputData = output;
//     return nullptr;
// }

// // IIR Filter
// void* applyIIRFilterThread(void* args) {
//     FilterArgs* data = (FilterArgs*)args;
//     size_t M = data->coefficientsB.size();
//     size_t N = data->coefficientsA.size();
//     size_t L = data->inputData->size();
//     vector<float> output(L, 0);
//     for (size_t n = 0; n < L; ++n) {
//         for (size_t k = 0; k < M; ++k) {
//             if (n >= k) output[n] += data->coefficientsB[k] * (*data->inputData)[n - k];
//         }
//         for (size_t j = 1; j < N; ++j) {
//             if (n >= j) output[n] -= data->coefficientsA[j] * output[n - j];
//         }
//     }
//     data->outputData = output;
//     return nullptr;
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

// void writeWavFile(const std::string& outputFile, const std::vector<float>& data,  SF_INFO& fileInfo) {
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


// struct BandPassParams {
//     vector<float>* audioData; // Pointer to the audio data
//     float deltaF;
//     float sampleRate;
//     size_t start;
//     size_t end;
// };

// void* bandPassWorker(void* args) {
//     BandPassParams* params = (BandPassParams*)args;

//     vector<float>& audioData = *(params->audioData);
//     float deltaF = params->deltaF;
//     float sampleRate = params->sampleRate;

//     for (size_t k = params->start; k < params->end; ++k) {
//         float freq = (float)k * sampleRate / audioData.size();
//         float Hf = (freq * freq) / (freq * freq + deltaF * deltaF);
//         audioData[k] *= Hf;
//     }

//     return nullptr;
// }

// vector<float> applyBandPassFilter(vector<float>& audioData, float deltaF, float sampleRate, int numThreads) {
//     size_t N = audioData.size();
//     pthread_t threads[numThreads];
//     BandPassParams params[numThreads];

//     // Split the data into chunks for each thread
//     size_t chunkSize = N / numThreads;

//     for (int i = 0; i < numThreads; ++i) {
//         params[i].audioData = &audioData;
//         params[i].deltaF = deltaF;
//         params[i].sampleRate = sampleRate;
//         params[i].start = i * chunkSize;
//         params[i].end = (i == numThreads - 1) ? N : (i + 1) * chunkSize;

//         pthread_create(&threads[i], nullptr, bandPassWorker, &params[i]);
//     }

//     // Join threads
//     for (int i = 0; i < numThreads; ++i) {
//         pthread_join(threads[i], nullptr);
//     }

//     return audioData;
// }







// int main(int argc, char* argv[]) {
//     if (argc != 2) {
//         printf("Invalid arguments\n");
//         return 1;
//     }

//     string inputFile = argv[1];
//     string outputFile = "filtered_output.wav";

//     SF_INFO fileInfo;
//     vector<float> audioData, a, b, h;
//     memset(&fileInfo, 0, sizeof(fileInfo));

//     generate_random_floats(a, 2000, -1.0f, 1.0f);
//     generate_random_floats(b, 2000, -1.0f, 1.0f);
//     generate_random_floats(h, 2000, -1.0f, 1.0f);

//     readWavFile(inputFile, audioData, fileInfo);

//     // Define thread data
//     pthread_t notchThread, bandThread, FIRThread, IIRThread;
//     FilterArgs notchArgs, bandArgs, FIRArgs, IIRArgs;

//     float deltaF = 30000.0f;
//     float sampleRate = fileInfo.samplerate;
//     float f0 = 945;
//     int n = 3;

//     // Initialize arguments for each thread
//     notchArgs = {&audioData, {}, {}, {}, {}, 0, sampleRate, f0, n};
//     bandArgs = {&audioData, {}, {}, {}, {}, deltaF, sampleRate, 0, 0};
//     FIRArgs = {&audioData, {}, {}, {}, h, 0, 0, 0, 0};
//     IIRArgs = {&audioData, {}, a, b, {}, 0, 0, 0, 0};

//     auto start = high_resolution_clock::now();
//     auto start_notch = high_resolution_clock::now();
//     auto start_band = high_resolution_clock::now();
//     auto start_IIR = high_resolution_clock::now();
//     auto start_FIR = high_resolution_clock::now();

//     pthread_create(&notchThread, nullptr, applyNotchFilterThread, &notchArgs);
//     pthread_create(&bandThread, nullptr, applyBandPassFilterThread, &bandArgs);
//     pthread_create(&FIRThread, nullptr, applyFIRFilterThread, &FIRArgs);
//     pthread_create(&IIRThread, nullptr, applyIIRFilterThread, &IIRArgs);

//     // Wait for threads to finish
//     pthread_join(notchThread, nullptr);
//     auto end_notch = high_resolution_clock::now();
//     pthread_join(bandThread, nullptr);
//     auto end_band = high_resolution_clock::now();
//     pthread_join(FIRThread, nullptr);
//     auto end_IIR = high_resolution_clock::now();
//     pthread_join(IIRThread, nullptr);
//     auto end_FIR = high_resolution_clock::now();

//     auto end = high_resolution_clock::now();
//     auto duration = duration_cast<milliseconds>(end - start);
//     auto notch_duration = duration_cast<milliseconds>(end_notch - start_notch);
//     auto band_duration = duration_cast<milliseconds>(end_band - start_band);
//     auto FIR_duration = duration_cast<milliseconds>(end_FIR - start_FIR);
//     auto IIR_duration = duration_cast<milliseconds>(end_IIR - start_IIR);
//     writeWavFile(outputFile, bandArgs.outputData, fileInfo);

//     cout << "Band-pass Filter: " << band_duration.count() << " ms" << endl;
//     cout << "Notch Filter: " << notch_duration.count() << " ms" << endl;
//     cout << "IRR Filter: " << IIR_duration.count() << " ms" << endl;
//     cout << "FIR Filter: " << FIR_duration.count() << " ms" << endl;
//     cout << "All filters applied in parallel: " << duration.count() << " ms" << endl;

//     return 0;
// }

#include <iostream>
#include <sndfile.h>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <pthread.h>
#include <chrono>
#include <random>

using namespace std;
using namespace std::chrono;

// Struct for thread arguments
struct FilterChunkArgs {
    vector<float>* inputData;
    vector<float>* outputData;
    vector<float> a; // Coefficients for IIR filter
    vector<float> b; // Coefficients for IIR filter
    vector<float> h; // Coefficients for FIR filter
    float deltaF;
    float sampleRate;
    float f0;
    int n;
    size_t startIdx; // Start index of chunk
    size_t endIdx;   // End index of chunk
};

// Function to generate random floats
void generate_random_floats(vector<float>& input, int size, float low, float high) {
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<float> distribution(low, high);
    for (int i = 0; i < size; ++i) {
        input.push_back(distribution(generator));
    }
}

// Read WAV file
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
// Notch filter function for a chunk
void* applyNotchFilterChunk(void* args) {
    FilterChunkArgs* params = (FilterChunkArgs*)args;
    for (size_t k = params->startIdx; k < params->endIdx; ++k) {
        float f = (k < params->inputData->size() / 2 ? k : k - params->inputData->size()) * params->sampleRate / params->inputData->size();
        float Hf = 1.0f / (pow(f / params->f0, 2 * params->n) + 1.0f);
        (*params->outputData)[k] = (*params->inputData)[k] * Hf;
    }
    return nullptr;
}

// Band-pass filter function for a chunk
void* applyBandPassFilterChunk(void* args) {
    FilterChunkArgs* params = (FilterChunkArgs*)args;
    for (size_t k = params->startIdx; k < params->endIdx; ++k) {
        float freq = (float)k * params->sampleRate / params->inputData->size();
        float Hf = (freq * freq) / (freq * freq + params->deltaF * params->deltaF);
        (*params->outputData)[k] = (*params->inputData)[k] * Hf;
    }
    return nullptr;
}

// FIR filter function for a chunk
void* applyFIRFilterChunk(void* args) {
    FilterChunkArgs* params = (FilterChunkArgs*)args;
    size_t M = params->h.size();
    for (size_t n = params->startIdx; n < params->endIdx; ++n) {
        (*params->outputData)[n] = 0;
        for (size_t k = 0; k < M; ++k) {
            if (n >= k) {
                (*params->outputData)[n] += params->h[k] * (*params->inputData)[n - k];
            }
        }
    }
    return nullptr;
}

// IIR filter function for a chunk
void* applyIIRFilterChunk(void* args) {
    FilterChunkArgs* params = (FilterChunkArgs*)args;
    size_t M = params->b.size();
    size_t N = params->a.size();
    for (size_t n = params->startIdx; n < params->endIdx; ++n) {
        (*params->outputData)[n] = 0;
        for (size_t k = 0; k < M; ++k) {
            if (n >= k) {
                (*params->outputData)[n] += params->b[k] * (*params->inputData)[n - k];
            }
        }
        for (size_t j = 1; j < N; ++j) {
            if (n >= j) {
                (*params->outputData)[n] -= params->a[j] * (*params->outputData)[n - j];
            }
        }
    }
    return nullptr;
}

// Function to divide work into chunks
void divideAndRunFilter(void* (*filterFunction)(void*), vector<float>& inputData, vector<float>& outputData, FilterChunkArgs& baseArgs, int numThreads) {
    pthread_t threads[numThreads];
    FilterChunkArgs args[numThreads];
    size_t chunkSize = inputData.size() / numThreads;

    for (int i = 0; i < numThreads; ++i) {
        args[i] = baseArgs;
        args[i].startIdx = i * chunkSize;
        args[i].endIdx = (i == numThreads - 1) ? inputData.size() : (i + 1) * chunkSize;
        args[i].inputData = &inputData;
        args[i].outputData = &outputData;
        pthread_create(&threads[i], nullptr, filterFunction, &args[i]);
    }
    for (int i = 0; i < numThreads; ++i) {
        pthread_join(threads[i], nullptr);
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Invalid arguments\n");
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = "filtered_output_parallel.wav";

    SF_INFO fileInfo;
    vector<float> audioData, a, b, h, outputData_notch, outputData_band, outputData_FIR, outputData_IIR;
    memset(&fileInfo, 0, sizeof(fileInfo));

    generate_random_floats(a, 2000, -1.0f, 1.0f);
    generate_random_floats(b, 2000, -1.0f, 1.0f);
    generate_random_floats(h, 2000, -1.0f, 1.0f);

    readWavFile(inputFile, audioData, fileInfo);
    outputData_notch.resize(audioData.size());
    outputData_band.resize(audioData.size());
    outputData_FIR.resize(audioData.size());
    outputData_IIR.resize(audioData.size());

    FilterChunkArgs notch_args = {&audioData, &outputData_notch, a, b, h, 10000.0f, (float)fileInfo.samplerate, 945, 3};
    FilterChunkArgs band_args = {&audioData, &outputData_band, a, b, h, 10000.0f, (float)fileInfo.samplerate, 945, 3};
    FilterChunkArgs FIR_args = {&audioData, &outputData_FIR, a, b, h, 10000.0f, (float)fileInfo.samplerate, 945, 3};
    FilterChunkArgs IIR_args = {&audioData, &outputData_IIR, a, b, h, 10000.0f, (float)fileInfo.samplerate, 945, 3};

    int numThreads = 4;

    auto start = high_resolution_clock::now();
    divideAndRunFilter(applyNotchFilterChunk, audioData, outputData_notch, notch_args, numThreads);
    divideAndRunFilter(applyBandPassFilterChunk, audioData, outputData_band, band_args, numThreads);
    divideAndRunFilter(applyFIRFilterChunk, audioData, outputData_FIR, FIR_args, numThreads);
    divideAndRunFilter(applyIIRFilterChunk, audioData, outputData_IIR, IIR_args, numThreads);
    auto end = high_resolution_clock::now();

    writeWavFile(outputFile, outputData_band, fileInfo);

    auto duration = duration_cast<milliseconds>(end - start);
    cout << "All filters applied in parallel with chunking: " << duration.count() << " ms" << endl;

    return 0;
}
