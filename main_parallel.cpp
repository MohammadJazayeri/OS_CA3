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

void* applyNotchFilterChunk(void* args) {
    FilterChunkArgs* params = (FilterChunkArgs*)args;
    for (size_t k = params->startIdx; k < params->endIdx; ++k) {
        float Hf = 1.0f / (pow((*params->inputData)[k] / params->f0, 2 * params->n) + 1.0f);
        (*params->outputData)[k] = Hf;
    }
    return nullptr;
}

void* applyBandPassFilterChunk(void* args) {
    FilterChunkArgs* params = (FilterChunkArgs*)args;
    for (size_t k = params->startIdx; k < params->endIdx; ++k) {
        float Hf = ((*params->inputData)[k] * (*params->inputData)[k]) / ((*params->inputData)[k] * (*params->inputData)[k] + params->deltaF * params->deltaF);
        (*params->outputData)[k] = Hf;
    }
    return nullptr;
}

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
    string outputFile_band = "band_parallel.wav";
    string outputFile_notch = "notch_parallel.wav";
    string outputFile_IIR = "IIR_parallel.wav";
    string outputFile_FIR = "FIR_parallel.wav";

    SF_INFO fileInfo;
    vector<float> audioData, a, b, h, outputData_notch, outputData_band, outputData_FIR, outputData_IIR;
    memset(&fileInfo, 0, sizeof(fileInfo));

    generate_random_floats(a, 10000, 0.0f, 4.0f);
    generate_random_floats(b, 10000, 0.0f, 4.0f);
    generate_random_floats(h, 10000, 0.0f, 4.0f);

    auto start = high_resolution_clock::now();

    readWavFile(inputFile, audioData, fileInfo);

    auto end = high_resolution_clock::now();

    auto read_duration = duration_cast<milliseconds>(end - start);

    outputData_notch.resize(audioData.size());
    outputData_band.resize(audioData.size());
    outputData_FIR.resize(audioData.size());
    outputData_IIR.resize(audioData.size());

    float deltaF = 0.0f;
    float f0 = 0.00001f;
    int n = 2;

    // Filter arguments setup
    FilterChunkArgs notch_args = {&audioData, &outputData_notch, a, b, h, deltaF, f0, n};
    FilterChunkArgs band_args = {&audioData, &outputData_band, a, b, h, deltaF, f0, n};
    FilterChunkArgs FIR_args = {&audioData, &outputData_FIR, a, b, h, deltaF, f0, n};
    FilterChunkArgs IIR_args = {&audioData, &outputData_IIR, a, b, h, deltaF, f0, n};

    int numThreads = 5;

    // Measure Notch Filter
    auto start_notch = high_resolution_clock::now();
    divideAndRunFilter(applyNotchFilterChunk, audioData, outputData_notch, notch_args, numThreads);
    auto end_notch = high_resolution_clock::now();
    auto notch_duration = duration_cast<milliseconds>(end_notch - start_notch);

    // Measure Band-Pass Filter
    auto start_band = high_resolution_clock::now();
    divideAndRunFilter(applyBandPassFilterChunk, audioData, outputData_band, band_args, numThreads);
    auto end_band = high_resolution_clock::now();
    auto band_duration = duration_cast<milliseconds>(end_band - start_band);

    // Measure FIR Filter
    auto start_FIR = high_resolution_clock::now();
    divideAndRunFilter(applyFIRFilterChunk, audioData, outputData_FIR, FIR_args, numThreads);
    auto end_FIR = high_resolution_clock::now();
    auto FIR_duration = duration_cast<milliseconds>(end_FIR - start_FIR);

    // Measure IIR Filter
    auto start_IIR = high_resolution_clock::now();
    divideAndRunFilter(applyIIRFilterChunk, audioData, outputData_IIR, IIR_args, numThreads);
    auto end_IIR = high_resolution_clock::now();
    auto IIR_duration = duration_cast<milliseconds>(end_IIR - start_IIR);

    cout << "Read: " << read_duration.count() << " ms" << endl;
    cout << "Notch Filter Time: " << notch_duration.count() << " ms" << endl;
    cout << "Band-Pass Filter Time: " << band_duration.count() << " ms" << endl;
    cout << "FIR Filter Time: " << FIR_duration.count() << " ms" << endl;
    cout << "IIR Filter Time: " << IIR_duration.count() << " ms" << endl;

// Total execution time


    writeWavFile(outputFile_band, outputData_band, fileInfo);
    writeWavFile(outputFile_notch, outputData_notch, fileInfo);
    writeWavFile(outputFile_FIR, outputData_FIR, fileInfo);
    writeWavFile(outputFile_IIR, outputData_IIR, fileInfo);

    auto total_duration = notch_duration + band_duration + FIR_duration + IIR_duration;
    cout << "Total Execution Time: " << total_duration.count() << " ms" << endl;
    return 0;
}
