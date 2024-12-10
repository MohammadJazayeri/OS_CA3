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

// Band-pass filter in the frequency domain
void applyBandPassFilter(vector<float>& freqData, float deltaF, float sampleRate) {
    size_t N = freqData.size();

    for (size_t k = 0; k < N; ++k) {
        float freq = (float)k * sampleRate / N;
        float Hf = (freq * freq) / (freq * freq + deltaF * deltaF);

        freqData[k] *= Hf;
    }
}

void applyNotchFilter(vector<float>& freqData, float f0, int n, float sampleRate) {
    size_t N = freqData.size();
    for (size_t k = 0; k < N; ++k) {
        float f = (k < N / 2 ? k : k - N) * sampleRate / N; // Handle negative frequencies
        float Hf = 1.0f / (pow(f / f0, 2 * n) + 1.0f);     // Notch filter response
        freqData[k] *= Hf;                                // Apply filter
    }
}

vector<float> applyFIRFilter(const vector<float>& input, const vector<float>& h) {
    size_t M = h.size();          // Number of filter coefficients
    size_t N = input.size();      // Length of input signal
    std::vector<float> output(N); // Output signal

    for (size_t n = 0; n < N; ++n) {
        output[n] = 0;
        for (size_t k = 0; k < M; ++k) {
            if (n >= k) { // Avoid accessing out-of-bounds indices
                output[n] += h[k] * input[n - k];
            }
        }
    }
    return output;
}

vector<float> applyIIRFilter(const vector<float>& input, const vector<float>& a, const vector<float>& b) {
    size_t M = b.size();          // Number of feedforward coefficients
    size_t N = a.size();          // Number of feedback coefficients
    size_t L = input.size();      // Length of input signal
    std::vector<float> output(L); // Output signal

    for (size_t n = 0; n < L; ++n) {
        output[n] = 0;

        // Feedforward contribution
        for (size_t k = 0; k < M; ++k) {
            if (n >= k) { // Avoid accessing out-of-bounds indices
                output[n] += b[k] * input[n - k];
            }
        }

        // Feedback contribution
        for (size_t j = 1; j < N; ++j) {
            if (n >= j) { // Avoid accessing out-of-bounds indices
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
    if(argc != 2){
        printf("Invalid arguments\n");
    }
    string inputFile = argv[1];
    string outputFile = "filtered_output.wav";

    SF_INFO fileInfo;
    vector<float> audioData, a, b, h;
    memset(&fileInfo, 0, sizeof(fileInfo));

    generate_random_floats(a, 100, -10.0f, 10.0f);
    generate_random_floats(b, 100, -5.0f, 5.0f);
    generate_random_floats(h, 100, -10.0f, 10.0f);

    // Step 1: Read the WAV file
    auto start = high_resolution_clock::now();
    readWavFile(inputFile, audioData, fileInfo);
    
    float deltaF = 30000.0f;
    float sampleRate = fileInfo.samplerate;
    float f0 = 50;
    int n = 2;
    // applyNotchFilter(audioData, f0, n, sampleRate);
    // applyBandPassFilter(audioData, deltaF, sampleRate);
    vector<float> filteredData(audioData.size());
    // filteredData = applyFIRFilter(audioData, h);
    filteredData = applyIIRFilter(audioData, a, b);


    // for (size_t i = 0; i < audioData.size(); ++i) {
    //     filteredData[i] = audioData[i];
    // }

    
    // Step 5: Write the filtered audio data back to a new file
    writeWavFile(outputFile, filteredData, fileInfo);

    cout << "Band-Pass Filter applied successfully. Output written to " << outputFile << endl;
    auto end = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(end - start);
    cout << "execution time: " << duration.count() << " ms" << endl;
    return 0;
}