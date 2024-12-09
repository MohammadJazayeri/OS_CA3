#include <iostream>
#include <sndfile.h>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <complex>
using namespace std;
// Discrete Fourier Transform (DFT)
void DFT(const vector<float>& input, vector<complex<float>>& output) {
    size_t N = input.size();
    output.resize(N);
    printf("%ld\n", N);
    for (size_t k = 0; k < N; ++k) {
        complex<float> sum(0.0, 0.0);
        for (size_t n = 0; n < N; ++n) {
            float angle = -2.0f * M_PI * k * n / N;
            sum += input[n] * complex<float>(cos(angle), sin(angle));
        }
        output[k] = sum;
        // printf("IN\n");
    }
    printf("OUT\n");
}

// Inverse Discrete Fourier Transform (IDFT)
void IDFT(const vector<complex<float>>& input, vector<float>& output) {
    size_t N = input.size();
    output.resize(N);

    for (size_t n = 0; n < N; ++n) {
        complex<float> sum(0.0, 0.0);
        for (size_t k = 0; k < N; ++k) {
            float angle = 2.0f * M_PI * k * n / N;
            sum += input[k] * complex<float>(cos(angle), sin(angle));
        }
        output[n] = sum.real() / N;  // Normalize the result
    }
}

void FFT(vector<complex<float>>& data, bool inverse) {
    size_t N = data.size();
    if (N <= 1) return;

    // Split even and odd indices
    vector<complex<float>> even(N / 2), odd(N / 2);
    for (size_t i = 0; i < N / 2; ++i) {
        even[i] = data[i * 2];
        odd[i] = data[i * 2 + 1];
    }

    // Recursively apply FFT
    FFT(even, inverse);
    FFT(odd, inverse);

    // Combine results
    float angle = (inverse ? 2.0f : -2.0f) * M_PI / N;
    complex<float> w(1.0, 0.0), wn(cos(angle), sin(angle));
    for (size_t i = 0; i < N / 2; ++i) {
        data[i] = even[i] + w * odd[i];
        data[i + N / 2] = even[i] - w * odd[i];
        if (inverse) {
            data[i] /= 2;
            data[i + N / 2] /= 2;
        }
        w *= wn;
    }
}


// Band-pass filter in the frequency domain
void applyBandPassFilter(vector<complex<float>>& freqData, float deltaF, float sampleRate) {
    size_t N = freqData.size();

    for (size_t k = 0; k < N; ++k) {
        float freq = (float)k * sampleRate / N;  // Frequency bin
        float Hf = (freq * freq) / (freq * freq + deltaF * deltaF);  // Filter response

        // Apply the filter to both real and imaginary parts
        freqData[k] *= Hf;
    }
}

void applyNotchFilter(vector<complex<float>>& freqData, float f0, int n, float sampleRate) {
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

void writeWavFile(const string& outputFile, const vector<float>& data, SF_INFO& fileInfo) {
    // printf("%ld ------ %ld\n", fileInfo.frames);
    SNDFILE* outFile = sf_open(outputFile.c_str(), SFM_WRITE, &fileInfo);
    if (!outFile) {
        cerr << "Error opening output file: " << sf_strerror(NULL) << endl;
        exit(1);
    }
    // cout << "Data size: " << data.size() << ", Channels: " << fileInfo.channels 
    //       << ", Frames: " << fileInfo.frames << endl;

    sf_count_t numFrames = sf_writef_float(outFile, data.data(), fileInfo.frames);
    if (numFrames != fileInfo.frames) {
        cerr << "Error writing frames to file." << endl;
        sf_close(outFile);
        exit(1);
    }

    sf_close(outFile);
    cout << "Successfully wrote " << numFrames << " frames to " << outputFile << endl;
}


int main(int argc, char* argv[]) {
    if(argc != 2){
        printf("Invalid arguments\n");
    }
    string inputFile = "input.wav";
    string outputFile = "filtered_output.wav";

    SF_INFO fileInfo;
    vector<float> audioData;

    memset(&fileInfo, 0, sizeof(fileInfo));

    // Step 1: Read the WAV file
    readWavFile(inputFile, audioData, fileInfo);

    // Step 2: Apply DFT to the audio data
    vector<complex<float>> freqData(audioData.begin(), audioData.end());
    FFT(freqData, false);

    // Step 3: Apply Band-Pass Filter in the frequency domain
    float deltaF = 3000000.0f;  // Bandwidth parameter
    float sampleRate = fileInfo.samplerate;
    float f0 = 50;
    int n = 2;
    // applyNotchFilter(freqData, f0, n, sampleRate);
    applyBandPassFilter(freqData, deltaF, sampleRate);

    FFT(freqData, true);

    // Step 4: Apply IDFT to get the filtered signal
    vector<float> filteredData(freqData.size());
    for (size_t i = 0; i < freqData.size(); ++i) {
        filteredData[i] = freqData[i].real();
        if(i < 50)
            cout << filteredData[i] << endl;
    }

    
    // Step 5: Write the filtered audio data back to a new file
    writeWavFile(outputFile, filteredData, fileInfo);

    cout << "Band-Pass Filter applied successfully. Output written to " << outputFile << endl;
    return 0;
}













// int main() {
//     string inputFile = "input.wav";
//     string outputFile = "output.wav";

//     SF_INFO fileInfo;
//     vector<float> audioData;

//     memset(&fileInfo, 0, sizeof(fileInfo));

//     readWavFile(inputFile, audioData, fileInfo);

//     writeWavFile(outputFile, audioData, fileInfo);

//     return 0;
// }
