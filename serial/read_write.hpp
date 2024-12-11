#ifndef READ_WRITE_H
#define READ_WRITE_H

#include <vector>
#include <sndfile.h>
#include <string>

using namespace std;

void generate_random_floats(vector <float>& input, int size, float low, float high);
void readWavFile(const string& inputFile, vector<float>& data, SF_INFO& fileInfo);
void writeWavFile(const string& outputFile, const vector<float>& data,  SF_INFO fileInfo);
void print_info(int read_duration, int notch_duration, int band_duration, int FIR_duration, int IIR_duration);

#endif
