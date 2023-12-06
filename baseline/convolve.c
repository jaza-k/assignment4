/**
JAZA KHAN (UCID 30119100)

Sme boilerplate code taken from Tutorial T03/T04
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

// WAV file header format
typedef struct {
    char chunk_id[4];
    int chunk_size;
    char format[4];
    char subchunk1_id[4];
    int subchunk1_size;
    short audio_format;
    short num_channels;
    int sample_rate;
    int byte_rate;
    short block_align;
    short bits_per_sample;
} WavHeader;

// Function definitions
float bytesToFloat(short s);
void readWavFile(char *fileName, float **audioData, WavHeader *header, int *num_samples, int *subchunk2_size);
void writeWavFile(char *fileName, float *audioData, WavHeader *header, int num_samples);
void convolve(float x[], int N, float h[], int M, float y[], int P);

// main line of execution
int main (int argc, char *argv[]){
   if (argc != 4) {
        fprintf(stderr, "Usage: %s inputfile IRfile outputfile\n", argv[0]);
        exit(-1);
    }

    char *inputFile = argv[1];
    char *IRFile = argv[2];
    char *outputFile = argv[3];

    float *inputData, *IRData, *outputData;
    WavHeader inputHeader, IRHeader;
    int numInputSamples, numIRSamples, numOutputSamples;

    int subchunk2_size_input, subchunk2_size_IR;
    readWavFile(inputFile, &inputData, &inputHeader, &numInputSamples, &subchunk2_size_input);
    readWavFile(IRFile, &IRData, &IRHeader, &numIRSamples, &subchunk2_size_IR);


    numOutputSamples = numInputSamples + numIRSamples - 1;
    outputData = (float *)calloc(numOutputSamples, sizeof(float));

    convolve(inputData, numInputSamples, IRData, numIRSamples, outputData, numOutputSamples);

    writeWavFile(outputFile, outputData, &inputHeader, numOutputSamples);

    // Free allocated memory
    free(inputData);
    free(IRData);
    free(outputData);

    return 0;
}

// Function to read WAV file
void readWavFile(char *fileName, float **audioData, WavHeader *header, int *num_samples, int *subchunk2_size) {
    FILE *file = fopen(fileName, "rb");
    if (!file) {
        fprintf(stderr, "Unable to open file %s\n", fileName);
        exit(-1);
    }

    fread(header, sizeof(WavHeader), 1, file);

    if (header->subchunk1_size != 16) {
        int remainder = header->subchunk1_size - 16;
        fseek(file, remainder, SEEK_CUR);
    }
    int subchunk2_size_value;
    char subchunk2_id[4];
    fread(subchunk2_id, sizeof(subchunk2_id), 1, file);
     fread(&subchunk2_size_value, sizeof(subchunk2_size_value), 1, file);
    *subchunk2_size = subchunk2_size_value;

    *num_samples = *subchunk2_size / (header->bits_per_sample / 8);
    *audioData = (float *)malloc(*num_samples * sizeof(float));

    short *tempBuffer = (short *)malloc(subchunk2_size_value);
    fread(tempBuffer, sizeof(short), *num_samples, file);

    for (int i = 0; i < *num_samples; ++i) {
        (*audioData)[i] = bytesToFloat(tempBuffer[i]);
    }

    free(tempBuffer);
    fclose(file);
}

// Function to write WAV file
void writeWavFile(char *fileName, float *audioData, WavHeader *header, int num_samples) {
    FILE *file = fopen(fileName, "wb");
    if (!file) {
        fprintf(stderr, "Unable to open file %s\n", fileName);
        exit(-1);
    }

    // Update header sizes based on convolution output
    int dataBytes = num_samples * sizeof(short);
    int chunkSize = 36 + dataBytes;
    fwrite(header, sizeof(WavHeader), 1, file);

    fwrite("data", sizeof(char), 4, file);
    fwrite(&dataBytes, sizeof(int), 1, file);

    short *outputBuffer = (short *)malloc(dataBytes);
    
    // Normalize and map back to short
    float maxValue = 0.0;
    for (int i = 0; i < num_samples; ++i) {
        if (fabs(audioData[i]) > maxValue) {
            maxValue = fabs(audioData[i]);
        }
    }

    for (int i = 0; i < num_samples; ++i) {
        outputBuffer[i] = (short)(audioData[i] / maxValue * 32767.0);
    }

    fwrite(outputBuffer, sizeof(short), num_samples, file);

    free(outputBuffer);
    fclose(file);
}

float bytesToFloat(short s) {
    return s / 32768.0;
}

/*
    The function convolve takes six arguments: 
        Two input arrays x[] and h[], their respective sizes N and M, and an output array y[] with size P.

    The first loop initializes the output array y[] to zero. 
        This is necessary because the convolution operation involves accumulating values in y[].

    The second loop (outer loop) iterates over each element of the input array x[].

    The third loop (inner loop) iterates over each element of the array h[]. 
        For each pair of elements x[n] and h[m], it adds their sum to the corresponding element in y[].
*/
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
    int n,m;

    /* Clear Output Buffer y[] */
    for (n=0; n < P; n++)
    {
        y[n] = 0.0;
    }

    /* Outer Loop: process each input value x[n] in turn */
    for (n=0; n<N; n++){
        /* Inner loop: process x[n] with each sample of h[n] */
        for (m=0; m<M; m++){
            y[n+m] += x[n] * h[m];
        }
    }
}