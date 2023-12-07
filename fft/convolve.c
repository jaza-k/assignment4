/**
JAZA KHAN (UCID 30119100)

Sme boilerplate code taken from Tutorial T03/T04
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

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
// void convolve(float x[], int N, float h[], int M, float y[], int P);
void four1(double data[], int nn, int isign);
int calculateNextPowerOfTwo(int size);
double* padData(float* data, int originalSize, int paddedSize);
void complexMultiply(double* data1, double* data2, int size);
void trimAndNormalize(double* data, int paddedSize, float* outputData, int outputSize);

void trimAndNormalize(double* data, int paddedSize, float* outputData, int outputSize) {
    // find maximum absolute value for normalization
    double maxVal = 0.0;
    for (int i = 0; i < paddedSize * 2; i += 2) {
        if (fabs(data[i]) > maxVal) {
            maxVal = fabs(data[i]);
        }
    }
    // avoid division by 0
    if (maxVal == 0) {
        maxVal = 1;
    }

    for (int i = 0; i < outputSize; ++i) {
        if (i < paddedSize) {
            outputData[i] = (float)(data[2 * i] / maxVal);
        } else {
            // if output size is larger than padded size, pad with 0s
            outputData[i] = 0.0f;
        }
    }
}

// taken from test.c //
void four1(double data[], int nn, int isign) {
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data[j], data[i]);
            SWAP(data[j+1], data[i+1]);
        }
        m = nn;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    mmax = 2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j+1];
                tempi = wr * data[j+1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}

int calculateNextPowerOfTwo(int size) {
    int power = 1;
    while (power < size) {
        power *= 2;
    }
    return power;
}

void complexMultiply(double* data1, double* data2, int size) {
    for (int i = 0; i < size; i += 2) {
        double a = data1[i]; // real, first array
        double b = data1[i + 1]; // imaginary, first array
        double c = data2[i];     // real, second array
        double d = data2[i + 1]; // imaginary, second array

        // minimize the number of multiplication
        double ac_minus_bd = a * c - b * d;
        double ad_plus_bc = a * d + b * c;

        data1[i] = ac_minus_bd;
        data1[i + 1] = ad_plus_bc;
    }
}

double* padData(float* data, int originalSize, int paddedSize) {
    // allocate memory for padded data
    double* paddedData = (double*) calloc(paddedSize * 2, sizeof(double));

    if (!paddedData) {
        fprintf(stderr, "Memory allocation failed for padded data\n");
        exit(-1);
    }

    // optimize data layout for better cache performance
    for (int i = 0, j = 0; i < originalSize; ++i, j += 2) {
        paddedData[j] = (double) data[i]; // Real part
    }
    return paddedData;
}


// main line of execution
int main (int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s inputfile IRfile outputfile\n", argv[0]);
        exit(-1);
    }

    char *inputFile = argv[1];
    char *IRFile = argv[2];
    char *outputFile = argv[3];

    float *inputData, *IRData, *outputData;
    WavHeader inputHeader, IRHeader;
    int numInputSamples, numIRSamples, numOutputSamples, paddedSize;

    int subchunk2_size_input, subchunk2_size_IR;
    readWavFile(inputFile, &inputData, &inputHeader, &numInputSamples, &subchunk2_size_input);
    readWavFile(IRFile, &IRData, &IRHeader, &numIRSamples,  &subchunk2_size_IR);

    paddedSize = calculateNextPowerOfTwo(numInputSamples + numIRSamples - 1);

    double *paddedInputData = padData(inputData, numInputSamples, paddedSize);
    double *paddedIRData = padData(IRData, numIRSamples, paddedSize);

    // fft on both signals
    four1(paddedInputData - 1, paddedSize / 2, 1);
    four1(paddedIRData - 1, paddedSize / 2, 1);

    complexMultiply(paddedInputData, paddedIRData, paddedSize);

    four1(paddedInputData - 1, paddedSize / 2, -1);
    numOutputSamples = numInputSamples + numIRSamples - 1;
    outputData = (float *)calloc(numOutputSamples, sizeof(float));

    trimAndNormalize(paddedInputData, paddedSize, outputData, numOutputSamples);
    writeWavFile(outputFile, outputData, &inputHeader, numOutputSamples);

    free(inputData);
    free(IRData);
    free(outputData);
    free(paddedInputData);
    free(paddedIRData);

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
    if (fread(tempBuffer, sizeof(short), *num_samples, file) != *num_samples) {
        fprintf(stderr, "Failed to read audio data from %s\n", fileName);
        free(tempBuffer);
        fclose(file);
        exit(-1);
    }

    for (int i = 0; i < *num_samples; ++i) {
        (*audioData)[i] = tempBuffer[i] / 32768.0f;
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

    // update header sizes based on convolution output
    int dataBytes = num_samples * sizeof(short);
    int chunkSize = 36 + dataBytes;
    fwrite(header, sizeof(WavHeader), 1, file);

    fwrite("data", sizeof(char), 4, file);
    fwrite(&dataBytes, sizeof(int), 1, file);

    short *outputBuffer = (short *)malloc(dataBytes);
    
    // normalize and map back to short
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

inline float bytesToFloat(short s) {
    return s / 32768.0f;
}