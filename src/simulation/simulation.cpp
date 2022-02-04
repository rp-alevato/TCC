/* Este exemplo mostra como a biblioteca é utilizada. Ex.: Como inicializar os estimadores, */
/* como carregar as amostras de IQ, como produzir uma estimativa de ângulo, etc. */

#include "../../libdoa_backup/libdoa/doa.h"

#include <chrono>
#include <iostream>
#include <stdio.h>  // fopen()

#define IQFILE "../../samples/iqsamples.txt"

// Allocate 2D double Buffer for IQ samples
int allocate2DFloatBuffer(double*** buf, int rows, int cols) {
    // TODO change to new
    *buf = (double**)malloc(sizeof(double*) * rows);
    if (*buf == NULL) {
        return 0;
    }

    for (int i = 0; i < rows; i++) {
        (*buf)[i] = (double*)malloc(sizeof(double) * cols);
        if ((*buf)[i] == NULL) {
            return 0;
        }
    }

    return 1;
}

int main(int argc, char const* argv[]) {

    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    // Load IQ samples from file
    double** i_samples;
    double** q_samples;

    // Allocate memory for the IQs
    allocate2DFloatBuffer(&i_samples, IQLENGTH, NUM_ANTENNAS);
    allocate2DFloatBuffer(&q_samples, IQLENGTH, NUM_ANTENNAS);

    // Open the file containing the IQ samples
    FILE* fp;

    fp = fopen(IQFILE, "r");
    if (fp == NULL) {
        printf("Error in opening file\n");
        return -1;
    }

    /** This section is used to read the samples as strings and the convert to double **/
    char iqsampleline[MAX_SAMPSIZE];
    // Two auxiliary variables because of strtok function
    char* i_samplech;
    char* q_samplech;
    char i_samplesch[MAX_SAMPSIZE];
    char q_samplesch[MAX_SAMPSIZE];

    int n = 0;
    int k = 0;

    while (fgets(iqsampleline, MAX_SAMPSIZE, fp) != NULL) {
        // Read one IQ sample as a strings
        i_samplech = strtok(iqsampleline, ",\n");
        q_samplech = strtok(NULL, ",\n");

        // Parse the strings into two strings containing I and Q
        strcpy(i_samplesch, i_samplech);
        strcpy(q_samplesch, q_samplech);

        // Convert these to double values
        i_samples[k][n] = atof(i_samplesch);
        q_samples[k][n] = atof(q_samplesch);

        n++;
        if (n == NUM_ANTENNAS) {
            n = 0;
            k++;
        }
    }

    fclose(fp);

    // Initialize the estimator with specific element distance
    auto doa_estimator = aoa_estimator();
    double elements_distance = 0.04;
    doa_estimator.initDoAEstimator(elements_distance);

    // Initialize selection matrices used in ESPRIT algorithm
    doa_estimator.initSelectionMatrices();

    double azimuth, elevation;

    auto t1 = high_resolution_clock::now();

    for (int i = 0; i < 100000; i++) {
        // Load iq samples into the estimator
        doa_estimator.load_x(i_samples, q_samples, NUM_ANTENNAS, IQLENGTH);
        // Compensate phase rotation from IQ samples
        double phase_rotation = -179.117203;
        doa_estimator.compensateRotation(phase_rotation);
        // Estimate covariance matrix
        doa_estimator.estimateRxx();
        doa_estimator.processESPRIT(2444000000);
        doa_estimator.getProcessed(1, &azimuth, &elevation);
    }

    // // MUSIC
    // doa_estimator.load_x(i_samples, q_samples, NUM_ANTENNAS, IQLENGTH);
    // doa_estimator.compensateRotation(-179.117203);
    // doa_estimator.estimateRxx();
    // doa_estimator.processMUSIC(360, 90, 2444000000);
    // doa_estimator.getProcessed(2, &azimuth, &elevation);
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> runtime = t2 - t1;
    std::cout << "Runtime: " << runtime.count() << "ms" << std::endl;

    printf("Processed value: az: %f, el: %f\n", azimuth, elevation);
    /* code */
    return 0;
}