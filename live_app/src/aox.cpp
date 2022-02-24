/***************************************************************************************************
 * @file   aox.c
 * @brief  Module responsible for processing IQ samples and calculate angle estimation from them
 *         using the AoX library
 ***************************************************************************************************
 * # License
 * <b>Copyright 2019 Silicon Laboratories Inc. www.silabs.com</b>
 ***************************************************************************************************
 * The licensor of this software is Silicon Laboratories Inc. Your use of this software is governed
 * by the terms of Silicon Labs Master Software License Agreement (MSLA) available at
 * www.silabs.com/about-us/legal/master-software-license-agreement. This software is distributed to
 * you in Source Code format and is governed by the sections of the MSLA applicable to Source Code.
 **************************************************************************************************/

// Standard library headers
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// BGAPI libraries
#include "aox.h"
#include "bg_types.h"
#include "doa_estimator.h"

extern "C" {
#include "sl_rtl_clib_api.h"
}

/***************************************************************************************************
 * Public Variable Declarations
 **************************************************************************************************/

iqSamples_t iqSamples, iqSamplesBuffer;
SamplesData samples_data;

MUTEX_T iqSamplesCriticalSection;
MUTEX_T bgBufferCriticalSection;
MUTEX_T printfCriticalSection;
#ifndef WINDOWS
MUTEX_T newSamplesAvailableMutex;
#endif

// Conditional variable
CONDITION_T newSamplesAvailable;

// Application state
eAOX_APP_CTRL eAppCtrl;

/***************************************************************************************************
 * Static Variable Declarations
 **************************************************************************************************/

// State variables
static uint8_t aox_initialized = 0;
static sl_rtl_aox_libitem libitem;
static sl_rtl_util_libitem util_libitem;

/***************************************************************************************************
 * Static Function Declarations
 **************************************************************************************************/

static void aoxWaitNewSamples();
static enum sl_rtl_error_code aoxProcessSamples(iqSamples_t* samples, float* azimuth, float* elevation, uint32_t* qaResult);
static float calcFrequencyFromChannel(uint8_t channel);
static uint32_t allocate2DFloatBuffer(float*** buf, int rows, int cols);

/***************************************************************************************************
 * Public Function Definitions
 **************************************************************************************************/

void aoxInit(void* args) {
    if (!aox_initialized) {
        // Make sure IQ sample buffers are not read/written while buffers are not initialized
        ENTER_MUTEX(&iqSamplesCriticalSection);

        printf("AoX library init...\r\n");

        // Initialize AoX library
        sl_rtl_aox_init(&libitem);
        // Set the number of snapshots - how many times the antennas are scanned during one measurement
        sl_rtl_aox_set_num_snapshots(&libitem, numSnapshots);
        // Set the antenna array type
        sl_rtl_aox_set_array_type(&libitem, AOX_ARRAY_TYPE);
        // Select mode (high speed/high accuracy/etc.)
        sl_rtl_aox_set_mode(&libitem, AOX_MODE);
        // Enable IQ sample quality analysis processing
        sl_rtl_aox_iq_sample_qa_configure(&libitem);
        // Create AoX estimator
        sl_rtl_aox_create_estimator(&libitem);

        // Initialize an util item
        sl_rtl_util_init(&util_libitem);
        sl_rtl_util_set_parameter(&util_libitem, SL_RTL_UTIL_PARAMETER_AMOUNT_OF_FILTERING, 0.6f);

        // // Allocate buffers for IQ samples
        allocate2DFloatBuffer(&(iqSamples.i_samples), numSnapshots, numArrayElements);
        allocate2DFloatBuffer(&(iqSamples.q_samples), numSnapshots, numArrayElements);
        allocate2DFloatBuffer(&(iqSamples.ref_i_samples), 1, ref_period_samples);
        allocate2DFloatBuffer(&(iqSamples.ref_q_samples), 1, ref_period_samples);
        allocate2DFloatBuffer(&(iqSamplesBuffer.i_samples), numSnapshots, numArrayElements);
        allocate2DFloatBuffer(&(iqSamplesBuffer.q_samples), numSnapshots, numArrayElements);
        allocate2DFloatBuffer(&(iqSamplesBuffer.ref_i_samples), 1, ref_period_samples);
        allocate2DFloatBuffer(&(iqSamplesBuffer.ref_q_samples), 1, ref_period_samples);

        aox_initialized = 1;

        // Now buffers can be read/written
        EXIT_MUTEX(&iqSamplesCriticalSection);
    }
}

THREAD_RETURN_T aoxMain(void* args) {
    double azimuth;
    double elevation;
    uint32_t qualityResult;

    // If not initialized, initialize
    if (!aox_initialized) {
        aoxInit(args);
    }

    // Define and initialize custom DoA estimator
    DoaEstimator estimator;
    DoaAngles angles;

    while (eAppCtrl != eAOX_SHUTDOWN) {
        // Wait for new IQ samples
        aoxWaitNewSamples();

        // Transform IQ samples to a Eigen matrix

        // for (int i = 0; i < ref_period_samples; ++i) {
        //     samples_reference(i) = Eigen::dcomplex(iqSamples.ref_i_samples[0][i], iqSamples.ref_q_samples[0][i]);
        // }
        // for (int i = 0; i < numArrayElements; ++i) {
        //     for (int j = 0; j < numSnapshots; ++j) {
        //         samples(i, j) = Eigen::dcomplex(iqSamples.i_samples[j][i], iqSamples.q_samples[j][i]);
        //     }
        // }
    }

    THREAD_EXIT;
}

/***************************************************************************************************
 * Custom Static Function Definitions
 **************************************************************************************************/

/***************************************************************************************************
 * Static Function Definitions
 **************************************************************************************************/

static void aoxWaitNewSamples() {
    static int n_samples_saved = 0;
    static constexpr int first_samples_discarded = 10;
    static constexpr int max_samples_saved = (10000 + first_samples_discarded);

    // Make sure IQ sample buffer is not read while BG thread is writing it
    ENTER_MUTEX(&iqSamplesCriticalSection);
    // Wait until new samples are available (signalled by BG thread)
    CONDITION_WAIT(&newSamplesAvailable, &iqSamplesCriticalSection);

    iqSamples.channel = iqSamplesBuffer.channel;
    iqSamples.rssi = iqSamplesBuffer.rssi;

    // Open file to append IQs

    samples_data.channel_frequency = calcFrequencyFromChannel(iqSamplesBuffer.channel);

    // Copy reference IQ samples from buffer into working copy
    for (int i = 0; i < ref_period_samples; ++i) {
        iqSamples.ref_i_samples[0][i] = iqSamplesBuffer.ref_i_samples[0][i];
        iqSamples.ref_q_samples[0][i] = iqSamplesBuffer.ref_q_samples[0][i];
        samples_data.samples_reference(i) = Eigen::dcomplex(iqSamples.ref_i_samples[0][i], iqSamples.ref_q_samples[0][i]);
    }

    // Copy IQ samples from buffer into working copy
    for (int i = 0; i < numArrayElements; ++i) {
        for (int j = 0; j < numSnapshots; ++j) {
            iqSamples.i_samples[j][i] = iqSamplesBuffer.i_samples[j][i];
            iqSamples.q_samples[j][i] = iqSamplesBuffer.q_samples[j][i];
            samples_data.samples(i, j) = Eigen::dcomplex(iqSamples.i_samples[j][i], iqSamples.q_samples[j][i]);
        }
    }

    // // Process new IQ Samples using the custom DoA estimator
    // static DoaEstimator estimator;
    // DoaAngles angles;
    // estimator.load_samples(samples_data);
    // angles = estimator.process_samples(DoaTechnique::music, MusicSearch::linear_grid_gradient, M_PI / 10);
    // printf("Azimuth: %6.1f  \tElevation: %6.1f  \trssi: %6.0f \tch: %2d\n",
    //        (angles.azimuth * 180 / M_PI), (angles.elevation * 180 / M_PI), iqSamples.rssi / 1.0, iqSamples.channel);

    // Process new IQ samples and calculate Angle of Arrival(azimuth, elevation)
    float phase_rotation, azimuth_float, elevation_float;
    uint32_t qualityResult;
    enum sl_rtl_error_code ret = aoxProcessSamples(&iqSamples, &azimuth_float, &elevation_float, &qualityResult);
    // sl_rtl_aox_process will return SL_RTL_ERROR_ESTIMATION_IN_PROGRESS until it has received enough packets for angle estimation
    if (ret == SL_RTL_ERROR_SUCCESS) {
        float distance;
        sl_rtl_aox_calculate_iq_sample_phase_rotation(&libitem, 2.0f, iqSamples.ref_i_samples[0], iqSamples.ref_q_samples[0], ref_period_samples, &phase_rotation);
        // Calculate distance from RSSI, and calculate a rough position estimation
        sl_rtl_util_rssi2distance(TAG_TX_POWER, iqSamples.rssi / 1.0, &distance);
        sl_rtl_util_filter(&util_libitem, distance, &distance);
    }

    if (qualityResult == 0 && n_samples_saved < max_samples_saved) {
        // Print out Silabs results
        printf("azimuth: %6.1f  \televation: %6.1f  \trssi: %6.0f \tch: %2d\n",
               azimuth_float, elevation_float, iqSamples.rssi / 1.0, iqSamples.channel);
        // Save samples data
        n_samples_saved++;
        if (n_samples_saved > first_samples_discarded) {
            std::ofstream iq_file;
            iq_file.open("iq_samples.txt", std::ios_base::app);

            iq_file << "channel frequency: " << samples_data.channel_frequency << "\n";
            iq_file << "rssi: " << iqSamples.rssi << "\n";
            iq_file << "phase rotation: " << phase_rotation << "\n";
            iq_file << "azimuth, elevation: (" << azimuth_float << ", " << elevation_float << ")\n";

            iq_file << "reference samples:\n";
            for (int i = 0; i < ref_period_samples; ++i) {
                iq_file << "(" << iqSamples.ref_i_samples[0][i] << "," << iqSamples.ref_q_samples[0][i] << ")\n";
            }

            iq_file << "samples:\n";
            for (int i = 0; i < numArrayElements; ++i) {
                for (int j = 0; j < numSnapshots; ++j) {
                    iq_file << "(" << iqSamples.i_samples[j][i] << "," << iqSamples.q_samples[j][i] << ")";
                    if (j < (numSnapshots - 1)) {
                        iq_file << ", ";
                    } else {
                        iq_file << "\n";
                    }
                }
            }
            iq_file << "\n";
        }
    }

    // Now IQ sample buffer can be written by BG thread again
    if (n_samples_saved >= max_samples_saved) {
        std::cout << "Collected " << (max_samples_saved - first_samples_discarded) << " samples. No more samples to collect.\n";
    }
    EXIT_MUTEX(&iqSamplesCriticalSection);
}

static enum sl_rtl_error_code aoxProcessSamples(iqSamples_t* samples, float* azimuth, float* elevation, uint32_t* qaResult) {
    float phase_rotation;

    // Calculate phase rotation from reference IQ samples
    sl_rtl_aox_calculate_iq_sample_phase_rotation(&libitem, 2.0f, samples->ref_i_samples[0], samples->ref_q_samples[0], ref_period_samples, &phase_rotation);

    // Provide calculated phase rotation to the estimator
    sl_rtl_aox_set_iq_sample_phase_rotation(&libitem, phase_rotation);

    // Estimate Angle of Arrival / Angle of Departure from IQ samples
    enum sl_rtl_error_code ret = sl_rtl_aox_process(&libitem, samples->i_samples, samples->q_samples, calcFrequencyFromChannel(samples->channel), azimuth, elevation);

    // fetch the quality results
    *qaResult = sl_rtl_aox_iq_sample_qa_get_results(&libitem);

    return ret;
}

static float calcFrequencyFromChannel(uint8_t channel) {
    static const uint8_t logical_to_physical_channel[40] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                                                            13, 14, 15, 16, 17, 18, 19, 20, 21,
                                                            22, 23, 24, 25, 26, 27, 28, 29, 30,
                                                            31, 32, 33, 34, 35, 36, 37, 38,
                                                            0, 12, 39};

    // return the center frequency of the given channel
    return 2402000000 + (2000000 * logical_to_physical_channel[channel]);
}

static uint32_t allocate2DFloatBuffer(float*** buf, int rows, int cols) {
    *buf = malloc(sizeof(float*) * rows);
    if (*buf == NULL) {
        return 0;
    }

    for (uint32_t i = 0; i < rows; i++) {
        (*buf)[i] = malloc(sizeof(float) * cols);
        if ((*buf)[i] == NULL) {
            return 0;
        }
    }

    return 1;
}
