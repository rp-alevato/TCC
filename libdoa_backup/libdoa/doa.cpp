#include "doa.h"

#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>

#define LIGHT_SPEED 299792458  // Speed of light in meters/second

/**
 * Direction of Arrival estimation library. MUSIC and ESPRIT algorithms for a 16-element Uniform Rectangular Array.
 *
 * @copyright   UFSC
 * @file        libaoa.cpp
 * @author      Pedro Lemos
 */

void aoa_estimator::initDoAEstimator(double elements_distance) {

    /* Distance between antenna array elements */
    this->distance = elements_distance;
}

double aoa_estimator::calculate_wavelength(double channel_frequency) {
    return LIGHT_SPEED / channel_frequency;
}

// The board gets 7 samples from the incoming signal to calculate carrier and CTE frequency.
// This allows to correct for phase differences (practical matter)
// reference_period_length is number of these samples (ideally 7)
double aoa_estimator::estimate_phase_rotation(double* i_samples, double* q_samples, int reference_period_length) {
    double sample_phase = 0;
    double next_sample_phase = 0;
    double phase_difference = 0;
    double phase_difference_acc = 0;
    for (int sample = 0; sample < (reference_period_length - 1); sample++) {
        sample_phase = std::atan2(q_samples[sample], i_samples[sample]);
        next_sample_phase = std::atan2(q_samples[sample + 1], i_samples[sample + 1]);
        phase_difference = std::atan2(std::sin(next_sample_phase - sample_phase), std::cos(next_sample_phase - sample_phase));
        phase_difference_acc += std::atan2(std::sin(phase_difference), std::cos(phase_difference));
    }
    phase_difference_acc /= (reference_period_length - 1);
    phase_difference_acc *= 2 * 180 / M_PI;

    return phase_difference_acc;
}

void aoa_estimator::compensateRotation(double phase_rotation) {

    std::complex<double> phase_rot;
    phase_rot = {0, M_PI * phase_rotation / 180};

    for (int antenna = 0; antenna < NUM_ANTENNAS; antenna++) {
        for (int snapshot = 0; snapshot < IQLENGTH; snapshot++) {
            this->x(antenna, snapshot) = (this->x(antenna, snapshot)) * (std::exp((std::complex<double>(antenna, 0)) * phase_rot));
        }
    }
}

// Transforms doubles into a complex
void aoa_estimator::load_x(double** i_samples, double** q_samples, int num_antennas, int num_samples) {

    // Auxiliary complex variable
    std::complex<double> iq_sample(0, 0);

    for (int n = 0; n < num_antennas; n++) {
        for (int k = 0; k < num_samples; k++) {
            // Transform to complex value
            iq_sample = {i_samples[k][n], q_samples[k][n]};
            // Assign to aoa_estimator class
            this->x(n, k) = iq_sample;
        }
    }
}

// Corresponds to equation 2.31 from Pedro's TCC. It's the autocorrelation matrix
void aoa_estimator::estimateRxx() {

    // Estimate covariance matrix
    this->Rxx = (this->x) * (this->x.adjoint()) / IQLENGTH;

    // Eigendecomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(this->Rxx);

    // Find the eigenvectors related to noise
    double maximum_value = 0;
    int index = 0;

    for (int n = 0; n < this->antenna_num; n++) {
        if (eigensolver.eigenvalues()(n) > maximum_value) {
            maximum_value = eigensolver.eigenvalues()(n);
            index = n;
        }
    }

    // Equation 2.33 (Pedro's TCC)
    this->Qn.block(0, 0, 16, index) = eigensolver.eigenvectors().block(0, 0, 16, index);
    this->Qn.block(0, index, 16, this->antenna_num - (index + 1)) = eigensolver.eigenvectors().block(0, index + 1, 16, this->antenna_num - (index + 1));

    // Assign Signal eigenvectors
    this->Qs = eigensolver.eigenvectors().block(0, index, 16, 1);

    this->Qn_Prd = Qn * Qn.adjoint();
}

void aoa_estimator::updateSteeringVector(double azimuth, double elevation, double channel_frequency) {

    this->x_st = (2 * M_PI * this->distance * std::sin(elevation) * std::cos(azimuth)) / calculate_wavelength(channel_frequency);
    this->y_st = (2 * M_PI * this->distance * std::sin(elevation) * std::sin(azimuth)) / calculate_wavelength(channel_frequency);

    this->phi_x = {0, this->x_st};
    this->phi_y = {0, this->y_st};

    for (int n = 0; n < 4; n++) {
        this->a_x(0, n) = std::exp((std::complex<double>(n, 0)) * this->phi_x);
        this->a_y(0, n) = std::exp((std::complex<double>(n, 0)) * this->phi_y);
    }

    for (int n = 0; n < 4; n++) {
        for (int k = 0; k < 4; k++) {
            this->st_vec(4 * n + k, 0) = this->a_y(0, n) * this->a_x(0, k);
        }
    }
}

double aoa_estimator::MUSICSpectrum(double azimuth, double elevation, double channel_frequency) {
    double product = 0;

    // Update the steering vector to be used right after
    updateSteeringVector(azimuth, elevation, channel_frequency);
    // Compute the MUSIC Spectrum denominator
    this->aux = this->st_vec.adjoint() * this->Qn_Prd * this->st_vec;

    // Just extract the real part and invert
    product = this->aux(0, 0).real();
    product = 1 / product;

    return product;
}

void aoa_estimator::processMUSIC(int azimuth_max, int elevation_max, double channel_frequency) {

    // Just some auxiliary variables
    double last_value = 0;
    double product;
    double azimuth, elevation;

    for (int n_azim = 0; n_azim < azimuth_max; n_azim++) {
        for (int n_elev = 0; n_elev < elevation_max; n_elev++) {
            azimuth = n_azim * (this->angle_step);
            elevation = n_elev * (this->angle_step);
            product = MUSICSpectrum(azimuth, elevation, channel_frequency);
            // Store the last max value and index
            if (product > last_value) {
                last_value = product;
                this->az_index = azimuth;
                this->el_index = elevation;
            }
        }
    }
}

void aoa_estimator::initSelectionMatrices() {

    this->Js_1y << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;

    this->Js_1x << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;

    this->Js_2y << 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;

    this->Js_2x << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
}

void aoa_estimator::processESPRIT(double channel_frequency) {

    // Equations 2.35 e 2.36 (Pedro's TCC)
    this->Qs_1x = (this->Js_1x) * (this->Qs);
    this->Qs_2x = (this->Js_2x) * (this->Qs);
    this->Qs_1y = (this->Js_1y) * (this->Qs);
    this->Qs_2y = (this->Js_2y) * (this->Qs);

    this->phi_X_LS = ((this->Qs_1x.adjoint() * this->Qs_1x).inverse()) * ((this->Qs_1x.adjoint()) * Qs_2x);
    this->phi_Y_LS = ((this->Qs_1y.adjoint() * this->Qs_1y).inverse()) * ((this->Qs_1y.adjoint()) * Qs_2y);

    double A, B, azimuth, elevation, wavelength, constant;
    A = std::arg(this->phi_Y_LS(0, 0));
    B = std::arg(this->phi_X_LS(0, 0));
    azimuth = std::atan2(A, B);
    wavelength = calculate_wavelength(channel_frequency);
    constant = (2 * M_PI * this->distance) / wavelength;
    elevation = std::asin(A / (constant * std::sin(azimuth)));

    this->az_esprit = 180 * azimuth / M_PI;
    this->el_esprit = 180 * elevation / M_PI;
}

void aoa_estimator::getProcessed(int algorithm, double* az, double* el) {
    if (algorithm == 1) {
        *az = this->az_esprit;
        *el = this->el_esprit;
    }
    if (algorithm == 2) {
        *az = 180 * this->az_index / M_PI;
        *el = 180 * this->el_index / M_PI;
    }
}
