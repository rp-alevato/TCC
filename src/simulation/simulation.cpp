#include "libdoa.h"

#include <Eigen/Dense>
#include <chrono>
#include <iostream>
#include <string>

static constexpr double deg_pi = M_PI / 180;
static constexpr double pi_deg = 180 / M_PI;

double runtime_tests(Eigen::Matrix<Eigen::dcomplex, 16, 4>& samples,
                     Eigen::Matrix<Eigen::dcomplex, 7, 1>& samples_reference,
                     double channel_frequency,
                     std::string technique_name,
                     DoaTechnique technique,
                     MusicSearchOptim search_optmization = MusicSearchOptim::simple_grid,
                     double grid_step = 2 * M_PI / 1000,
                     GradientOptimSpecs gradient_specs = {0.1, 1e-5, 1e-8},
                     int n_iterations = 10000) {

    DoaEstimator estimator;
    DoaAngles angles;
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n_iterations; i++) {
        estimator.load_samples(samples, samples_reference, channel_frequency);
        angles = estimator.process_samples(technique, search_optmization, grid_step, gradient_specs);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> runtime = t2 - t1;
    std::cout << "Technique: " << technique_name << "\n"
              << "Runtime: " << runtime.count() << "ms\n"
              << "azimuth: " << angles.azimuth * pi_deg
              << ", elevation: " << angles.elevation * pi_deg << "\n\n";
    return runtime.count();
}

int main() {
    Eigen::Matrix<Eigen::dcomplex, 16, 4> samples;
    Eigen::Matrix<Eigen::dcomplex, 7, 1> samples_reference;
    samples << Eigen::dcomplex(-0.047244, 0.551181), Eigen::dcomplex(0.551181, 0.086614), Eigen::dcomplex(0.110236, -0.559055), Eigen::dcomplex(-0.519685, -0.157480),
        Eigen::dcomplex(-0.740157, -0.614173), Eigen::dcomplex(-0.653543, 0.716535), Eigen::dcomplex(0.677165, 0.677165), Eigen::dcomplex(0.779528, -0.590551),
        Eigen::dcomplex(0.787402, 0.519685), Eigen::dcomplex(0.559055, -0.771654), Eigen::dcomplex(-0.732283, -0.637795), Eigen::dcomplex(-0.700787, 0.637795),
        Eigen::dcomplex(-0.031496, 0.094488), Eigen::dcomplex(0.094488, 0.047244), Eigen::dcomplex(0.070866, -0.086614), Eigen::dcomplex(-0.086614, -0.094488),
        Eigen::dcomplex(-0.700787, 0.244094), Eigen::dcomplex(0.173228, 0.708661), Eigen::dcomplex(0.700787, -0.157480), Eigen::dcomplex(-0.118110, -0.716535),
        Eigen::dcomplex(0.440945, -0.535433), Eigen::dcomplex(-0.511811, -0.472441), Eigen::dcomplex(-0.496063, 0.511811), Eigen::dcomplex(0.480315, 0.527559),
        Eigen::dcomplex(-0.141732, 0.685039), Eigen::dcomplex(0.685039, 0.165354), Eigen::dcomplex(0.188976, -0.685039), Eigen::dcomplex(-0.669291, -0.244094),
        Eigen::dcomplex(0.000000, -0.322835), Eigen::dcomplex(-0.322835, -0.007874), Eigen::dcomplex(-0.015748, 0.346457), Eigen::dcomplex(0.338583, 0.047244),
        Eigen::dcomplex(-0.433071, -0.267717), Eigen::dcomplex(-0.299213, 0.401575), Eigen::dcomplex(0.377953, 0.307087), Eigen::dcomplex(0.370079, -0.322835),
        Eigen::dcomplex(0.716535, 0.125984), Eigen::dcomplex(0.181102, -0.692913), Eigen::dcomplex(-0.685039, -0.244094), Eigen::dcomplex(-0.338583, 0.629921),
        Eigen::dcomplex(-0.535433, 0.023622), Eigen::dcomplex(-0.055118, 0.535433), Eigen::dcomplex(0.535433, 0.055118), Eigen::dcomplex(0.102362, -0.527559),
        Eigen::dcomplex(0.204724, -0.165354), Eigen::dcomplex(-0.141732, -0.228346), Eigen::dcomplex(-0.236220, 0.141732), Eigen::dcomplex(0.133858, 0.259843),
        Eigen::dcomplex(0.110236, -0.370079), Eigen::dcomplex(-0.354331, -0.165354), Eigen::dcomplex(-0.133858, 0.362205), Eigen::dcomplex(0.330709, 0.157480),
        Eigen::dcomplex(0.220472, 0.511811), Eigen::dcomplex(0.527559, -0.188976), Eigen::dcomplex(-0.212598, -0.527559), Eigen::dcomplex(-0.535433, 0.181102),
        Eigen::dcomplex(-0.283465, -0.456693), Eigen::dcomplex(-0.464567, 0.244094), Eigen::dcomplex(0.220472, 0.472441), Eigen::dcomplex(0.488189, -0.196850),
        Eigen::dcomplex(-0.007874, 0.055118), Eigen::dcomplex(0.062992, 0.007874), Eigen::dcomplex(0.000000, -0.062992), Eigen::dcomplex(-0.062992, 0.000000);
    samples_reference << std::polar<double>(1, 4.703218 * deg_pi),
        std::polar<double>(1, -174.413985 * deg_pi), std::polar<double>(1, 6.468812 * deg_pi),
        std::polar<double>(1, -172.648391 * deg_pi), std::polar<double>(1, 8.234406 * deg_pi),
        std::polar<double>(1, -170.882797 * deg_pi), std::polar<double>(1, 10.0 * deg_pi);

    runtime_tests(samples, samples_reference, 2444000000.0, "ESPRIT", DoaTechnique::esprit);
    runtime_tests(samples, samples_reference, 2444000000.0, "MUSIC Coarse Gradient", DoaTechnique::music, MusicSearchOptim::coarse_grid_gradient, 2 * M_PI / 20);
    runtime_tests(samples, samples_reference, 2444000000.0, "MUSIC Linear Gradient", DoaTechnique::music, MusicSearchOptim::linear_grid_gradient, 2 * M_PI / 20);
}