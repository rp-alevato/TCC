#include "aoa/estimator.h"
#include "misc/read_data_files.h"

#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

static constexpr double deg_pi = M_PI / 180;
static constexpr double pi_deg = 180 / M_PI;

int main() {
    std::vector<SamplesData> samples_data;
    std::vector<AoaAngles> music_results;
    read_files::get_iq_samples(samples_data, "data/iq_samples/close.txt");
    read_files::get_music_result_angles(music_results, "data/music_result_angles/close.txt");
    AoaEstimator estimator;
    AoaAngles angles;
    GradientSpecs gradient_specs = {1e-5, 1e-8, 0.1, 0.9};
    int sample_ind = 50 * 5;
    // for (std::size_t sample_ind = 0; sample_ind < 1000; sample_ind += 5) {
    angles = estimator.process_samples(samples_data[sample_ind], AoaTechnique::music,
                                       MusicSearch::coarse_grid, M_PI / 1800,
                                       MusicOptimization::gradient_adapt_lr, (6 * M_PI / 180), gradient_specs);
    std::cout << angles.azimuth * 180 / M_PI << ", " << angles.elevation * 180 / M_PI << "\n";
    std::cout << music_results[sample_ind].azimuth * 180 / M_PI << ", " << music_results[sample_ind].elevation * 180 / M_PI << "\n";
    //     if (estimator.iterations > 500) {
    //         std::cout << "sample_index: " << sample_ind << "\n";
    //     }
    // }
    return 0;
}