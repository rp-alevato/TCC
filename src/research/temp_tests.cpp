#include "libdoa/doa_estimator.h"
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
    std::vector<DoaAngles> music_results;
    read_files::get_iq_samples(samples_data, "data/iq_samples/close.txt");
    read_files::get_music_result_angles(music_results, "data/music_result_angles/close.txt");
    DoaEstimator estimator;
    DoaAngles angles;
    GradientSpecs gradient_specs = {1e-5, 1e-6, 0, 0};
    int sample_ind = 50 * 5;
    angles = estimator.process_samples(samples_data[sample_ind], DoaTechnique::music,
                                       MusicSearch::coarse_grid, M_PI / 1800,
                                       MusicOptimization::finer_grid_search, (1 * M_PI / 180), gradient_specs);
    std::cout << angles.azimuth * 180 / M_PI << ", " << angles.elevation * 180 / M_PI << "\n";
    std::cout << music_results[sample_ind].azimuth * 180 / M_PI << ", " << music_results[sample_ind].elevation * 180 / M_PI << "\n";
    return 0;
}