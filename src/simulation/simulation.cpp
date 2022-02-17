#include "libdoa/libdoa.h"
#include "misc/read_data_files.h"

#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

static constexpr double deg_pi = M_PI / 180;
static constexpr double pi_deg = 180 / M_PI;

// double runtime_tests(SamplesData& samples_data,
//                      std::string technique_name,
//                      DoaTechnique technique,
//                      MusicSearch search_method = MusicSearch::simple_grid,
//                      double grid_step = 2 * M_PI / 1000,
//                      Gradient gradient = Gradient::momentum,
//                      GradientSpecs gradient_specs = {0.2, 1e-5, 1e-8, 0.3},
//                      int n_iterations = 1) {

//     DoaEstimator estimator;
//     DoaAngles angles;

//     // auto t1 = std::chrono::high_resolution_clock::now();
//     for (int i = 0; i < n_iterations; i++) {
//         estimator.load_samples(samples_data);
//         angles = estimator.process_samples(technique, search_method, grid_step, gradient, gradient_specs);
//     }
//     // auto t2 = std::chrono::high_resolution_clock::now();
//     // std::chrono::duration<double, std::milli> runtime = t2 - t1;
//     // // std::cout << "Technique: " << technique_name << "\n"
//     // //           << "Runtime: " << runtime.count() << "ms\n";
//     std::cout << "azimuth: " << angles.azimuth * pi_deg
//               << ", elevation: " << angles.elevation * pi_deg << "\n";
//     // std::cout << "SL: azimuth: " << samples_data.sl_azimuth
//     //           << ", elevation: " << samples_data.sl_elevation << "\n\n";
//     // return runtime.count();
//     return 0;
// }

void save_music_results(std::string iq_file_name, std::string music_results_file_name) {
    std::vector<SamplesData> samples_data_list;
    read_files::iq_samples(samples_data_list, "data/iq_samples/" + iq_file_name);

    DoaEstimator estimator;
    DoaAngles angles;
    std::ofstream music_results;
    music_results.open("data/music_results/" + music_results_file_name);
    auto double_precision = std::numeric_limits<long double>::digits10;
    std::cout << "Saving music results for " << music_results_file_name << "\n";
    for (unsigned int i = 0; i < samples_data_list.size(); i++) {
        estimator.load_samples(samples_data_list[i]);
        angles = estimator.process_samples(DoaTechnique::music, MusicSearch::coarse_grid, M_PI / 18000, M_PI / 360, Optimization::finer_grid_search);
        music_results << std::setprecision(double_precision) << "("
                      << angles.azimuth << "," << angles.elevation << ")\n";
        if ((i % 100) == 0) {
            std::cout << "Counter: " << i << "\n";
        }
    }
    music_results.close();
}

int main() {
    std::vector<SamplesData> samples_data_list;
    read_files::iq_samples(samples_data_list, "data/iq_samples/iq_samples_close.txt");
    DoaEstimator estimator;
    DoaAngles angles;
    for (int i = 1; i < 10; i++) {
        int sample = i * 300;
        std::cout << sample << "\n";
        estimator.load_samples(samples_data_list[sample]);
        angles = estimator.process_samples(DoaTechnique::music, MusicSearch::save_spectrum, M_PI / 360);
    }

    // save_music_results("iq_samples_close.txt", "close_18000.txt");
    // save_music_results("iq_samples_office_walk.txt", "office_walk.txt");
    // save_music_results("iq_samples_fixed_1.txt", "fixed_1.txt");
    // save_music_results("iq_samples_fixed_2.txt", "fixed_2.txt");
    // save_music_results("iq_samples_fixed_3.txt", "fixed_3.txt");
    // save_music_results("iq_samples_fixed_4.txt", "fixed_4.txt");
    // save_music_results("iq_samples_fixed_5.txt", "fixed_5.txt");
}