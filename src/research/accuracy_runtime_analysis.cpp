#include "doa/estimator.h"
#include "misc/read_data_files.h"
#include "misc/statistics.h"
#include "misc/utility.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

const std::string iq_samples_dir = "data/iq_samples/";
const std::string music_results_dir = "data/music_result_angles/";
const std::string walk_filename = "office_walk.txt";
const std::string output_filename = "data/experimental_results/accuracy_runtime.csv";

static constexpr double fine_step = M_PI / 900;

void analysis(const std::vector<SamplesData>& samples_data, const std::vector<DoaAngles>& correct_results,
              const std::string technique_name, const double coarse_step, const GradientSpecs gradient_specs);
void save_csv_info_for_every_sample(std::ofstream& output_csv, const std::vector<SamplesData>& samples_data,
                                    const std::vector<DoaAngles>& correct_results, const std::string technique_name,
                                    const DoaTechnique doa_technique, const MusicSearch music_search,
                                    const MusicOptimization optimization, const double coarse_step,
                                    const GradientSpecs gradient_specs);

int main() {
    // Get Data
    std::vector<SamplesData> samples_data;
    std::vector<DoaAngles> correct_results;
    GradientSpecs gradient_specs;
    double coarse_step;
    std::string technique;

    std::ofstream output_csv;
    if (std::filesystem::exists(output_filename)) {
        throw std::runtime_error("File " + output_filename + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_filename, std::ios_base::app);
    output_csv << "technique,coarse_step,learning_rate,momentum,runtime,mae_len,rmse_len\n";
    output_csv.close();

    read_files::get_iq_samples(samples_data, (iq_samples_dir + walk_filename));
    read_files::get_music_result_angles(correct_results, (music_results_dir + walk_filename));

    technique = "esprit";
    coarse_step = 0;
    gradient_specs = {0, 0, 0, 0};
    analysis(samples_data, correct_results, technique, coarse_step, gradient_specs);

    technique = "music_coarse_fine_grid_search";
    coarse_step = 3;
    gradient_specs = {0, 0, 0, 0};
    analysis(samples_data, correct_results, technique, coarse_step, gradient_specs);

    technique = "music_gradient_simple";
    coarse_step = 5;
    gradient_specs = {1e-5, 1.5e-8, 0.1, 0};
    analysis(samples_data, correct_results, technique, coarse_step, gradient_specs);

    technique = "music_gradient_adapt_lr";
    coarse_step = 5;
    gradient_specs = {1e-5, 1.5e-8, 0.5, 0};
    analysis(samples_data, correct_results, technique, coarse_step, gradient_specs);

    technique = "music_gradient_momentum";
    coarse_step = 5;
    gradient_specs = {1e-5, 1.5e-8, 0.07, 0.85};
    analysis(samples_data, correct_results, technique, coarse_step, gradient_specs);

    technique = "music_gradient_nesterov";
    coarse_step = 5;
    gradient_specs = {1e-5, 1.5e-8, 0.05, 0.9};
    analysis(samples_data, correct_results, technique, coarse_step, gradient_specs);

    technique = "music_fine_grid_search";
    coarse_step = 0;
    gradient_specs = {0, 0, 0, 0};
    analysis(samples_data, correct_results, technique, coarse_step, gradient_specs);

    return 0;
}

void analysis(const std::vector<SamplesData>& samples_data, const std::vector<DoaAngles>& correct_results,
              const std::string technique_name, const double coarse_step, const GradientSpecs gradient_specs) {
    std::ofstream output_csv;
    std::vector<DoaAngles> results;
    std::vector<double> errors_len;
    DoaEstimator estimator;
    DoaTechnique doa_technique = DoaTechnique::esprit;
    MusicSearch music_search = MusicSearch::coarse_grid;
    MusicOptimization music_optimization = MusicOptimization::fine_grid_search;

    output_csv.open(output_filename, std::ios_base::app);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_filename);
    }

    if (technique_name == "esprit") {
        doa_technique = DoaTechnique::esprit;
    } else if (technique_name == "music_fine_grid_search") {
        doa_technique = DoaTechnique::music;
        music_search = MusicSearch::simple_grid;
    } else if (technique_name == "music_coarse_fine_grid_search") {
        doa_technique = DoaTechnique::music;
        music_search = MusicSearch::coarse_grid;
        music_optimization = MusicOptimization::fine_grid_search;
    } else if (technique_name == "music_gradient_simple") {
        doa_technique = DoaTechnique::music;
        music_search = MusicSearch::coarse_grid;
        music_optimization = MusicOptimization::gradient_simple;
    } else if (technique_name == "music_gradient_adapt_lr") {
        doa_technique = DoaTechnique::music;
        music_search = MusicSearch::coarse_grid;
        music_optimization = MusicOptimization::gradient_adapt_lr;
    } else if (technique_name == "music_gradient_momentum") {
        doa_technique = DoaTechnique::music;
        music_search = MusicSearch::coarse_grid;
        music_optimization = MusicOptimization::gradient_momentum;
    } else if (technique_name == "music_gradient_nesterov") {
        doa_technique = DoaTechnique::music;
        music_search = MusicSearch::coarse_grid;
        music_optimization = MusicOptimization::gradient_nesterov;
    } else {
        throw std::runtime_error("Wrong name for technique");
    }

    std::cout << "Calculating " << technique_name << "...\n";
    save_csv_info_for_every_sample(output_csv, samples_data, correct_results, technique_name, doa_technique,
                                   music_search, music_optimization, coarse_step, gradient_specs);
    output_csv.close();

    return;
}

void save_csv_info_for_every_sample(std::ofstream& output_csv, const std::vector<SamplesData>& samples_data,
                                    const std::vector<DoaAngles>& correct_results, const std::string technique_name,
                                    const DoaTechnique doa_technique, const MusicSearch music_search,
                                    const MusicOptimization optimization, const double coarse_step,
                                    const GradientSpecs gradient_specs) {
    auto double_precision = std::numeric_limits<double>::digits10;
    std::vector<double> errors;
    std::vector<DoaAngles> results;
    DoaEstimator estimator;
    double coarse_step_pi = utility::angle_to_pi(coarse_step);

    // Estimate angles
    auto t1 = std::chrono::high_resolution_clock::now();
    for (std::size_t sample_index = 0; sample_index < samples_data.size(); sample_index++) {
        results.push_back(estimator.process_samples(samples_data[sample_index], doa_technique,
                                                    music_search, fine_step,
                                                    optimization, coarse_step_pi, gradient_specs));
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float, std::milli> delta_t = t2 - t1;

    // Calculate errors
    for (std::size_t i = 0; i < results.size(); i++) {
        double result_az = utility::angle_to_degree(results[i].azimuth);
        double result_el = utility::angle_to_degree(results[i].elevation);
        double correct_result_az = utility::angle_to_degree(correct_results[i].azimuth);
        double correct_result_el = utility::angle_to_degree(correct_results[i].elevation);
        double error_az = correct_result_az - result_az;
        double error_el = correct_result_el - result_el;
        double error_len = std::hypot(error_az, error_el);
        errors.push_back(error_len);
    }

    // Calculate statistics
    double runtime = delta_t.count();
    double mae_len = stats::mae_double(errors);
    double rmse_len = stats::rmse_double(stats::mse_double(errors));
    if (technique_name == "esprit") {
        mae_len = 0;
        rmse_len = 0;
    }

    // Save values to CSV.
    // Columns: technique,coarse_step,learning_rate,momentum,runtime,mae_len,rmse_len
    output_csv << technique_name << ",";
    output_csv << coarse_step << ",";
    output_csv << gradient_specs.learning_rate << ",";
    output_csv << gradient_specs.momentum << ",";
    output_csv << std::setprecision(double_precision) << runtime << ",";
    output_csv << std::setprecision(double_precision) << mae_len << ",";
    output_csv << std::setprecision(double_precision) << rmse_len << "\n";

    return;
}
