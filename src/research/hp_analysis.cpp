#include "libdoa/doa_estimator.h"
#include "misc/progress_bar.h"
#include "misc/read_data_files.h"
#include "misc/statistics.h"
#include "misc/utility.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Equation for best runtime coarse_step based on finer_step (desired precision):
// min(((n_coarse^2)/4) + (2*n_finer/n_coarse)^2

const std::string iq_samples_dir = "data/iq_samples/";
const std::string music_results_dir = "data/music_result_angles/";
const std::string close_filename = "close.txt";
const std::string walk_filename = "office_walk.txt";
const std::string output_dir = "data/experimental_results/hyperparameters/";

static constexpr std::size_t training_stride = 5;
static constexpr double finer_step = M_PI / 1800;

// Finer grid functions
void complete_finer_grid_analysis(const std::vector<double>& coarse_steps);
void finer_grid_analysis(const std::string output_filename, const std::vector<double>& coarse_steps,
                         const std::vector<SamplesData>& samples_data, const std::vector<DoaAngles>& correct_results_vector);
//  Gradient simple functions
void complete_gradient_simple_analysis(const std::vector<double>& coarse_steps,
                                       const std::vector<double>& learning_rates);
void gradient_simple_analysis(const std::string output_filename, const std::vector<double>& coarse_steps,
                              const std::vector<double>& learning_rates, const std::vector<SamplesData>& samples_data,
                              const std::vector<DoaAngles>& correct_results_vector);
//  Gradient momentum functions
void complete_gradient_momentum_analysis(const std::vector<double>& coarse_steps,
                                         const std::vector<double>& learning_rates,
                                         const std::vector<double>& momentums);
void gradient_momentum_analysis(const std::string output_filename, const std::vector<double>& coarse_steps,
                                const std::vector<double>& learning_rates, const std::vector<double>& momentums,
                                const std::vector<SamplesData>& samples_data,
                                const std::vector<DoaAngles>& correct_results_vector);
// Utility functions
void get_training_data(std::vector<SamplesData>& training_samples_close, std::vector<DoaAngles>& training_results_close,
                       std::vector<SamplesData>& training_samples_walk, std::vector<DoaAngles>& training_results_walk,
                       std::vector<SamplesData>& training_samples_both, std::vector<DoaAngles>& training_results_both);
void make_csv_columns(std::ofstream& output_csv, const std::string analysis_method);
void save_csv_info_for_every_sample(std::ofstream& output_csv, const std::string analysis_method,
                                    const std::vector<SamplesData>& samples_data,
                                    const std::vector<DoaAngles>& correct_results_vector,
                                    const MusicOptimization optimization, const double coarse_step,
                                    double learning_rate, double momentum);

int main() {
    // Create hyperparameters vectors
    std::vector<double> coarse_steps;
    for (std::size_t i = 1; i <= 10; i++) {
        coarse_steps.push_back(i);
    }
    std::vector<double> learning_rates;
    for (double i = 0.1; i <= 0.9; i += 0.1) {
        learning_rates.push_back(i);
    }
    std::vector<double> momentums;
    for (double i = 0.1; i <= 0.9; i += 0.1) {
        momentums.push_back(i);
    }

    complete_finer_grid_analysis(coarse_steps);

    complete_gradient_simple_analysis(coarse_steps, learning_rates);

    complete_gradient_momentum_analysis(coarse_steps, learning_rates, momentums);

    return 0;
}

// ****************************************************************************
// **************************  FINER GRID FUNCTIONS ***************************
// ****************************************************************************

void complete_finer_grid_analysis(const std::vector<double>& coarse_steps) {
    std::vector<SamplesData> training_samples_close, training_samples_walk, training_samples_both;
    std::vector<DoaAngles> training_results_close, training_results_walk, training_results_both;

    get_training_data(training_samples_close, training_results_close,
                      training_samples_walk, training_results_walk,
                      training_samples_both, training_results_both);

    std::cout << "Finer grid analysis: close.txt:\n";
    finer_grid_analysis("finer_grid_close.csv", coarse_steps, training_samples_close, training_results_close);
    std::cout << "Finer grid analysis: walk.txt:\n";
    finer_grid_analysis("finer_grid_walk.csv", coarse_steps, training_samples_walk, training_results_walk);
    std::cout << "Finer grid analysis: both:\n";
    finer_grid_analysis("finer_grid_both.csv", coarse_steps, training_samples_both, training_results_both);

    return;
}

void finer_grid_analysis(const std::string output_filename, const std::vector<double>& coarse_steps,
                         const std::vector<SamplesData>& samples_data, const std::vector<DoaAngles>& correct_results_vector) {
    std::ofstream output_csv;
    DoaEstimator doa_estimator;
    const std::string output_name = output_dir + output_filename;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_name);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    // Add columns to csv file
    make_csv_columns(output_csv, "finer_grid");

    std::cout << "Total number of steps: " << coarse_steps.size() << "\n";

    for (std::size_t coarse_index = 0; coarse_index < coarse_steps.size(); coarse_index++) {
        double coarse_step = coarse_steps[coarse_index];
        std::cout << "Coarse step: " << coarse_step << "\n";
        save_csv_info_for_every_sample(output_csv, "finer_grid", samples_data,
                                       correct_results_vector, MusicOptimization::finer_grid_search,
                                       coarse_step, 0, 0);
    }

    output_csv.close();
    std::cout << "\n";
    return;
}

// ****************************************************************************
// ************************  GRADIENT SIMPLE FUNCTIONS ************************
// ****************************************************************************

void complete_gradient_simple_analysis(const std::vector<double>& coarse_steps,
                                       const std::vector<double>& learning_rates) {
    std::vector<SamplesData> training_samples_close, training_samples_walk, training_samples_both;
    std::vector<DoaAngles> training_results_close, training_results_walk, training_results_both;

    get_training_data(training_samples_close, training_results_close,
                      training_samples_walk, training_results_walk,
                      training_samples_both, training_results_both);

    std::cout << "Gradient simple analysis: close.txt:\n";
    gradient_simple_analysis("gradient_simple_close.csv", coarse_steps, learning_rates,
                             training_samples_close, training_results_close);

    std::cout << "Gradient simple analysis: walk.txt:\n";
    gradient_simple_analysis("gradient_simple_walk.csv", coarse_steps, learning_rates,
                             training_samples_walk, training_results_walk);

    std::cout << "Gradient simple analysis: both:\n";
    gradient_simple_analysis("gradient_simple_both.csv", coarse_steps, learning_rates,
                             training_samples_both, training_results_both);

    return;
}

void gradient_simple_analysis(const std::string output_filename, const std::vector<double>& coarse_steps,
                              const std::vector<double>& learning_rates, const std::vector<SamplesData>& samples_data,
                              const std::vector<DoaAngles>& correct_results_vector) {
    std::ofstream output_csv;
    DoaEstimator doa_estimator;
    const std::string output_name = output_dir + output_filename;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_name);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    // Add columns to csv file
    make_csv_columns(output_csv, "gradient_simple");

    std::cout << "Total number of coarse_steps: " << coarse_steps.size() << "\n";
    std::cout << "Total number of learning_rates: " << learning_rates.size() << "\n";

    for (std::size_t coarse_index = 0; coarse_index < coarse_steps.size(); coarse_index++) {
        double coarse_step = coarse_steps[coarse_index];
        std::cout << "Coarse step: " << coarse_step << "\n";
        for (std::size_t lr_index = 0; lr_index < learning_rates.size(); lr_index++) {
            double learning_rate = learning_rates[lr_index];
            std::cout << "Learning rate: " << learning_rate << "\n";
            save_csv_info_for_every_sample(output_csv, "gradient_simple", samples_data,
                                           correct_results_vector, MusicOptimization::gradient_simple,
                                           coarse_step, learning_rate, 0);
        }
    }

    output_csv.close();
    std::cout << "\n";
    return;
}

// ****************************************************************************
// ***********************  GRADIENT MOMENTUM FUNCTIONS ***********************
// ****************************************************************************

void complete_gradient_momentum_analysis(const std::vector<double>& coarse_steps,
                                         const std::vector<double>& learning_rates,
                                         const std::vector<double>& momentums) {
    std::vector<SamplesData> training_samples_close, training_samples_walk, training_samples_both;
    std::vector<DoaAngles> training_results_close, training_results_walk, training_results_both;

    get_training_data(training_samples_close, training_results_close,
                      training_samples_walk, training_results_walk,
                      training_samples_both, training_results_both);

    std::cout << "Gradient momentum analysis: close.txt:\n";
    gradient_momentum_analysis("gradient_momentum_close.csv", coarse_steps, learning_rates,
                               momentums, training_samples_close, training_results_close);

    std::cout << "Gradient momentum analysis: walk.txt:\n";
    gradient_momentum_analysis("gradient_momentum_walk.csv", coarse_steps, learning_rates,
                               momentums, training_samples_walk, training_results_walk);

    std::cout << "Gradient momentum analysis: both:\n";
    gradient_momentum_analysis("gradient_momentum_both.csv", coarse_steps, learning_rates,
                               momentums, training_samples_both, training_results_both);

    return;
}

void gradient_momentum_analysis(const std::string output_filename, const std::vector<double>& coarse_steps,
                                const std::vector<double>& learning_rates, const std::vector<double>& momentums,
                                const std::vector<SamplesData>& samples_data,
                                const std::vector<DoaAngles>& correct_results_vector) {
    std::ofstream output_csv;
    DoaEstimator doa_estimator;
    const std::string output_name = output_dir + output_filename;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_name);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    // Add columns to csv file
    make_csv_columns(output_csv, "gradient_momentum");

    std::cout << "Total number of coarse_steps: " << coarse_steps.size() << "\n";
    std::cout << "Total number of learning_rates: " << learning_rates.size() << "\n";
    std::cout << "Total number of momentums: " << momentums.size() << "\n";

    for (std::size_t coarse_index = 0; coarse_index < coarse_steps.size(); coarse_index++) {
        double coarse_step = coarse_steps[coarse_index];
        std::cout << "Coarse step: " << coarse_step << "\n";
        for (std::size_t lr_index = 0; lr_index < learning_rates.size(); lr_index++) {
            double learning_rate = learning_rates[lr_index];
            std::cout << "Learning rate: " << learning_rate << "\n";
            for (std::size_t m_index = 0; m_index < momentums.size(); m_index++) {
                double momentum = momentums[m_index];
                std::cout << "Momentum: " << learning_rate << "\n";
                save_csv_info_for_every_sample(output_csv, "gradient_momentum", samples_data,
                                               correct_results_vector, MusicOptimization::gradient_momentum,
                                               coarse_step, learning_rate, momentum);
            }
        }
    }

    output_csv.close();
    std::cout << "\n";
    return;
}

// ****************************************************************************
// ****************************  UTILITY FUNCTONS *****************************
// ****************************************************************************

void get_training_data(std::vector<SamplesData>& training_samples_close, std::vector<DoaAngles>& training_results_close,
                       std::vector<SamplesData>& training_samples_walk, std::vector<DoaAngles>& training_results_walk,
                       std::vector<SamplesData>& training_samples_both, std::vector<DoaAngles>& training_results_both) {

    std::vector<SamplesData> samples_data_close, samples_data_walk;
    std::vector<DoaAngles> correct_results_close, correct_results_walk;
    read_files::get_iq_samples(samples_data_close, (iq_samples_dir + close_filename));
    read_files::get_iq_samples(samples_data_walk, (iq_samples_dir + walk_filename));
    read_files::get_music_result_angles(correct_results_close, (music_results_dir + close_filename));
    read_files::get_music_result_angles(correct_results_walk, (music_results_dir + walk_filename));

    for (std::size_t i = 0; i < samples_data_close.size(); i += training_stride) {
        training_samples_close.push_back(samples_data_close[i]);
        training_samples_walk.push_back(samples_data_walk[i]);
        training_results_close.push_back(correct_results_close[i]);
        training_results_walk.push_back(correct_results_walk[i]);
    }

    training_samples_both.insert(training_samples_both.end(), training_samples_close.begin(), training_samples_close.end());
    training_samples_both.insert(training_samples_both.end(), training_samples_walk.begin(), training_samples_walk.end());
    training_results_both.insert(training_results_both.end(), training_results_close.begin(), training_results_close.end());
    training_results_both.insert(training_results_both.end(), training_results_walk.begin(), training_results_walk.end());

    return;
}

void make_csv_columns(std::ofstream& output_csv, const std::string analysis_method) {
    output_csv << "coarse_step,";
    if (analysis_method == "gradient_simple") {
        output_csv << "learning_rate,";
    } else if (analysis_method == "gradient_momentum") {
        output_csv << "learning_rate,momentum,";
    }
    output_csv << "number_accurates,total_iterations,"
               << "accuracy,runtime,mae,mse,rmse\n";
}

void save_csv_info_for_every_sample(std::ofstream& output_csv, const std::string analysis_method,
                                    const std::vector<SamplesData>& samples_data,
                                    const std::vector<DoaAngles>& correct_results_vector,
                                    const MusicOptimization optimization, const double coarse_step,
                                    double learning_rate, double momentum) {
    auto double_precision = std::numeric_limits<long double>::digits10;
    std::vector<double> error_lengths;
    DoaEstimator doa_estimator;
    GradientSpecs gradient_specs = {1e-5, 1e-6, learning_rate, momentum};
    double coarse_step_pi = utility::angle_to_pi(coarse_step);
    double n_accurates = 0;
    int iterations = 0;
    ProgressBar progress_bar(samples_data.size());

    auto t1 = std::chrono::high_resolution_clock::now();

    for (std::size_t sample_index = 0; sample_index < samples_data.size(); sample_index++) {
        DoaAngles result_angles, correct_angles;
        correct_angles = correct_results_vector[sample_index];
        result_angles = doa_estimator.process_samples(samples_data[sample_index], DoaTechnique::music,
                                                      MusicSearch::coarse_grid, M_PI / 1800,
                                                      optimization, coarse_step_pi, gradient_specs);
        if (utility::is_equal_angles(correct_angles, result_angles, 4 * finer_step)) {
            n_accurates++;
        } else {
            correct_angles = utility::angles_to_degree(correct_angles);
            result_angles = utility::angles_to_degree(result_angles);
            double error = std::hypot((correct_angles.azimuth - result_angles.azimuth),
                                      (correct_angles.elevation - result_angles.elevation));
            error_lengths.push_back(error);
        }
        progress_bar.update();
        iterations++;
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float, std::milli> runtime = t2 - t1;

    double mae = 0, mse = 0, rmse = 0;
    if (error_lengths.size() > 0) {
        mae = stats::mae_double(error_lengths);
        mse = stats::mse_double(error_lengths);
        rmse = stats::rmse_double(mse);
    }

    // Save values to CSV.
    // Columns: index, number_accurates, total_iterations, accuracy, runtime, mae, mse, rmse
    output_csv << coarse_step << ",";
    if (analysis_method == "gradient_simple") {
        output_csv << learning_rate << ",";
    } else if (analysis_method == "gradient_momentum") {
        output_csv << learning_rate << ",";
        output_csv << momentum << ",";
    }
    output_csv << n_accurates << ",";
    output_csv << iterations << ",";
    output_csv << std::setprecision(double_precision) << (n_accurates / iterations) << ",";
    output_csv << runtime.count() << ",";
    output_csv << std::setprecision(double_precision) << mae << ",";
    output_csv << std::setprecision(double_precision) << mse << ",";
    output_csv << std::setprecision(double_precision) << rmse << "\n";

    std::cout << "\n";
    return;
}
