#include "libdoa/doa_estimator.h"
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
void complete_finer_grid_analysis();
void finer_grid_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                         const std::vector<DoaAngles>& correct_results_vector);
//  Gradient simple functions
void complete_gradient_simple_analysis();
void gradient_simple_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                              const std::vector<DoaAngles>& correct_results_vector);
//  Gradient momentum functions
void complete_gradient_momentum_analysis();
void gradient_momentum_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
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
                                    const double learning_rate, const double momentum);

int main() {
    complete_finer_grid_analysis();
    complete_gradient_simple_analysis();
    complete_gradient_momentum_analysis();
    return 0;
}

// ****************************************************************************
// **************************  FINER GRID FUNCTIONS ***************************
// ****************************************************************************

void complete_finer_grid_analysis() {
    std::vector<SamplesData> training_samples_close, training_samples_walk, training_samples_both;
    std::vector<DoaAngles> training_results_close, training_results_walk, training_results_both;

    get_training_data(training_samples_close, training_results_close,
                      training_samples_walk, training_results_walk,
                      training_samples_both, training_results_both);

    std::cout << "Finer grid analysis: close.txt:\n";
    finer_grid_analysis("finer_grid_close.csv", training_samples_close, training_results_close);
    std::cout << "Finer grid analysis: walk.txt:\n";
    finer_grid_analysis("finer_grid_walk.csv", training_samples_walk, training_results_walk);
    std::cout << "Finer grid analysis: both:\n";
    finer_grid_analysis("finer_grid_both.csv", training_samples_both, training_results_both);

    return;
}

void finer_grid_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                         const std::vector<DoaAngles>& correct_results_vector) {
    std::ofstream output_csv;
    DoaEstimator doa_estimator;
    std::vector<double> coarse_steps;
    const std::string output_name = output_dir + output_filename;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_name);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    for (std::size_t i = 1; i <= 10; i++) {
        coarse_steps.push_back(i);
    }

    make_csv_columns(output_csv, "finer_grid");

    std::cout << "Total number of coarse_steps: " << coarse_steps.size() << "\n";

    for (std::size_t coarse_index = 0; coarse_index < coarse_steps.size(); coarse_index++) {
        double coarse_step = coarse_steps[coarse_index];
        std::cout << "cs: " << std::setw(2) << coarse_step << "\n";
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

void complete_gradient_simple_analysis() {
    std::vector<SamplesData> training_samples_close, training_samples_walk, training_samples_both;
    std::vector<DoaAngles> training_results_close, training_results_walk, training_results_both;

    get_training_data(training_samples_close, training_results_close,
                      training_samples_walk, training_results_walk,
                      training_samples_both, training_results_both);

    std::cout << "Gradient simple analysis: close.txt:\n";
    gradient_simple_analysis("gradient_simple_close.csv", training_samples_close, training_results_close);
    std::cout << "Gradient simple analysis: walk.txt:\n";
    gradient_simple_analysis("gradient_simple_walk.csv", training_samples_walk, training_results_walk);
    std::cout << "Gradient simple analysis: both:\n";
    gradient_simple_analysis("gradient_simple_both.csv", training_samples_both, training_results_both);

    return;
}

void gradient_simple_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                              const std::vector<DoaAngles>& correct_results_vector) {
    std::ofstream output_csv;
    DoaEstimator doa_estimator;
    std::vector<double> coarse_steps, learning_rates;
    const std::string output_name = output_dir + output_filename;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_name);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    for (std::size_t i = 3; i <= 10; i++) {
        coarse_steps.push_back(i);
    }
    for (double i = 0.008; i <= 0.0095; i += 0.001) {
        learning_rates.push_back(i);
    }
    for (double i = 0.01; i <= 0.105; i += 0.01) {
        learning_rates.push_back(i);
    }
    learning_rates.push_back(0.2);

    make_csv_columns(output_csv, "gradient_simple");

    std::cout << "Total number of coarse_steps: " << coarse_steps.size() << "\n";
    std::cout << "Total number of learning_rates: " << learning_rates.size() << "\n";

    for (std::size_t coarse_index = 0; coarse_index < coarse_steps.size(); coarse_index++) {
        double coarse_step = coarse_steps[coarse_index];
        for (std::size_t lr_index = 0; lr_index < learning_rates.size(); lr_index++) {
            double learning_rate = learning_rates[lr_index];
            std::cout << "cs: " << std::setw(2) << coarse_step << "    "
                      << "lr: " << std::setw(5) << learning_rate << "\n";
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

void complete_gradient_momentum_analysis() {
    std::vector<SamplesData> training_samples_close, training_samples_walk, training_samples_both;
    std::vector<DoaAngles> training_results_close, training_results_walk, training_results_both;

    get_training_data(training_samples_close, training_results_close,
                      training_samples_walk, training_results_walk,
                      training_samples_both, training_results_both);

    std::cout << "Gradient momentum analysis: close.txt:\n";
    gradient_momentum_analysis("gradient_momentum_close.csv", training_samples_close, training_results_close);
    std::cout << "Gradient momentum analysis: walk.txt:\n";
    gradient_momentum_analysis("gradient_momentum_walk.csv", training_samples_walk, training_results_walk);
    std::cout << "Gradient momentum analysis: both:\n";
    gradient_momentum_analysis("gradient_momentum_both.csv", training_samples_both, training_results_both);

    return;
}

void gradient_momentum_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                                const std::vector<DoaAngles>& correct_results_vector) {
    std::ofstream output_csv;
    DoaEstimator doa_estimator;
    std::vector<double> coarse_steps, learning_rates, momentums;
    const std::string output_name = output_dir + output_filename;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_name);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    for (std::size_t i = 3; i <= 10; i++) {
        coarse_steps.push_back(i);
    }
    for (double i = 0.01; i <= 0.075; i += 0.01) {
        learning_rates.push_back(i);
    }
    for (double i = 0.35; i <= 0.955; i += 0.1) {
        momentums.push_back(i);
    }

    make_csv_columns(output_csv, "gradient_momentum");

    std::cout << "Total number of coarse_steps: " << coarse_steps.size() << "\n";
    std::cout << "Total number of learning_rates: " << learning_rates.size() << "\n";
    std::cout << "Total number of momentums: " << momentums.size() << "\n";

    for (std::size_t coarse_index = 0; coarse_index < coarse_steps.size(); coarse_index++) {
        double coarse_step = coarse_steps[coarse_index];
        for (std::size_t lr_index = 0; lr_index < learning_rates.size(); lr_index++) {
            double learning_rate = learning_rates[lr_index];
            for (std::size_t m_index = 0; m_index < momentums.size(); m_index++) {
                double momentum = momentums[m_index];
                std::cout << "cs: " << std::setw(2) << coarse_step << "    "
                          << "lr: " << std::setw(4) << learning_rate << "    "
                          << "mm: " << std::setw(4) << momentum << "\n";
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

void save_csv_info_for_every_sample(std::ofstream& output_csv, const std::string analysis_method,
                                    const std::vector<SamplesData>& samples_data,
                                    const std::vector<DoaAngles>& correct_results,
                                    const MusicOptimization optimization, const double coarse_step,
                                    const double learning_rate, const double momentum) {
    auto double_precision = std::numeric_limits<long double>::digits10;
    std::vector<double> errors_len, errors_az, errors_el;
    std::vector<DoaAngles> results;
    DoaEstimator estimator;
    GradientSpecs gradient_specs = {1e-9, 1e-9, learning_rate, momentum};
    double coarse_step_pi = utility::angle_to_pi(coarse_step);

    // Estimate angles
    auto t1 = std::chrono::high_resolution_clock::now();
    for (std::size_t sample_index = 0; sample_index < samples_data.size(); sample_index++) {
        results.push_back(estimator.process_samples(samples_data[sample_index], DoaTechnique::music,
                                                    MusicSearch::coarse_grid, finer_step,
                                                    optimization, coarse_step_pi, gradient_specs));
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float, std::milli> delta_t = t2 - t1;

    // Calculate errors for length, azimuth and elevation
    int n_accurates = 0;
    for (std::size_t i = 0; i < results.size(); i++) {
        double result_az = utility::angle_to_degree(results[i].azimuth);
        double result_el = utility::angle_to_degree(results[i].elevation);
        double correct_result_az = utility::angle_to_degree(correct_results[i].azimuth);
        double correct_result_el = utility::angle_to_degree(correct_results[i].elevation);
        double error_az = correct_result_az - result_az;
        double error_el = correct_result_el - result_el;
        double error_len = std::hypot(error_az, error_el);
        errors_az.push_back(error_az);
        errors_el.push_back(error_el);
        errors_len.push_back(error_len);
        if (error_len < 0.2) {
            n_accurates++;
        }
    }

    // Get first 99% of errors and last 1% of errors
    std::sort(errors_len.begin(), errors_len.end());
    std::sort(errors_az.begin(), errors_az.end());
    std::sort(errors_el.begin(), errors_el.end());
    auto end_len_99p = errors_len.begin() + (0.99 * errors_len.size());
    std::vector<double> errors_len_99p(errors_len.begin(), end_len_99p);
    std::vector<double> errors_len_1p(end_len_99p, errors_len.end());
    auto end_az_99p = errors_az.begin() + (0.99 * errors_az.size());
    std::vector<double> errors_az_99p(errors_az.begin(), end_az_99p);
    std::vector<double> errors_az_1p(end_az_99p, errors_az.end());
    auto end_el_99p = errors_el.begin() + (0.99 * errors_el.size());
    std::vector<double> errors_el_99p(errors_el.begin(), end_el_99p);
    std::vector<double> errors_el_1p(end_el_99p, errors_el.end());

    // Calculate statistics
    double runtime = delta_t.count();
    double accuracy = ((double)n_accurates / results.size());
    double mae_len = stats::mae_double(errors_len);
    double mae_len_99p = stats::mae_double(errors_len_99p);
    double mae_len_1p = stats::mae_double(errors_len_1p);
    double mae_az = stats::mae_double(errors_az);
    double mae_az_99p = stats::mae_double(errors_az_99p);
    double mae_az_1p = stats::mae_double(errors_az_1p);
    double mae_el = stats::mae_double(errors_el);
    double mae_el_99p = stats::mae_double(errors_el_99p);
    double mae_el_1p = stats::mae_double(errors_el_1p);
    double rmse_len = stats::rmse_double(stats::mse_double(errors_len));
    double rmse_len_99p = stats::rmse_double(stats::mse_double(errors_len_99p));
    double rmse_len_1p = stats::rmse_double(stats::mse_double(errors_len_1p));
    double rmse_az = stats::rmse_double(stats::mse_double(errors_az));
    double rmse_az_99p = stats::rmse_double(stats::mse_double(errors_az_99p));
    double rmse_az_1p = stats::rmse_double(stats::mse_double(errors_az_1p));
    double rmse_el = stats::rmse_double(stats::mse_double(errors_el));
    double rmse_el_99p = stats::rmse_double(stats::mse_double(errors_el_99p));
    double rmse_el_1p = stats::rmse_double(stats::mse_double(errors_el_1p));

    // Save values to CSV.
    output_csv << coarse_step << ",";
    if (analysis_method == "gradient_simple") {
        output_csv << learning_rate << ",";
    } else if (analysis_method == "gradient_momentum") {
        output_csv << learning_rate << ",";
        output_csv << momentum << ",";
    }

    // Columns: runtime, accuracy,
    //          mae_len, rmse_len, mae_len_99p, rmse_len_99p, mae_len_1p, rmse_len_1p,
    //          mae_az, rmse_az, mae_az_99p, rmse_az_99p, mae_az_1p, rmse_az_1p,
    //          mae_el, rmse_el, mae_el_99p, rmse_el_99p, mae_el_1p, rmse_el_1p,
    output_csv << std::setprecision(double_precision) << runtime << ",";
    output_csv << std::setprecision(double_precision) << accuracy << ",";

    output_csv << std::setprecision(double_precision) << mae_len << ",";
    output_csv << std::setprecision(double_precision) << rmse_len << ",";
    output_csv << std::setprecision(double_precision) << mae_len_99p << ",";
    output_csv << std::setprecision(double_precision) << rmse_len_99p << ",";
    output_csv << std::setprecision(double_precision) << mae_len_1p << ",";
    output_csv << std::setprecision(double_precision) << rmse_len_1p << ",";

    output_csv << std::setprecision(double_precision) << mae_az << ",";
    output_csv << std::setprecision(double_precision) << rmse_az << ",";
    output_csv << std::setprecision(double_precision) << mae_az_99p << ",";
    output_csv << std::setprecision(double_precision) << rmse_az_99p << ",";
    output_csv << std::setprecision(double_precision) << mae_az_1p << ",";
    output_csv << std::setprecision(double_precision) << rmse_az_1p << ",";

    output_csv << std::setprecision(double_precision) << mae_el << ",";
    output_csv << std::setprecision(double_precision) << rmse_el << ",";
    output_csv << std::setprecision(double_precision) << mae_el_99p << ",";
    output_csv << std::setprecision(double_precision) << rmse_el_99p << ",";
    output_csv << std::setprecision(double_precision) << mae_el_1p << ",";
    output_csv << std::setprecision(double_precision) << rmse_el_1p << "\n";

    return;
}

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
    output_csv << "runtime,accuracy,"
               << "mae_len,rmse_len,mae_len_99p,rmse_len_99p,mae_len_1p,rmse_len_1p,"
               << "mae_az,rmse_az,mae_az_99p,rmse_az_99p,mae_az_1p,rmse_az_1p,"
               << "mae_el,rmse_el,mae_el_99p,rmse_el_99p,mae_el_1p,rmse_el_1p\n";
}
