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
const std::string close_filename = "close.txt";
const std::string walk_filename = "office_walk.txt";
const std::string output_dir = "data/experimental_results/hyperparameters/";

static constexpr std::size_t training_stride = 5;
static constexpr double fine_step = M_PI / 900;

void coarse_fine_search_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                                 const std::vector<DoaAngles>& correct_results);
void gradient_simple_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                              const std::vector<DoaAngles>& correct_results);
void gradient_adapt_lr_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                                const std::vector<DoaAngles>& correct_results);
void gradient_momentum_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                                const std::vector<DoaAngles>& correct_results);
void gradient_nesterov_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                                const std::vector<DoaAngles>& correct_results);
// Utility functions
void get_training_data(std::vector<SamplesData>& training_samples, std::vector<DoaAngles>& training_results);
void make_csv_columns(std::ofstream& output_csv, const std::string analysis_method);
void save_csv_info_for_every_sample(std::ofstream& output_csv, const std::string analysis_method,
                                    const std::vector<SamplesData>& samples_data,
                                    const std::vector<DoaAngles>& correct_results,
                                    const MusicOptimization optimization, const double coarse_step,
                                    const double learning_rate, const double momentum);

int main() {
    std::vector<SamplesData> training_samples;
    std::vector<DoaAngles> training_results;

    get_training_data(training_samples, training_results);

    std::cout << "Coarse-fine search analysis:\n";
    coarse_fine_search_analysis("coarse_fine_search.csv", training_samples, training_results);

    std::cout << "Gradient simple analysis:\n";
    gradient_simple_analysis("gradient_simple.csv", training_samples, training_results);

    std::cout << "Gradient adapt lr analysis:\n";
    gradient_adapt_lr_analysis("gradient_adapt_lr.csv", training_samples, training_results);

    std::cout << "Gradient momentum analysis:\n";
    gradient_momentum_analysis("gradient_momentum.csv", training_samples, training_results);

    std::cout << "Gradient nesterov analysis:\n";
    gradient_nesterov_analysis("gradient_nesterov.csv", training_samples, training_results);

    return 0;
}

void coarse_fine_search_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                                 const std::vector<DoaAngles>& correct_results) {
    std::ofstream output_csv;
    DoaEstimator estimator;
    std::vector<double> coarse_steps;
    const std::string output_name = output_dir + output_filename;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_name);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    for (std::size_t i = 2; i <= 9; i++) {
        coarse_steps.push_back(i);
    }

    make_csv_columns(output_csv, "coarse_fine_grid");

    std::cout << "Total number of coarse_steps: " << coarse_steps.size() << "\n";

    for (std::size_t coarse_index = 0; coarse_index < coarse_steps.size(); coarse_index++) {
        double coarse_step = coarse_steps[coarse_index];
        std::cout << "cs: " << std::setw(2) << coarse_step << "\n";
        save_csv_info_for_every_sample(output_csv, "coarse_fine_grid", samples_data,
                                       correct_results, MusicOptimization::fine_grid_search,
                                       coarse_step, 0, 0);
    }

    output_csv.close();
    std::cout << "\n";
    return;
}

void gradient_simple_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                              const std::vector<DoaAngles>& correct_results) {
    std::ofstream output_csv;
    DoaEstimator estimator;
    std::vector<double> coarse_steps, learning_rates;
    const std::string output_name = output_dir + output_filename;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_name);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    for (std::size_t i = 2; i <= 9; i += 1) {
        coarse_steps.push_back(i);
    }

    for (double i = 0.05; i <= 0.091; i += 0.02) {
        learning_rates.push_back(i);
    }
    for (double i = 0.1; i <= 0.51; i += 0.1) {
        learning_rates.push_back(i);
    }

    make_csv_columns(output_csv, "gradient_simple");

    std::cout << "Total number of coarse_steps: " << coarse_steps.size() << "\n";
    std::cout << "Total number of learning_rates: " << learning_rates.size() << "\n";

    for (std::size_t coarse_index = 0; coarse_index < coarse_steps.size(); coarse_index++) {
        double coarse_step = coarse_steps[coarse_index];
        for (std::size_t lr_index = 0; lr_index < learning_rates.size(); lr_index++) {
            double learning_rate = learning_rates[lr_index];
            std::cout << "cs: " << std::setw(2) << coarse_step << "    "
                      << "lr: " << std::setw(4) << learning_rate << "\n";
            save_csv_info_for_every_sample(output_csv, "gradient_simple", samples_data,
                                           correct_results, MusicOptimization::gradient_simple,
                                           coarse_step, learning_rate, 0);
        }
    }

    output_csv.close();
    std::cout << "\n";
    return;
}

void gradient_adapt_lr_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                                const std::vector<DoaAngles>& correct_results) {
    std::ofstream output_csv;
    DoaEstimator estimator;
    std::vector<double> coarse_steps, learning_rates;
    const std::string output_name = output_dir + output_filename;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_name);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    for (std::size_t i = 2; i <= 9; i += 1) {
        coarse_steps.push_back(i);
    }

    for (double i = 0.05; i <= 0.091; i += 0.02) {
        learning_rates.push_back(i);
    }
    for (double i = 0.1; i <= 0.91; i += 0.1) {
        learning_rates.push_back(i);
    }

    make_csv_columns(output_csv, "gradient_adapt_lr");

    std::cout << "Total number of coarse_steps: " << coarse_steps.size() << "\n";
    std::cout << "Total number of learning_rates: " << learning_rates.size() << "\n";

    for (std::size_t coarse_index = 0; coarse_index < coarse_steps.size(); coarse_index++) {
        double coarse_step = coarse_steps[coarse_index];
        for (std::size_t lr_index = 0; lr_index < learning_rates.size(); lr_index++) {
            double learning_rate = learning_rates[lr_index];
            std::cout << "cs: " << std::setw(2) << coarse_step << "    "
                      << "lr: " << std::setw(4) << learning_rate << "\n";
            save_csv_info_for_every_sample(output_csv, "gradient_adapt_lr", samples_data,
                                           correct_results, MusicOptimization::gradient_adapt_lr,
                                           coarse_step, learning_rate, 0);
        }
    }

    output_csv.close();
    std::cout << "\n";
    return;
}

void gradient_momentum_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                                const std::vector<DoaAngles>& correct_results) {
    std::ofstream output_csv;
    DoaEstimator estimator;
    std::vector<double> coarse_steps, learning_rates, momentums;
    const std::string output_name = output_dir + output_filename;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_name);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    for (std::size_t i = 2; i <= 9; i += 1) {
        coarse_steps.push_back(i);
    }

    for (double i = 0.01; i <= 0.091; i += 0.02) {
        learning_rates.push_back(i);
    }
    for (double i = 0.1; i <= 0.21; i += 0.1) {
        learning_rates.push_back(i);
    }

    for (double i = 0.70; i <= 0.951; i += 0.05) {
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
                std::cout << "cs: " << std::setw(2) << std::right << coarse_step << "    "
                          << "lr: " << std::setw(5) << std::left << learning_rate << "    "
                          << "mm: " << std::setw(4) << std::left << momentum << "\n";
                save_csv_info_for_every_sample(output_csv, "gradient_momentum", samples_data,
                                               correct_results, MusicOptimization::gradient_momentum,
                                               coarse_step, learning_rate, momentum);
            }
        }
    }

    output_csv.close();
    std::cout << "\n";
    return;
}

void gradient_nesterov_analysis(const std::string output_filename, const std::vector<SamplesData>& samples_data,
                                const std::vector<DoaAngles>& correct_results) {
    std::ofstream output_csv;
    DoaEstimator estimator;
    std::vector<double> coarse_steps, learning_rates, momentums;
    const std::string output_name = output_dir + output_filename;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_name);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    for (std::size_t i = 2; i <= 9; i += 1) {
        coarse_steps.push_back(i);
    }

    for (double i = 0.01; i <= 0.091; i += 0.02) {
        learning_rates.push_back(i);
    }
    for (double i = 0.1; i <= 0.21; i += 0.1) {
        learning_rates.push_back(i);
    }

    for (double i = 0.70; i <= 0.951; i += 0.05) {
        momentums.push_back(i);
    }

    make_csv_columns(output_csv, "gradient_nesterov");

    std::cout << "Total number of coarse_steps: " << coarse_steps.size() << "\n";
    std::cout << "Total number of learning_rates: " << learning_rates.size() << "\n";
    std::cout << "Total number of momentums: " << momentums.size() << "\n";

    for (std::size_t coarse_index = 0; coarse_index < coarse_steps.size(); coarse_index++) {
        double coarse_step = coarse_steps[coarse_index];
        for (std::size_t lr_index = 0; lr_index < learning_rates.size(); lr_index++) {
            double learning_rate = learning_rates[lr_index];
            for (std::size_t m_index = 0; m_index < momentums.size(); m_index++) {
                double momentum = momentums[m_index];
                std::cout << "cs: " << std::setw(2) << std::right << coarse_step << "    "
                          << "lr: " << std::setw(5) << std::left << learning_rate << "    "
                          << "mm: " << std::setw(4) << std::left << momentum << "\n";
                save_csv_info_for_every_sample(output_csv, "gradient_nesterov", samples_data,
                                               correct_results, MusicOptimization::gradient_nesterov,
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
    auto double_precision = std::numeric_limits<double>::digits10;
    std::vector<double> errors_len, errors_az, errors_el;
    std::vector<DoaAngles> results;
    DoaEstimator estimator;
    GradientSpecs gradient_specs = {1e-5, 1.5e-8, learning_rate, momentum};
    double coarse_step_pi = utility::angle_to_pi(coarse_step);

    // Estimate angles
    auto t1 = std::chrono::high_resolution_clock::now();
    for (std::size_t sample_index = 0; sample_index < samples_data.size(); sample_index++) {
        results.push_back(estimator.process_samples(samples_data[sample_index], DoaTechnique::music,
                                                    MusicSearch::coarse_grid, fine_step,
                                                    optimization, coarse_step_pi, gradient_specs));
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float, std::milli> delta_t = t2 - t1;

    // Calculate errors for length, azimuth and elevation
    int n_accurates = 0;
    int n_max_iterations = 0;
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
        if (estimator.get_was_max_iterations()) {
            n_max_iterations++;
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
    if (analysis_method == "gradient_simple" || analysis_method == "gradient_adapt_lr") {
        output_csv << learning_rate << ",";
    } else if (analysis_method == "gradient_momentum" || analysis_method == "gradient_nesterov") {
        output_csv << learning_rate << ",";
        output_csv << momentum << ",";
    }

    // Columns: runtime, n_max_iterations, accuracy,
    //          mae_len, rmse_len, mae_len_99p, rmse_len_99p, mae_len_1p, rmse_len_1p,
    //          mae_az, rmse_az, mae_az_99p, rmse_az_99p, mae_az_1p, rmse_az_1p,
    //          mae_el, rmse_el, mae_el_99p, rmse_el_99p, mae_el_1p, rmse_el_1p,
    output_csv << std::setprecision(double_precision) << runtime << ",";
    output_csv << n_max_iterations << ",";
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

void get_training_data(std::vector<SamplesData>& training_samples, std::vector<DoaAngles>& training_results) {

    std::vector<SamplesData> samples_data_walk;
    std::vector<DoaAngles> correct_results_walk;
    read_files::get_iq_samples(samples_data_walk, (iq_samples_dir + walk_filename));
    read_files::get_music_result_angles(correct_results_walk, (music_results_dir + walk_filename));

    for (std::size_t i = 0; i < samples_data_walk.size(); i += training_stride) {
        training_samples.push_back(samples_data_walk[i]);
        training_results.push_back(correct_results_walk[i]);
    }

    return;
}

void make_csv_columns(std::ofstream& output_csv, const std::string analysis_method) {
    output_csv << "coarse_step,";
    if (analysis_method == "gradient_simple" || analysis_method == "gradient_adapt_lr") {
        output_csv << "learning_rate,";
    } else if (analysis_method == "gradient_momentum" || analysis_method == "gradient_nesterov") {
        output_csv << "learning_rate,momentum,";
    }
    output_csv << "runtime,n_max_iterations,accuracy,"
               << "mae_len,rmse_len,mae_len_99p,rmse_len_99p,mae_len_1p,rmse_len_1p,"
               << "mae_az,rmse_az,mae_az_99p,rmse_az_99p,mae_az_1p,rmse_az_1p,"
               << "mae_el,rmse_el,mae_el_99p,rmse_el_99p,mae_el_1p,rmse_el_1p\n";
    return;
}
