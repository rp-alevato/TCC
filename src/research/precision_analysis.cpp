#include "doa/estimator.h"
#include "misc/progress_bar.h"
#include "misc/read_data_files.h"
#include "misc/statistics.h"
#include "misc/utility.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// Step 1: Load all samples and music results
// Step 2: Get all results for ESPRIT and Silabs API
// Step 3: Calculate parameters

std::string fixed_filename(unsigned int i);
void precision_analysis(const std::string& filename);
void calc_save_output(const std::string& technique, std::ofstream& output_csv, std::vector<double>& results_len,
                      std::vector<double>& results_az, std::vector<double>& results_el);

int main() {
    for (int i = 1; i <= 5; i++) {
        precision_analysis(fixed_filename(i));
    }

    return 0;
}

void precision_analysis(const std::string& filename) {
    const std::string iq_samples_filename = "data/iq_samples/" + filename + ".txt";
    const std::string music_results_filename = "data/music_result_angles/" + filename + ".txt";
    const std::string output_filename = "data/experimental_results/precision/" + filename + ".csv";

    std::vector<SamplesData> samples_data;
    std::vector<DoaAngles> music_results;
    std::vector<double> music_results_len, music_results_az, music_results_el;
    std::vector<double> esprit_results_len, esprit_results_az, esprit_results_el;
    std::vector<double> sl_results_len, sl_results_az, sl_results_el;
    read_files::get_iq_samples(samples_data, iq_samples_filename);
    read_files::get_music_result_angles(music_results, music_results_filename);

    std::ofstream output_csv;
    if (std::filesystem::exists(output_filename)) {
        throw std::runtime_error("File " + output_filename + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_csv.open(output_filename);
    if (!output_csv.is_open()) {
        throw std::runtime_error("Error opening file " + output_filename);
    }

    std::cout << "Saving file " << output_filename << " ...\n";
    DoaEstimator estimator;
    for (auto result : music_results) {
        double azimuth = utility::angle_to_degree(result.azimuth);
        double elevation = utility::angle_to_degree(result.elevation);
        music_results_len.push_back(std::hypot(azimuth, elevation));
        music_results_az.push_back(azimuth);
        music_results_el.push_back(elevation);
    }
    for (auto sample : samples_data) {
        DoaAngles esprit_result = estimator.process_samples(sample, DoaTechnique::esprit);
        double azimuth = utility::angle_to_degree(esprit_result.azimuth);
        double elevation = utility::angle_to_degree(esprit_result.elevation);
        esprit_results_len.push_back(std::hypot(azimuth, elevation));
        esprit_results_az.push_back(utility::angle_to_degree(esprit_result.azimuth));
        esprit_results_el.push_back(utility::angle_to_degree(esprit_result.elevation));
        azimuth = sample.sl_azimuth;
        elevation = sample.sl_elevation;
        sl_results_len.push_back(std::hypot(azimuth, elevation));
        sl_results_az.push_back(azimuth);
        sl_results_el.push_back(elevation);
    }
    output_csv << "technique,"
               << "mean_len,mad_len,std_dev_len,mad_len_99p,std_dev_len_99p,mad_len_1p,std_dev_len_1p,"
               << "mean_az,mad_az,std_dev_az,mad_az_99p,std_dev_az_99p,mad_az_1p,std_dev_az_1p,"
               << "mean_el,mad_el,std_dev_el,mad_el_99p,std_dev_el_99p,mad_el_1p,std_dev_el_1p\n";
    calc_save_output("music", output_csv, music_results_len, music_results_az, music_results_el);
    calc_save_output("esprit", output_csv, esprit_results_len, esprit_results_az, esprit_results_el);
    calc_save_output("silabs", output_csv, sl_results_len, sl_results_az, sl_results_el);

    return;
}

void calc_save_output(const std::string& technique, std::ofstream& output_csv, std::vector<double>& results_len,
                      std::vector<double>& results_az, std::vector<double>& results_el) {

    auto double_precision = std::numeric_limits<double>::digits10;
    std::vector<double> errors_len, errors_az, errors_el;

    std::sort(results_len.begin(), results_len.end());
    std::sort(results_az.begin(), results_az.end());
    std::sort(results_el.begin(), results_el.end());

    double mean_len = stats::mean_double(results_len);
    double mean_az = stats::mean_double(results_az);
    double mean_el = stats::mean_double(results_el);

    for (std::size_t i = 0; i < results_len.size(); i++) {
        errors_len.push_back(results_len[i] - mean_len);
        errors_az.push_back(results_az[i] - mean_az);
        errors_el.push_back(results_el[i] - mean_el);
    }
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

    double mad_len = stats::mae_double(errors_len);
    double mad_az = stats::mae_double(errors_az);
    double mad_el = stats::mae_double(errors_el);
    double std_dev_len = stats::rmse_double(stats::mse_double(errors_len));
    double std_dev_az = stats::rmse_double(stats::mse_double(errors_az));
    double std_dev_el = stats::rmse_double(stats::mse_double(errors_el));
    double mad_len_99p = stats::mae_double(errors_len_99p);
    double mad_az_99p = stats::mae_double(errors_az_99p);
    double mad_el_99p = stats::mae_double(errors_el_99p);
    double std_dev_len_99p = stats::rmse_double(stats::mse_double(errors_len_99p));
    double std_dev_az_99p = stats::rmse_double(stats::mse_double(errors_az_99p));
    double std_dev_el_99p = stats::rmse_double(stats::mse_double(errors_el_99p));
    double mad_len_1p = stats::mae_double(errors_len_1p);
    double mad_az_1p = stats::mae_double(errors_az_1p);
    double mad_el_1p = stats::mae_double(errors_el_1p);
    double std_dev_len_1p = stats::rmse_double(stats::mse_double(errors_len_1p));
    double std_dev_az_1p = stats::rmse_double(stats::mse_double(errors_az_1p));
    double std_dev_el_1p = stats::rmse_double(stats::mse_double(errors_el_1p));

    // Columns: technique, mean_len, median_len, mad_len, std_dev_len, mad_len_99p, std_dev_len_99p, mad_len_1p, std_dev_len_1p
    //                     mean_az, median_az, mad_az, std_dev_az, mad_az_99p, std_dev_az_99p, mad_az_1p, std_dev_az_1p
    //                     mean_el, median_el, mad_el, std_dev_el, mad_el_99p, std_dev_el_99p, mad_el_1p, std_dev_el_1p
    output_csv << technique << ",";
    output_csv << std::setprecision(double_precision) << mean_len << ",";
    output_csv << std::setprecision(double_precision) << mad_len << ",";
    output_csv << std::setprecision(double_precision) << std_dev_len << ",";
    output_csv << std::setprecision(double_precision) << mad_len_99p << ",";
    output_csv << std::setprecision(double_precision) << std_dev_len_99p << ",";
    output_csv << std::setprecision(double_precision) << mad_len_1p << ",";
    output_csv << std::setprecision(double_precision) << std_dev_len_1p << ",";

    output_csv << std::setprecision(double_precision) << mean_az << ",";
    output_csv << std::setprecision(double_precision) << mad_az << ",";
    output_csv << std::setprecision(double_precision) << std_dev_az << ",";
    output_csv << std::setprecision(double_precision) << mad_az_99p << ",";
    output_csv << std::setprecision(double_precision) << std_dev_az_99p << ",";
    output_csv << std::setprecision(double_precision) << mad_az_1p << ",";
    output_csv << std::setprecision(double_precision) << std_dev_az_1p << ",";

    output_csv << std::setprecision(double_precision) << mean_el << ",";
    output_csv << std::setprecision(double_precision) << mad_el << ",";
    output_csv << std::setprecision(double_precision) << std_dev_el << ",";
    output_csv << std::setprecision(double_precision) << mad_el_99p << ",";
    output_csv << std::setprecision(double_precision) << std_dev_el_99p << ",";
    output_csv << std::setprecision(double_precision) << mad_el_1p << ",";
    output_csv << std::setprecision(double_precision) << std_dev_el_1p << "\n";

    return;
}

std::string fixed_filename(unsigned int i) {
    return "fixed_" + std::to_string(i);
}
