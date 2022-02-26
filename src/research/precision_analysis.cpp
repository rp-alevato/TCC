#include "libdoa/doa_estimator.h"
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
    DoaEstimator doa_estimator;
    for (auto result : music_results) {
        double azimuth = utility::angle_to_degree(result.azimuth);
        double elevation = utility::angle_to_degree(result.elevation);
        music_results_len.push_back(std::hypot(azimuth, elevation));
        music_results_az.push_back(azimuth);
        music_results_el.push_back(elevation);
    }
    for (auto sample : samples_data) {
        DoaAngles esprit_result = doa_estimator.process_samples(sample, DoaTechnique::esprit);
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

    output_csv << "technique,mean_len,median_len,variance_len,std_dev_len,mean_az,median_az,variance_az,std_dev_az,mean_el,median_el,variance_el,std_dev_el\n";
    calc_save_output("music", output_csv, music_results_len, music_results_az, music_results_el);
    calc_save_output("esprit", output_csv, esprit_results_len, esprit_results_az, esprit_results_el);
    calc_save_output("sl", output_csv, sl_results_len, sl_results_az, sl_results_el);

    return;
}

void calc_save_output(const std::string& technique, std::ofstream& output_csv, std::vector<double>& results_len,
                      std::vector<double>& results_az, std::vector<double>& results_el) {

    auto double_precision = std::numeric_limits<long double>::digits10;

    std::sort(results_len.begin(), results_len.end());
    std::sort(results_az.begin(), results_az.end());
    std::sort(results_el.begin(), results_el.end());
    double mean_len = stats::mean_double(results_len);
    double mean_az = stats::mean_double(results_az);
    double mean_el = stats::mean_double(results_el);
    double median_len = stats::median_sorted_double(results_len);
    double median_az = stats::median_sorted_double(results_az);
    double median_el = stats::median_sorted_double(results_el);
    double variance_len = stats::variance_double(results_len, mean_len);
    double variance_az = stats::variance_double(results_az, mean_az);
    double variance_el = stats::variance_double(results_el, mean_el);
    double std_dev_len = stats::std_deviation_double(variance_len);
    double std_dev_az = stats::std_deviation_double(variance_az);
    double std_dev_el = stats::std_deviation_double(variance_el);

    // Columns: technique, mean_len, median_len, variance_len, std_dev_len, mean_az, median_az, variance_az, std_dev_az, mean_el, median_el, variance_el, std_dev_el
    output_csv << technique << ",";
    output_csv << std::setprecision(double_precision) << mean_len << ",";
    output_csv << std::setprecision(double_precision) << median_len << ",";
    output_csv << std::setprecision(double_precision) << variance_len << ",";
    output_csv << std::setprecision(double_precision) << std_dev_len << ",";
    output_csv << std::setprecision(double_precision) << mean_az << ",";
    output_csv << std::setprecision(double_precision) << median_az << ",";
    output_csv << std::setprecision(double_precision) << variance_az << ",";
    output_csv << std::setprecision(double_precision) << std_dev_az << ",";
    output_csv << std::setprecision(double_precision) << mean_el << ",";
    output_csv << std::setprecision(double_precision) << median_el << ",";
    output_csv << std::setprecision(double_precision) << variance_el << ",";
    output_csv << std::setprecision(double_precision) << std_dev_el << "\n";

    return;
}

std::string fixed_filename(unsigned int i) {
    return "fixed_" + std::to_string(i);
}
