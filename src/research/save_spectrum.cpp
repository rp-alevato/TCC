#include "libdoa/doa_estimator.h"
#include "misc/read_data_files.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

void music_save_spectrum(SamplesData& samples_data, double grid_step, std::string file_name);

int main() {
    std::vector<SamplesData> samples_data_list;
    read_files::get_iq_samples(samples_data_list, "data/iq_samples/close.txt");
    // music_save_spectrum(samples_data_list[0], M_PI / 180, "close_sample_0000.csv");
    // music_save_spectrum(samples_data_list[21], M_PI / 180, "close_sample_0021.csv");
    // music_save_spectrum(samples_data_list[29], M_PI / 180, "close_sample_0029.csv");
    std::cout << samples_data_list.size() << "\n";
    music_save_spectrum(samples_data_list[300], M_PI / 180, "close_sample_0300.csv");
    // music_save_spectrum(samples_data_list[900], M_PI / 180, "close_sample_0900.csv");
    // music_save_spectrum(samples_data_list[1500], M_PI / 180, "close_sample_1500.csv");
    // music_save_spectrum(samples_data_list[2100], M_PI / 180, "close_sample_2100.csv");
}

void music_save_spectrum(SamplesData& samples_data, double const grid_step, std::string file_name) {
    // Data structure to export a csv
    std::vector<double> rows_labels;
    std::vector<double> cols_labels;
    std::vector<std::vector<double>> music_spectrum;
    DoaEstimator doa_estimator;
    int azimuth_max_iter = (int)(2 * M_PI / grid_step);
    int elevation_max_iter = (int)(M_PI / grid_step);

    std::cout << "Saving " + file_name + "...";

    // Load samples and calculate noise eigenvectors
    doa_estimator.process_samples(samples_data, DoaTechnique::music, MusicSearch::simple_grid, M_PI / 2);

    for (int elevation_index = -elevation_max_iter; elevation_index < elevation_max_iter; elevation_index++) {
        cols_labels.push_back(elevation_index * grid_step);
    }
    for (int azimuth_index = 0; azimuth_index < azimuth_max_iter; azimuth_index++) {
        double azimuth = azimuth_index * grid_step;
        rows_labels.push_back(azimuth);
        music_spectrum.push_back({});
        for (int elevation_index = -elevation_max_iter; elevation_index < elevation_max_iter; elevation_index++) {
            double elevation = elevation_index * grid_step;
            double result = doa_estimator.estimate_music_result({azimuth, elevation});
            music_spectrum.back().push_back(result);
        }
    }

    // Print to file
    std::ofstream spectrum_csv;
    spectrum_csv.open("data/music_spectrum_data/" + file_name);
    auto double_precision = std::numeric_limits<long double>::digits10 + 1;

    for (auto it = cols_labels.begin(); it != cols_labels.end(); it++) {
        spectrum_csv << std::setprecision(double_precision) << *it;
        if (std::next(it) != cols_labels.end()) {
            spectrum_csv << ",";
        }
    }
    spectrum_csv << "\n";

    for (std::size_t i = 0; i < rows_labels.size(); i++) {
        spectrum_csv << std::setprecision(double_precision) << rows_labels[i] << ",";
        for (std::size_t j = 0; j < music_spectrum[0].size(); j++) {
            spectrum_csv << std::setprecision(double_precision) << music_spectrum[i][j];
            if ((j + 1) < music_spectrum[0].size()) {
                spectrum_csv << ",";
            }
        }
        spectrum_csv << "\n";
    }
    spectrum_csv.close();

    std::cout << "  Done!\n";

    return;
}
