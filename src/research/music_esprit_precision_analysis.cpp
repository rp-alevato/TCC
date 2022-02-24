#include "libdoa/doa_estimator.h"
#include "misc/progress_bar.h"
#include "misc/read_data_files.h"
#include "misc/statistics.h"
#include "misc/utility.h"

#include <algorithm>
#include <chrono>
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

int main() {
    for (int i = 1; i <= 5; i++) {
        precision_analysis(fixed_filename(i));
    }

    return 0;
}

void precision_analysis(const std::string& filename) {
    const std::string iq_samples_filename = "data/iq_samples/" + filename;
    const std::string music_results_filename = "data/music_result_angles/" + filename;
    const std::string output_filename = "data/experimental_results/music_esprit_precision/" + filename;

    std::vector<SamplesData> samples_data;
    std::vector<DoaAngles> music_results;
    read_files::get_iq_samples(samples_data, iq_samples_filename);
    read_files::get_music_result_angles(music_results, iq_samples_filename);

    std::ofstream output_csv;

    output_csv << "technique,mean,median,variance,std_dev,mae_sl_api,mse_sl_api,rmse_sl_api,\n";
}

std::string fixed_filename(unsigned int i) {
    return "fixed_" + std::to_string(i) + ".txt";
}
