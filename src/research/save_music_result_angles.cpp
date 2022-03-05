#include "aoa/estimator.h"
#include "misc/progress_bar.h"
#include "misc/read_data_files.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

void save_music_result_angles(std::string iq_file_name, std::string music_result_angles_file_name);

int main() {
    save_music_result_angles("close.txt", "close_finer.txt");
    save_music_result_angles("office_walk.txt", "office_walk_finer.txt");
    save_music_result_angles("fixed_1.txt", "fixed_1_finer.txt");
    save_music_result_angles("fixed_2.txt", "fixed_2_finer.txt");
    save_music_result_angles("fixed_3.txt", "fixed_3_finer.txt");
    save_music_result_angles("fixed_4.txt", "fixed_4_finer.txt");
    save_music_result_angles("fixed_5.txt", "fixed_5_finer.txt");
    return 0;
}

void save_music_result_angles(std::string iq_file_name, std::string music_result_angles_file_name) {
    std::vector<SamplesData> samples_data_vector;
    read_files::get_iq_samples(samples_data_vector, "data/iq_samples/" + iq_file_name);

    ProgressBar progress_bar(samples_data_vector.size());
    AoaEstimator estimator;
    AoaAngles angles;
    std::ofstream output_file;
    auto double_precision = std::numeric_limits<double>::digits10;
    const std::string output_name = "data/music_result_angles/" + music_result_angles_file_name;

    if (std::filesystem::exists(output_name)) {
        throw std::runtime_error("File " + output_name + " exists.\n\t   Please remove file if you want to overwrite it.\n");
    }
    output_file.open(output_name);
    if (!output_file.is_open()) {
        throw std::runtime_error("Error opening file " + output_name);
    }

    std::cout << "Saving music result angles for " << music_result_angles_file_name << "\n";
    for (auto samples_data : samples_data_vector) {
        progress_bar.update();
        angles = estimator.process_samples(samples_data, AoaTechnique::music, MusicSearch::simple_grid, M_PI / 1800);
        output_file << std::setprecision(double_precision) << "("
                    << angles.azimuth << "," << angles.elevation << ")\n";
    }
    output_file.close();

    std::cout << "\nDone!\n\n";

    return;
}
