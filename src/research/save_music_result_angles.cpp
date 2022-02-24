#include "libdoa/doa_estimator.h"
#include "misc/progress_bar.h"
#include "misc/read_data_files.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

void save_music_result_angles(std::string iq_file_name, std::string music_result_angles_file_name);

int main() {
    // save_music_result_angles("close.txt", "close.txt");
    // save_music_result_angles("office_walk.txt", "office_walk.txt");
    // save_music_result_angles("fixed_1.txt", "fixed_1.txt");
    // save_music_result_angles("fixed_2.txt", "fixed_2.txt");
    // save_music_result_angles("fixed_3.txt", "fixed_3.txt");
    // save_music_result_angles("fixed_4.txt", "fixed_4.txt");
    // save_music_result_angles("fixed_5.txt", "fixed_5.txt");
    return 0;
}

void save_music_result_angles(std::string iq_file_name, std::string music_result_angles_file_name) {
    std::vector<SamplesData> samples_data_vector;
    read_files::get_iq_samples(samples_data_vector, "data/iq_samples/" + iq_file_name);

    ProgressBar progress_bar(samples_data_vector.size());
    DoaEstimator estimator;
    DoaAngles angles;
    std::ofstream music_result_angles_file;
    music_result_angles_file.open("data/music_result_angles/" + music_result_angles_file_name);
    auto double_precision = std::numeric_limits<long double>::digits10;
    std::cout << "Saving music result angles for " << music_result_angles_file_name << "\n";
    for (auto samples_data : samples_data_vector) {
        progress_bar.update();
        angles = estimator.process_samples(samples_data, DoaTechnique::music, MusicSearch::simple_grid, M_PI / 3600);
        music_result_angles_file << std::setprecision(double_precision) << "("
                                 << angles.azimuth << "," << angles.elevation << ")\n";
    }
    music_result_angles_file.close();

    std::cout << "\nDone!\n\n";

    return;
}
