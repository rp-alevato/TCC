#include "doa/estimator.h"

#include <fstream>
#include <iostream>
#include <vector>

#ifndef MISC_READ_DATA_FILES_H
#define MISC_READ_DATA_FILES_H

namespace read_files {
void get_iq_samples(std::vector<SamplesData>& samples_data_list, const std::string& file_name);
void get_music_result_angles(std::vector<DoaAngles>& music_results, const std::string& file_name);
}  // namespace read_files

#endif  // MISC_READ_DATA_FILES_H
