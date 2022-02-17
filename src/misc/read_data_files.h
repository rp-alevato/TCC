#include "libdoa/libdoa.h"

#include <fstream>
#include <iostream>

namespace read_files {
void iq_samples(std::vector<SamplesData>& samples_data_list, const std::string& file_name);
}  // namespace read_files