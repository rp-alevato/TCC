#include "read_data_files.h"

#include <stdexcept>

int get_iq_samples_data_block(SamplesData& samples_data, std::ifstream& iq_file);
int get_iq_samples_header(SamplesData& samples_data, std::ifstream& iq_file);
int get_iq_samples_reference(SamplesData& samples_data, std::ifstream& iq_file);
int get_iq_samples_main(SamplesData& samples_data, std::ifstream& iq_file);
int string_to_complex(Eigen::dcomplex& complex, const std::string& complex_string);

void read_files::iq_samples(std::vector<SamplesData>& samples_data_list, const std::string& file_name) {
    std::ifstream iq_file;
    SamplesData current_data;

    iq_file.open(file_name);
    if (!iq_file.is_open()) {
        throw std::runtime_error("Error reading file");
    }

    while (get_iq_samples_data_block(current_data, iq_file) >= 0) {
        samples_data_list.push_back(current_data);
    }

    return;
}

int get_iq_samples_data_block(SamplesData& samples_data, std::ifstream& iq_file) {
    std::string line;

    if (get_iq_samples_header(samples_data, iq_file) < 0) {
        return -1;
    }
    if (get_iq_samples_reference(samples_data, iq_file) < 0) {
        return -1;
    }
    if (get_iq_samples_main(samples_data, iq_file) < 0) {
        return -1;
    }

    std::getline(iq_file, line);  // Throw away newline at the end of samples data block

    return 0;
}

int get_iq_samples_header(SamplesData& samples_data, std::ifstream& iq_file) {
    std::string line;

    // Get channel_frequency
    std::getline(iq_file, line);
    if (iq_file.eof()) {
        return -1;
    }
    samples_data.channel_frequency = std::stod(line.substr(19));

    // Get rssi
    std::getline(iq_file, line);
    if (iq_file.eof()) {
        return -1;
    }
    samples_data.rssi = std::stod(line.substr(6));

    // Get phase_rotation
    std::getline(iq_file, line);
    if (iq_file.eof()) {
        return -1;
    }
    samples_data.sl_phase_rotation = std::stod(line.substr(16));

    // Get azimuth and elevation
    std::getline(iq_file, line);
    if (iq_file.eof()) {
        return -1;
    }
    Eigen::dcomplex azimuth_elevation;
    if (string_to_complex(azimuth_elevation, line.substr(20)) < 0) {
        return -1;
    }
    samples_data.sl_azimuth = azimuth_elevation.real();
    samples_data.sl_elevation = azimuth_elevation.imag();

    return 0;
}

int get_iq_samples_reference(SamplesData& samples_data, std::ifstream& iq_file) {
    std::string line;

    std::getline(iq_file, line);  // Throw away line with "reference samples:"

    for (auto i = 0; i < libdoa_const::n_samples_ref; i++) {
        std::getline(iq_file, line);
        if (iq_file.eof()) {
            return -1;
        }
        Eigen::dcomplex complex;
        if (string_to_complex(complex, line) < 0) {
            return -1;
        }
        samples_data.samples_reference(i) = complex;
    }

    return 0;
}

int get_iq_samples_main(SamplesData& samples_data, std::ifstream& iq_file) {
    std::string line;

    std::getline(iq_file, line);  // Throw away line with "samples:"

    for (auto i = 0; i < libdoa_const::n_antennas; i++) {
        std::getline(iq_file, line);
        if (iq_file.eof()) {
            return -1;
        }
        std::string::size_type begin_position = 0;
        std::string::size_type end_position = 0;
        for (auto j = 0; j < libdoa_const::n_samples; j++) {
            end_position = line.find(")", begin_position);
            if (end_position == std::string::npos) {
                return -1;
            }
            auto complex_string = line.substr(begin_position, end_position);
            Eigen::dcomplex complex;
            if (string_to_complex(complex, complex_string) < 0) {
                return -1;
            }
            samples_data.samples(i, j) = complex;
            begin_position = end_position + 3;
        }
    }

    return 0;
}

int string_to_complex(Eigen::dcomplex& complex, const std::string& complex_string) {
    std::string::size_type comma_position;

    comma_position = complex_string.find(',');
    if (comma_position == std::string::npos) {
        return -1;
    }
    complex = {std::stod(complex_string.substr(1, comma_position)),
               std::stod(complex_string.substr(comma_position + 1))};
    return 0;
}
