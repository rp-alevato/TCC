#include "libdoa/doa_estimator.h"
#include "misc/read_data_files.h"

#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

static constexpr double deg_pi = M_PI / 180;
static constexpr double pi_deg = 180 / M_PI;

int main() {
    std::vector<SamplesData> samples_data_list;
    read_files::get_iq_samples(samples_data_list, "data/iq_samples/close.txt");
    // DoaEstimator estimator;
    // DoaAngles angles;
    // GradientSpecs gradient_specs = {1e-4, 1e-6, 0.5, 0.3};
}