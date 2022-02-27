#include "misc/statistics.h"

#include <cmath>

double stats::mean_double(const std::vector<double>& in_vector) {
    double accumulator = 0;
    for (auto v : in_vector) {
        accumulator += v;
    }
    return (accumulator / in_vector.size());
}

double stats::median_sorted_double(const std::vector<double>& in_vector) {
    std::size_t size = in_vector.size();

    // Check for even case
    if (size % 2 != 0)
        return in_vector[(size - 1) / 2];

    return (in_vector[(size - 2) / 2] + in_vector[size / 2]) / 2.0;
}

double stats::mad_double(const std::vector<double>& in_vector, const double mean) {
    double accumulator = 0;
    for (auto v : in_vector) {
        accumulator += std::abs(v - mean);
    }
    return (accumulator / in_vector.size());
}

double stats::variance_double(const std::vector<double>& in_vector, const double mean) {
    double accumulator = 0;
    for (auto v : in_vector) {
        accumulator += (v - mean) * (v - mean);
    }
    return (accumulator / in_vector.size());
}

double stats::std_deviation_double(const double variance) {
    return std::sqrt(variance);
}

double stats::lowest_sorted_double(const std::vector<double>& in_vector) {
    return in_vector.front();
}

double stats::highest_sorted_double(const std::vector<double>& in_vector) {
    return in_vector.back();
}

double stats::range_double_sorted_double(const std::vector<double>& in_vector) {
    return (stats::highest_sorted_double(in_vector)
            - stats::lowest_sorted_double(in_vector));
}

double stats::mae_double(const std::vector<double>& error_vector) {
    double accumulator = 0;
    for (auto error : error_vector) {
        accumulator += std::abs(error);
    }
    return (accumulator / error_vector.size());
}

double stats::mse_double(const std::vector<double>& error_vector) {
    double accumulator = 0;
    for (auto error : error_vector) {
        accumulator += std::pow(error, 2);
    }
    return (accumulator / error_vector.size());
}

double stats::rmse_double(const double mse) {
    return std::sqrt(mse);
}
