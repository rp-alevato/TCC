#include <vector>

namespace stats {
double mean_double(const std::vector<double>& in_vector);
double median_sorted_double(const std::vector<double>& in_vector);
double variance_double(const std::vector<double>& in_vector, const double mean);
double std_deviation_double(const double variance);
double lowest_sorted_double(const std::vector<double>& in_vector);
double highest_sorted_double(const std::vector<double>& in_vector);
double range_double_sorted_double(const std::vector<double>& in_vector);
double mae_double(const std::vector<double>& error_vector);
double mse_double(const std::vector<double>& error_vector);
double rmse_double(const double mse);
}  // namespace stats