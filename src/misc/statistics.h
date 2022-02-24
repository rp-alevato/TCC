#include <vector>

namespace stats {
double mean_double(std::vector<double>& in_vector);
double median_sorted_double(std::vector<double>& in_vector);
double variance_double(std::vector<double>& in_vector, double mean);
double std_deviation_double(double variance);
double lowest_sorted_double(std::vector<double>& in_vector);
double highest_sorted_double(std::vector<double>& in_vector);
double range_double_sorted_double(std::vector<double>& in_vector);
double mae_double(std::vector<double>& error_vector);
double mse_double(std::vector<double>& error_vector);
double rmse_double(double mse);
}  // namespace stats