#include "libdoa/doa_estimator.h"

#ifndef MISC_UTILITY_H
#define MISC_UTILITY_H

namespace utility {
double normalize_angle_180(double angle);
double normalize_angle_360(double angle);
double normalize_angle_pi(double angle);
double normalize_angle_2pi(double angle);
bool is_equal_double(const double x, const double y, const double epsilon = 1e-12);
bool is_equal_angles(const DoaAngles x, const DoaAngles y, const double epsilon);
double angle_to_degree(const double angle);
double angle_to_pi(const double angle);
DoaAngles angles_to_degree(DoaAngles angles);
DoaAngles angles_to_pi(DoaAngles angles);
}  // namespace utility

#endif  // MISC_UTILITY_H
