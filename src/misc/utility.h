#include "libdoa/doa_estimator.h"

#ifndef MISC_UTILITY_H
#define MISC_UTILITY_H

namespace utility {
double normalize_angle_180(double angle);
double normalize_angle_360(double angle);
double normalize_angle_pi(double angle);
double normalize_angle_2pi(double angle);
bool is_equal_double(double x, double y, double epsilon = 1e-12);
bool is_equal_angles(DoaAngles x, DoaAngles y, double epsilon);
double transfom_angle_to_degree(double angle);
double transfom_angle_to_pi(double angle);
DoaAngles transfom_angles_to_degree(DoaAngles angles);
DoaAngles transfom_angles_to_pi(DoaAngles angles);
}  // namespace utility

#endif  // MISC_UTILITY_H
