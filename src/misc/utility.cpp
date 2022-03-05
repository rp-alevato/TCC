#include "misc/utility.h"

#include <cmath>
#include <limits>
#include <vector>

double utility::normalize_angle_180(double angle) {
    angle = std::fmod(angle + 180, 360);
    if (angle < 0)
        angle += 360;
    return angle - 180;
}

double utility::normalize_angle_360(double angle) {
    angle = std::fmod(angle, 360);
    if (angle < 0)
        angle += 360;
    return angle;
}

double utility::normalize_angle_pi(double angle) {
    angle = std::fmod(angle + M_PI, (2 * M_PI));
    if (angle < 0)
        angle += (2 * M_PI);
    return angle - M_PI;
}

double utility::normalize_angle_2pi(double angle) {
    angle = std::fmod(angle, (2 * M_PI));
    if (angle < 0)
        angle += (2 * M_PI);
    return angle;
}

bool utility::is_equal_double(const double x, const double y, const double epsilon) {
    return (std::abs(x - y) < epsilon);
}

bool utility::is_equal_angles(const AoaAngles x, const AoaAngles y, const double epsilon) {
    bool azimuth_equal = utility::is_equal_double(x.azimuth, y.azimuth, epsilon);
    bool elevation_equal = utility::is_equal_double(x.elevation, y.elevation, epsilon);
    return (azimuth_equal && elevation_equal);
}

double utility::angle_to_degree(const double angle) {
    return angle * 180 / M_PI;
}

double utility::angle_to_pi(const double angle) {
    return angle * M_PI / 180;
}

AoaAngles utility::angles_to_degree(AoaAngles angles) {
    angles.azimuth = angle_to_degree(angles.azimuth);
    angles.elevation = angle_to_degree(angles.elevation);
    return angles;
}

AoaAngles utility::angles_to_pi(AoaAngles angles) {
    angles.azimuth = angle_to_pi(angles.azimuth);
    angles.elevation = angle_to_pi(angles.elevation);
    return angles;
}
