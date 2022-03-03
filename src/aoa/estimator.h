#include <Eigen/Dense>
#include <cmath>
#include <vector>

// #ifndef EIGEN_NO_DEBUG
// #define EIGEN_NO_DEBUG
// #endif

#ifndef AOA_ESTIMATOR_H
#define AOA_ESTIMATOR_H

// Constants
namespace aoa_const {
static constexpr int speed_of_light = 299792458;
static constexpr int n_samples_ref = 7;
static constexpr int n_samples = 4;
static constexpr int n_antennas_axis = 4;
static constexpr int n_antennas = 16;
static constexpr int n_subvector = 12;
static constexpr double antenna_gap_size = 0.04;  // Size in meters
}  // namespace aoa_const

struct SamplesData {
    Eigen::Matrix<Eigen::dcomplex, aoa_const::n_samples_ref, 1> samples_reference;
    Eigen::Matrix<Eigen::dcomplex, aoa_const::n_antennas, aoa_const::n_samples> samples;
    double channel_frequency;
    double rssi;
    double sl_phase_rotation;
    double sl_azimuth;
    double sl_elevation;
};

// Angles in radians
struct AoaAngles {
    double azimuth;
    double elevation;
};

struct GradientSpecs {
    double threshold;
    double diff_step;
    double learning_rate;
    double momentum;
};

enum class AoaTechnique {
    esprit,
    music
};

enum class MusicSearch {
    simple_grid,
    coarse_grid
};

enum class MusicOptimization {
    finer_grid_search,
    gradient_simple,
    gradient_adapt_lr,
    gradient_momentum,
    gradient_nesterov
};

class AoaEstimator {
  public:
    // Research purposes
    bool get_was_max_iterations() {
        return was_max_iterations;
    };
    AoaEstimator() {
        Eigen::setNbThreads(4);
    };
    // Main methods
    double estimate_music_result(const AoaAngles in_angles);  // Public for research purposes
    AoaAngles process_samples(const SamplesData& in_samples,
                              const AoaTechnique technique,
                              const MusicSearch search_method = MusicSearch::simple_grid,
                              const double grid_step = 2 * M_PI / 720,
                              const MusicOptimization optimization = MusicOptimization::gradient_simple,
                              const double coarse_step = 2 * M_PI / 60,
                              const GradientSpecs gradient_specs = {1e-5, 1e-8, 0.5, 0.95});

  private:
    bool was_max_iterations;  // Research purposes
    static constexpr int max_iterations = 1e4;
    static constexpr int speed_of_light = aoa_const::speed_of_light;
    static constexpr int n_samples_ref = aoa_const::n_samples_ref;
    static constexpr int n_samples = aoa_const::n_samples;
    static constexpr int n_antennas_axis = aoa_const::n_antennas_axis;
    static constexpr int n_antennas = aoa_const::n_antennas;
    static constexpr int n_subvector = aoa_const::n_subvector;
    static constexpr double antenna_gap_size = aoa_const::antenna_gap_size;  // Size in meters
    static constexpr std::array<int, n_subvector> signal_subvector_x_1_index = {0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14};
    static constexpr std::array<int, n_subvector> signal_subvector_x_2_index = {1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15};
    static constexpr std::array<int, n_subvector> signal_subvector_y_1_index = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    static constexpr std::array<int, n_subvector> signal_subvector_y_2_index = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    double channel_frequency;
    double phase_constant;
    Eigen::Matrix<Eigen::dcomplex, n_antennas, n_samples> samples;
    Eigen::Matrix<Eigen::dcomplex, n_antennas, n_antennas> autocorrelation_matrix;

    // Music only variables
    Eigen::Matrix<Eigen::dcomplex, n_antennas, (n_antennas - 1)> noise_eigenvectors;
    Eigen::Matrix<Eigen::dcomplex, n_antennas, n_antennas> noise_eigenvectors_product;
    Eigen::Matrix<Eigen::dcomplex, n_antennas, 1> steering_vector;
    Eigen::Matrix<Eigen::dcomplex, n_antennas_axis, 1> steering_vector_x;
    Eigen::Matrix<Eigen::dcomplex, n_antennas_axis, 1> steering_vector_y;

    // Esprit only variables
    Eigen::Matrix<Eigen::dcomplex, n_antennas, 1> signal_eigenvector;
    Eigen::Matrix<Eigen::dcomplex, n_subvector, 1> signal_subvector_x_1;
    Eigen::Matrix<Eigen::dcomplex, n_subvector, 1> signal_subvector_x_2;
    Eigen::Matrix<Eigen::dcomplex, n_subvector, 1> signal_subvector_y_1;
    Eigen::Matrix<Eigen::dcomplex, n_subvector, 1> signal_subvector_y_2;

    // Main methods
    void load_samples(const SamplesData& in_samples);
    AoaAngles process_esprit();
    AoaAngles process_music(const MusicSearch search_method, const double grid_step,
                            const MusicOptimization optimization, const double coarse_step,
                            const GradientSpecs gradient_specs);
    // double estimate_music_result(AoaAngles in_angles);
    // Music search algorithms
    AoaAngles music_simple_grid_search(const double grid_step);
    AoaAngles music_coarse_grid_search(const double finer_step, const double coarse_step, const MusicOptimization optimization, const GradientSpecs gradient_specs);
    // MusicOptimization algorithms
    AoaAngles music_finer_grid_search(const AoaAngles coarse_angles, const double finer_step, double coarse_step);
    AoaAngles music_gradient_simple(const AoaAngles coarse_angles, const GradientSpecs gradient_specs);
    AoaAngles music_gradient_simple_adapt_lr(const AoaAngles coarse_angles, const GradientSpecs gradient_specs);
    AoaAngles music_gradient_momentum(const AoaAngles coarse_angles, const GradientSpecs gradient_specs);
    AoaAngles music_gradient_nesterov(const AoaAngles coarse_angles, const GradientSpecs gradient_specs);
    // Utility
    AoaAngles shift_result_angles(AoaAngles result_angles);
    double normalize_angle_2pi(double angle);
};

#endif  // AOA_ESTIMATOR_H
