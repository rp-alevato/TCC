#include <Eigen/Dense>
#include <cmath>
#include <vector>

// #ifndef NDEBUG
// #define NDEBUG
// #endif

#ifndef LIBDOA_H
#define LIBDOA_H

// Constants
namespace libdoa_const {
static constexpr int speed_of_light = 299792458;
static constexpr int n_samples_ref = 7;
static constexpr int n_samples = 4;
static constexpr int n_antennas_axis = 4;
static constexpr int n_antennas = 16;
static constexpr int n_subvector = 12;
static constexpr double antenna_gap_size = 0.04;  // Size in meters
}  // namespace libdoa_const

// Angles in radians
struct DoaAngles {
    double azimuth;
    double elevation;
};

struct GradientSpecs {
    double learning_rate;
    double accuracy;
    double diff_step;
    double momentum;
};

struct SamplesData {
    Eigen::Matrix<Eigen::dcomplex, libdoa_const::n_samples_ref, 1> samples_reference;
    Eigen::Matrix<Eigen::dcomplex, libdoa_const::n_antennas, libdoa_const::n_samples> samples;
    double channel_frequency;
    double rssi;
    double sl_phase_rotation;
    double sl_azimuth;
    double sl_elevation;
};

enum class DoaTechnique {
    esprit,
    music
};

enum class MusicSearch {
    simple_grid,
    coarse_grid,
    save_spectrum
};

enum class Optimization {
    none,
    finer_grid_search,
    gradient_simple,
    gradient_momentum
};

class DoaEstimator {
  public:
    // Main methods
    DoaEstimator(){};
    void load_samples(SamplesData& in_samples);
    DoaAngles process_samples(DoaTechnique technique,
                              MusicSearch search_method = MusicSearch::simple_grid,
                              double grid_step = 2 * M_PI / 1440,
                              double coarse_step = 2 * M_PI / 45,
                              Optimization optimization = Optimization::gradient_momentum,
                              GradientSpecs gradient_specs = {0.1, 1e-5, 1e-6, 0.3});

  private:
    static constexpr int speed_of_light = libdoa_const::speed_of_light;
    static constexpr int n_samples_ref = libdoa_const::n_samples_ref;
    static constexpr int n_samples = libdoa_const::n_samples;
    static constexpr int n_antennas_axis = libdoa_const::n_antennas_axis;
    static constexpr int n_antennas = libdoa_const::n_antennas;
    static constexpr int n_subvector = libdoa_const::n_subvector;
    static constexpr double antenna_gap_size = libdoa_const::antenna_gap_size;  // Size in meters
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
    DoaAngles process_esprit();
    DoaAngles process_music(MusicSearch search_method, double grid_step, double coarse_step, Optimization optimization, GradientSpecs gradient_specs);
    double estimate_music_result(DoaAngles in_angles);
    // Music search algorithms
    DoaAngles music_simple_grid_search(double grid_step);
    DoaAngles music_coarse_grid_search(double finer_step, double coarse_step, Optimization optimization, GradientSpecs gradient_specs);
    // Optimization algorithms
    DoaAngles music_finer_grid_search(DoaAngles coarse_angles, double finer_step, double coarse_step);
    DoaAngles music_gradient(DoaAngles coarse_angles, GradientSpecs gradient_specs);
    DoaAngles music_gradient_momentum(DoaAngles coarse_angles, GradientSpecs gradient_specs);
    // Miscellaneous
    void music_save_spectrum(double grid_step);
    void music_get_spectrum(std::vector<double>& rows_names, std::vector<double>& cols_names,
                            std::vector<std::vector<double>>& music_spectrum, double grid_step);
};

#endif
