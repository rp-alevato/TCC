#include <Eigen/Dense>
#include <cmath>

// #ifndef NDEBUG
// #define NDEBUG
// #endif

// Angles in radians
struct DoaAngles {
    double azimuth;
    double elevation;
};

struct GradientOptimSpecs {
    double learning_rate;
    double accuracy;
    double diff_step;
};

enum class DoaTechnique {
    esprit,
    music
};

enum class MusicSearchOptim {
    simple_grid,
    linear_grid_gradient,
    coarse_grid_gradient,
    music_mapping,
};

class DoaEstimator {
  private:
    static constexpr int speed_of_light = 299792458;
    static constexpr int n_samples_ref = 7;
    static constexpr int n_samples = 4;
    static constexpr int n_antennas_axis = 4;
    static constexpr int n_antennas = 16;
    static constexpr int n_subvector = 12;
    static constexpr double antenna_gap_size = 0.04;  // Size in meters
    double channel_frequency;
    double phase_constant;
    Eigen::Matrix<Eigen::dcomplex, n_antennas, n_samples> samples;
    Eigen::Matrix<Eigen::dcomplex, n_antennas, n_antennas> autocorrelation_matrix;

    // MUSIC only variables
    Eigen::Matrix<Eigen::dcomplex, n_antennas, (n_antennas - 1)> noise_eigenvectors;
    Eigen::Matrix<Eigen::dcomplex, n_antennas, n_antennas> noise_eigenvectors_product;
    Eigen::Matrix<Eigen::dcomplex, n_antennas, 1> steering_vector;
    Eigen::Matrix<Eigen::dcomplex, n_antennas_axis, 1> steering_vector_x;
    Eigen::Matrix<Eigen::dcomplex, n_antennas_axis, 1> steering_vector_y;

    // ESPRIT only variables
    Eigen::Matrix<Eigen::dcomplex, n_antennas, 1> signal_eigenvector;
    Eigen::Matrix<Eigen::dcomplex, n_subvector, 1> signal_subvector_x_1;
    Eigen::Matrix<Eigen::dcomplex, n_subvector, 1> signal_subvector_x_2;
    Eigen::Matrix<Eigen::dcomplex, n_subvector, 1> signal_subvector_y_1;
    Eigen::Matrix<Eigen::dcomplex, n_subvector, 1> signal_subvector_y_2;
    static constexpr std::array<int, n_subvector> signal_subvector_x_1_index = {0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14};
    static constexpr std::array<int, n_subvector> signal_subvector_x_2_index = {1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15};
    static constexpr std::array<int, n_subvector> signal_subvector_y_1_index = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    static constexpr std::array<int, n_subvector> signal_subvector_y_2_index = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    DoaAngles process_music(MusicSearchOptim search_optmization, double grid_step, GradientOptimSpecs gradient_specs);
    double estimate_music_result(DoaAngles in_angles);
    DoaAngles music_simple_grid_search(double grid_step);
    DoaAngles music_linear_grid_gradient_search(double grid_step, GradientOptimSpecs gradient_specs);
    DoaAngles music_coarse_grid_gradient_search(double grid_step, GradientOptimSpecs gradient_specs);
    DoaAngles music_result_mapping(double grid_step);
    DoaAngles music_gradient_search(DoaAngles coarse_angles, GradientOptimSpecs gradient_specs);
    DoaAngles process_esprit();

  public:
    DoaEstimator(){};
    void load_samples(Eigen::Matrix<Eigen::dcomplex, n_antennas, n_samples>& in_samples,
                      Eigen::Matrix<Eigen::dcomplex, n_samples_ref, 1>& samples_reference,
                      double channel_frequency = 2444000000.0);
    DoaAngles process_samples(DoaTechnique technique,
                              MusicSearchOptim search_optmization = MusicSearchOptim::simple_grid,
                              double grid_step = 2 * M_PI / 1000,
                              GradientOptimSpecs gradient_specs = {0.1, 1e-5, 1e-8});
};
