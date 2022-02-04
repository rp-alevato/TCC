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

enum class DoaTechnique {
    esprit = 0,
    music
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

  public:
    double estimate_music_result(DoaAngles in_angles);
    DoaAngles simple_search();
    DoaAngles process_music();
    DoaAngles process_esprit();

    DoaEstimator(){};

    void load_samples(Eigen::Matrix<Eigen::dcomplex, n_antennas, n_samples>& in_samples,
                      Eigen::Matrix<Eigen::dcomplex, n_samples_ref, 1>& samples_reference,
                      double channel_frequency = 2444000000.0);

    DoaAngles process_samples(DoaTechnique technique);
};
