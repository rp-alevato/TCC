#include <Eigen/Dense>
#include <math.h>
#include <stdlib.h>

#define NUM_ANTENNAS 16
#define IQLENGTH     4  // The Number of samples used to estimate the DoA
#define MAX_SAMPSIZE 30

class aoa_estimator {

  private:
    // Signal related
    int num_samples = IQLENGTH;

    // Algorithm constants
    int antenna_num = 16;
    double distance = 0.04;  // Distância em m
    double product;
    double angle_step = 1 * M_PI / 180;  // 1º in radians

    // Last estimated direction
    double az_index = 0;
    double el_index = 0;
    double az_esprit;  // Azimuth from ESPRIT algorithm
    double el_esprit;  // Elevation from ESPRIT algorithm

    // Antenna array output. It is the algorithm input
    Eigen::Matrix<std::complex<double>, 16, 4> x;

    // Input signal autocorrelation matrix
    Eigen::Matrix<std::complex<double>, 16, 16> Rxx;
    Eigen::Matrix<std::complex<double>, 16, 15> Qn;  // Noise eigenvectors. MUSIC related
    Eigen::Matrix<std::complex<double>, 16, 1> Qs;   // Signal eigenvector

    double x_st;
    double y_st;
    std::complex<double> phi_x;
    std::complex<double> phi_y;
    Eigen::Matrix<std::complex<double>, 1, 4> a_x;       // steering vector x component
    Eigen::Matrix<std::complex<double>, 1, 4> a_y;       // steering vector y component
    Eigen::Matrix<std::complex<double>, 16, 1> st_vec;   // steering vector
    Eigen::Matrix<std::complex<double>, 16, 16> Qn_Prd;  // Qn times Qn.adjoint(). It should be stored before executing the MUSIC spectrum scanning
    Eigen::Matrix<std::complex<double>, 1, 1> aux;       // Used for the MUSIC spectrum denominator

    // ESPRIT selection matrices
    // Equações 2.45 e 2.46
    Eigen::Matrix<std::complex<double>, 12, 16> Js_1y;
    Eigen::Matrix<std::complex<double>, 12, 16> Js_2y;
    Eigen::Matrix<std::complex<double>, 12, 16> Js_1x;
    Eigen::Matrix<std::complex<double>, 12, 16> Js_2x;
    // ESPRIT signal matrices
    Eigen::Matrix<std::complex<double>, 12, 1> Qs_1x;
    Eigen::Matrix<std::complex<double>, 12, 1> Qs_2x;
    Eigen::Matrix<std::complex<double>, 12, 1> Qs_1y;
    Eigen::Matrix<std::complex<double>, 12, 1> Qs_2y;
    // Auxiliary variables for ESPRIT algorithm
    Eigen::Matrix<std::complex<double>, 1, 1> phi_X_LS;
    Eigen::Matrix<std::complex<double>, 1, 1> phi_Y_LS;

  public:
    // Init DoA estimator
    void initDoAEstimator(double elements_distance);
    // Estimate phase rotation based on CTE reference period IQ data
    double estimate_phase_rotation(double* i_samples, double* q_samples, int reference_period_length);
    // Load antenna array output matrix. Converte de double para complexo.
    void load_x(double** i_samples, double** q_samples, int num_antennas, int num_samples);
    // Compensate the phase rotation due to switching/CTE frequncy offset
    void compensateRotation(double phase_rotation);
    // Estimate and decomposes covariance matrix
    void estimateRxx();
    // Update the steering vector
    void updateSteeringVector(double azimuth, double elevation, double channel_frequency);
    // MUSIC spectrum function. It computes the MUSIC Spatial Spectrum power in the given direction using the Steering Vecotor and Qn matrix
    double MUSICSpectrum(double azimuth, double elevation, double channel_frequency);
    // Process MUSIC Spectrum
    void processMUSIC(int azimuth_max, int elevation_max, double channel_frequency);
    // Initialize Selecion Matrices for ESPRIT algorithm
    void initSelectionMatrices();
    // Process ESPRIT algorithm
    void processESPRIT(double channel_frequency);
    // Get processed values
    void getProcessed(int algorithm, double* x, double* y);
    // Calculate wavelength based on channel central frequency
    double calculate_wavelength(double channel_frequency);
};
