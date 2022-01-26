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
    float distance = 0.04;  // Distância em m
    float product;
    float angle_step = 1 * M_PI / 180;

    // Last estimated direction
    float az_index = 0;
    float el_index = 0;
    float az_esprit;  // Azimuth from ESPRIT algorithm
    float el_esprit;  // Elevation from ESPRIT algorithm

    // Antenna array output. It is the algorithm input
    Eigen::Matrix<std::complex<double>, 16, 4> x;

    // Input signal autocorrelation matrix
    Eigen::Matrix<std::complex<double>, 16, 16> Rxx;
    Eigen::Matrix<std::complex<double>, 16, 15> Qn;  // Noise eigenvectors. MUSIC related
    Eigen::Matrix<std::complex<double>, 16, 1> Qs;   // Signal eigenvector

    float x_st;
    float y_st;
    std::complex<float> phi_x;
    std::complex<float> phi_y;
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
    void initDoAEstimator(float elements_distance, int debug_flag);
    // Estimate phase rotation based on CTE reference period IQ data
    float estimate_phase_rotation(float* i_samples, float* q_samples, int reference_period_length);
    // Load antenna array output matrix. Converte de float para complexo.
    void load_x(float** i_samples, float** q_samples, int num_antennas, int num_samples, int debug_flag);
    // Compensate the phase rotation due to switching/CTE frequncy offset
    void compensateRotation(float phase_rotation, int debug_flag);
    // Estimate and decomposes covariance matrix
    void estimateRxx(int debug_flag);
    // Update the steering vector
    void updateSteeringVector(float azimuth, float elevation, int debug_flag);
    // MUSIC spectrum function. It computes the MUSIC Spatial Spectrum power in the given direction using the Steering Vecotor and Qn matrix
    float MUSICSpectrum(float azimuth, float elevation, int debug_flag);
    // Process MUSIC Spectrum
    void processMUSIC(int azimuth_max, int elevation_max, int debug_flag);
    // Initialize Selecion Matrices for ESPRIT algorithm
    void initSelectionMatrices(int debug_flag);
    // Process ESPRIT algorithm
    void processESPRIT(float channel_frequency, int debug_flag);
    // Get processed values
    void getProcessed(int algorithm, float* x, float* y);
    // Calculate wavelength based on channel central frequency
    float calculate_wavelength(float channel_frequency);
};
