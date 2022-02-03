#include "libdoa.h"

void DoaEstimator::load_samples(Eigen::Matrix<Eigen::dcomplex, n_antennas, n_samples>& in_samples,
                                Eigen::Matrix<Eigen::dcomplex, n_samples_ref, 1>& samples_reference) {

    // Calculate phase shift due to sampling with samples_ref
    double phase_difference = 0;
    for (int i = 1; i < samples_reference.size(); i++) {
        phase_difference += std::arg(samples_reference(i) * std::conj(samples_reference(i - 1)));
    }
    phase_difference /= (samples_reference.size() - 1);

    // Load samples and compensate for phase difference calculated above
    for (int i = 0; i < in_samples.rows(); i++) {
        for (int j = 0; j < in_samples.cols(); j++) {
            samples(i, j) = in_samples(i, j) * std::polar<double>(1, i * phase_difference);
        }
    }
    return;
}

DoaAngles DoaEstimator::process_samples(DoaTechnique technique, double channel_frequency = 2444000000) {
    this->autocorrelation_matrix = (this->samples * this->samples.adjoint()) / n_samples;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(this->autocorrelation_matrix);
    this->noise_eigenvectors = eigensolver.eigenvectors().block(0, 0, n_antennas, (n_antennas - 1));
    this->signal_eigenvector = eigensolver.eigenvectors().block(0, (n_antennas - 1), n_antennas, 1);

    if (technique == DoaTechnique::music) {
        return this->process_music();
    } else if (technique == DoaTechnique::esprit) {
        return this->process_esprit(channel_frequency);
    }

    return {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
}

DoaAngles DoaEstimator::process_music() {
    this->noise_eigenvectors_product = this->noise_eigenvectors * this->noise_eigenvectors.adjoint();
    return simple_search();
}

DoaAngles DoaEstimator::simple_search() {
    static constexpr double azimuth_max = 2 * M_PI;
    static constexpr double elevation_max = M_PI / 2;
    static constexpr double step = M_PI / 180;  // 1 degree / step
    DoaAngles result_doa_angles = {0, 0};
    double maximum_result = 0;

    for (double azimuth = step; azimuth < azimuth_max; azimuth += step) {
        for (double elevation = 0; elevation < elevation_max; elevation += step) {
            double result = this->estimate_music_result({azimuth, elevation});
            if (result > maximum_result) {
                maximum_result = result;
                result_doa_angles = {azimuth, elevation};
            }
        }
    }
    return result_doa_angles;
}

double DoaEstimator::estimate_music_result(DoaAngles angles) {
    Eigen::Matrix<Eigen::dcomplex, n_antennas_axis, 1> steering_vector_x;
    Eigen::Matrix<Eigen::dcomplex, n_antennas_axis, 1> steering_vector_y;
    Eigen::dcomplex music_result_complex;
    double constant = 2 * M_PI * this->antenna_gap_size;
    double phase_x = constant * std::cos(angles.azimuth) * std::sin(angles.elevation);
    double phase_y = constant * std::sin(angles.azimuth) * std::sin(angles.elevation);

    for (int i = 0; i < steering_vector_x.size(); i++) {
        steering_vector_x(i) = std::polar<double>(1.0, i * phase_x);
        steering_vector_y(i) = std::polar<double>(1.0, i * phase_y);
    }
    for (int i = 0; i < steering_vector_x.size(); i++) {
        this->steering_vector.block(steering_vector_x.size() * i, 0, steering_vector_x.size(), 1)
            = steering_vector_x * steering_vector_y(i);
    }

    music_result_complex = this->steering_vector.adjoint()
                           * this->noise_eigenvectors_product
                           * this->steering_vector;

    return (1 / music_result_complex.real());
}

DoaAngles DoaEstimator::process_esprit(double channel_frequency = 2444000000) {
    double wavelength = speed_of_light / channel_frequency;
    return {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
}
