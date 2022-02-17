#include "libdoa.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

// ****************************************************************************
// ******************************  MAIN METHODS *******************************
// ****************************************************************************
void DoaEstimator::load_samples(SamplesData& samples_data) {

    this->channel_frequency = samples_data.channel_frequency;

    // Calculate phase shift due to sampling with samples_ref
    double phase_difference = 0;
    for (int i = 1; i < n_samples_ref; i++) {
        double curr_argument = std::arg(samples_data.samples_reference(i));
        double prev_argument = std::arg(samples_data.samples_reference(i - 1));
        double curr_phase_diff = curr_argument - prev_argument;
        if (curr_phase_diff > M_PI) {
            curr_phase_diff -= (2 * M_PI);
        } else if (curr_phase_diff < (-M_PI)) {
            curr_phase_diff += (2 * M_PI);
        }
        phase_difference += curr_phase_diff;
    }
    phase_difference /= (n_samples_ref - 1);
    phase_difference *= -2;

    // Load samples and compensate for phase difference calculated above
    for (int i = 0; i < samples_data.samples.rows(); i++) {
        for (int j = 0; j < samples_data.samples.cols(); j++) {
            samples(i, j) = samples_data.samples(i, j) * std::polar<double>(1, i * phase_difference);
        }
    }
    return;
}

DoaAngles DoaEstimator::process_samples(DoaTechnique technique,
                                        MusicSearch search_method,
                                        double grid_step,
                                        double coarse_step,
                                        Optimization optimization,
                                        GradientSpecs gradient_specs) {

    this->autocorrelation_matrix = (this->samples * this->samples.adjoint()) / n_samples;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(this->autocorrelation_matrix);
    this->noise_eigenvectors = eigensolver.eigenvectors().block(0, 0, n_antennas, (n_antennas - 1));
    this->signal_eigenvector = eigensolver.eigenvectors().block(0, (n_antennas - 1), n_antennas, 1);
    this->phase_constant = (2 * M_PI * this->antenna_gap_size * this->channel_frequency) / speed_of_light;

    if (technique == DoaTechnique::music) {
        return this->process_music(search_method, grid_step, coarse_step, optimization, gradient_specs);
    } else if (technique == DoaTechnique::esprit) {
        return this->process_esprit();
    }

    return {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
}

DoaAngles DoaEstimator::process_esprit() {
    double phase_x, phase_y;
    DoaAngles angles = {0, 0};

    this->signal_subvector_x_1 = this->signal_eigenvector(this->signal_subvector_x_1_index);
    this->signal_subvector_x_2 = this->signal_eigenvector(this->signal_subvector_x_2_index);
    this->signal_subvector_y_1 = this->signal_eigenvector(this->signal_subvector_y_1_index);
    this->signal_subvector_y_2 = this->signal_eigenvector(this->signal_subvector_y_2_index);

    phase_x = std::arg(((this->signal_subvector_x_1.adjoint() * this->signal_subvector_x_1).inverse()
                        * (this->signal_subvector_x_1.adjoint() * this->signal_subvector_x_2))(0));
    phase_y = std::arg(((this->signal_subvector_y_1.adjoint() * this->signal_subvector_y_1).inverse()
                        * (this->signal_subvector_y_1.adjoint() * this->signal_subvector_y_2))(0));

    angles.azimuth = std::atan2(phase_y, phase_x);

    // Sometimes if the value is close to 90º it can fail the asin(x) (x > 1 || x < -1)
    // To correct for this we take the real part of the complex asin
    angles.elevation = std::real(std::asin(std::complex(phase_y / (this->phase_constant * std::sin(angles.azimuth)))));

    return angles;
}

DoaAngles DoaEstimator::process_music(MusicSearch search_method, double grid_step, double coarse_step,
                                      Optimization optimization, GradientSpecs gradient_specs) {
    DoaAngles angles;
    this->noise_eigenvectors_product = this->noise_eigenvectors * this->noise_eigenvectors.adjoint();
    switch (search_method) {
        case MusicSearch::simple_grid:
            if (grid_step < 1e-5) {
                throw std::invalid_argument("Grid step too small");
            }
            angles = music_simple_grid_search(grid_step);
            break;
        case MusicSearch::coarse_grid:
            angles = music_coarse_grid_search(grid_step, coarse_step, optimization, gradient_specs);
            break;
        case MusicSearch::save_spectrum:
            music_save_spectrum(grid_step);
            break;
        default:
            angles = {0, 0};
    }

    return angles;
}

double DoaEstimator::estimate_music_result(DoaAngles angles) {
    Eigen::dcomplex music_result_complex;
    double phase_x = this->phase_constant * std::cos(angles.azimuth) * std::sin(angles.elevation);
    double phase_y = this->phase_constant * std::sin(angles.azimuth) * std::sin(angles.elevation);

    for (int i = 0; i < n_antennas_axis; i++) {
        this->steering_vector_x(i) = std::polar<double>(1.0, i * phase_x);
        this->steering_vector_y(i) = std::polar<double>(1.0, i * phase_y);
    }
    for (int i = 0; i < n_antennas_axis; i++) {
        this->steering_vector.block(n_antennas_axis * i, 0, n_antennas_axis, 1)
            = this->steering_vector_x * this->steering_vector_y(i);
    }

    music_result_complex = this->steering_vector.adjoint()
                           * this->noise_eigenvectors_product
                           * this->steering_vector;

    return (1 / music_result_complex.real());
}

// ****************************************************************************
// *************************  MUSIC SEARCH ALGORITHMS *************************
// ****************************************************************************
DoaAngles DoaEstimator::music_simple_grid_search(double grid_step) {
    int azimuth_max_iter = (int)(2 * M_PI / grid_step);
    int elevation_max_iter = (int)((M_PI / 2) / grid_step);
    DoaAngles result_angles = {0, 0};
    double maximum_result = 0;

    for (int azimuth_index = 0; azimuth_index < azimuth_max_iter; azimuth_index++) {
        double azimuth = azimuth_index * grid_step;
        for (int elevation_index = 0; elevation_index < elevation_max_iter; elevation_index++) {
            double elevation = elevation_index * grid_step;
            double result = this->estimate_music_result({azimuth, elevation});
            if (result > maximum_result) {
                maximum_result = result;
                result_angles = {azimuth, elevation};
            }
        }
    }
    return result_angles;
}

DoaAngles DoaEstimator::music_coarse_grid_search(double finer_step, double coarse_step,
                                                 Optimization optimization, GradientSpecs gradient_specs) {
    int azimuth_max_iter = (int)(2 * M_PI / coarse_step);
    int elevation_max_iter = (int)((M_PI / 2) / coarse_step);
    DoaAngles coarse_angles = {0, 0};
    double maximum_result = 0;

    int counter = 0;

    for (int azimuth_index = 0; azimuth_index < azimuth_max_iter; azimuth_index++) {
        double azimuth = azimuth_index * coarse_step;
        for (int elevation_index = 0; elevation_index < elevation_max_iter; elevation_index++) {
            counter++;
            double elevation = elevation_index * coarse_step;
            double result = this->estimate_music_result({azimuth, elevation});
            if (result > maximum_result) {
                maximum_result = result;
                coarse_angles = {azimuth, elevation};
            }
        }
    }

    switch (optimization) {
        case Optimization::finer_grid_search:
            if (finer_step < 1e-5) {
                throw std::invalid_argument("Grid step too small");
            }
            return music_finer_grid_search(coarse_angles, finer_step, coarse_step);
            break;
        case Optimization::gradient_simple:
            return music_gradient(coarse_angles, gradient_specs);
            break;
        case Optimization::gradient_momentum:
            return music_gradient_momentum(coarse_angles, gradient_specs);
            break;
        default:
            break;
    }

    return coarse_angles;
}

// ****************************************************************************
// *************************  OPTIMIZATION ALGORITHMS *************************
// ****************************************************************************
DoaAngles DoaEstimator::music_finer_grid_search(DoaAngles coarse_angles, double finer_step, double coarse_step) {
    int azimuth_init_iter = (int)((coarse_angles.azimuth - coarse_step) / finer_step);
    int elevation_init_iter = (int)((coarse_angles.elevation - coarse_step) / finer_step);
    int azimuth_max_iter = (int)((coarse_angles.azimuth + coarse_step) / finer_step);
    int elevation_max_iter = (int)((coarse_angles.elevation + coarse_step) / finer_step);
    DoaAngles result_angles = {0, 0};
    double maximum_result = 0;

    for (int azimuth_index = azimuth_init_iter; azimuth_index < azimuth_max_iter; azimuth_index++) {
        double azimuth = azimuth_index * finer_step;
        for (int elevation_index = elevation_init_iter; elevation_index < elevation_max_iter; elevation_index++) {
            double elevation = elevation_index * finer_step;
            double result = this->estimate_music_result({azimuth, elevation});
            if (result > maximum_result) {
                maximum_result = result;
                result_angles = {azimuth, elevation};
            }
        }
    }
    return result_angles;
}

DoaAngles DoaEstimator::music_gradient(DoaAngles coarse_angles, GradientSpecs gradient_specs) {
    double learning_rate = gradient_specs.learning_rate;
    double accuracy = gradient_specs.accuracy;
    double diff_step = gradient_specs.diff_step;
    bool continue_azimuth = true;
    bool continue_elevation = true;
    DoaAngles result_angles = coarse_angles;
    int iter = 0;
    do {
        iter++;
        double gradient_azimuth, gradient_elevation;
        double curr_azimuth = result_angles.azimuth;
        double music_result = estimate_music_result(result_angles);
        gradient_azimuth = estimate_music_result({curr_azimuth + diff_step, result_angles.elevation});
        gradient_azimuth = (gradient_azimuth - music_result) / diff_step;
        double step = learning_rate * gradient_azimuth;
        result_angles.azimuth += step;
        continue_azimuth = std::abs(step) > accuracy;
        gradient_elevation = estimate_music_result({curr_azimuth, result_angles.elevation + diff_step});
        gradient_elevation = (gradient_elevation - music_result) / diff_step;
        step = learning_rate * gradient_elevation;
        result_angles.elevation += step;
        result_angles.elevation = std::abs(result_angles.elevation);
        continue_elevation = std::abs(step) > accuracy;
    } while (continue_azimuth && continue_elevation && iter < 100);

    std::cout << "iter: " << iter << "\n";

    return result_angles;
}

DoaAngles DoaEstimator::music_gradient_momentum(DoaAngles coarse_angles, GradientSpecs gradient_specs) {
    double learning_rate = gradient_specs.learning_rate;
    double accuracy = gradient_specs.accuracy;
    double diff_step = gradient_specs.diff_step;
    double momentum = gradient_specs.momentum;
    bool continue_azimuth = true;
    bool continue_elevation = true;
    DoaAngles result_angles = coarse_angles;
    double prev_step_azimuth = 0;
    double prev_step_elevation = 0;
    int iter = 0;

    do {
        iter++;
        double gradient_azimuth, gradient_elevation;
        double curr_azimuth = result_angles.azimuth;
        double music_result = estimate_music_result(result_angles);
        gradient_azimuth = estimate_music_result({curr_azimuth + diff_step, result_angles.elevation});
        gradient_azimuth = (gradient_azimuth - music_result) / diff_step;
        double step = learning_rate * gradient_azimuth + momentum * prev_step_azimuth;
        result_angles.azimuth += step;
        continue_azimuth = std::abs(step) > accuracy;
        prev_step_azimuth = step;
        gradient_elevation = estimate_music_result({curr_azimuth, result_angles.elevation + diff_step});
        gradient_elevation = (gradient_elevation - music_result) / diff_step;
        step = learning_rate * gradient_elevation + momentum * prev_step_elevation;
        result_angles.elevation += step;
        result_angles.elevation = std::abs(result_angles.elevation);
        continue_elevation = std::abs(step) > accuracy;
        prev_step_elevation = step;
        // std::cout << result_angles.azimuth << ", " << result_angles.elevation << "\n";
    } while (continue_azimuth && continue_elevation && iter < 100);

    // std::cout << "Iterations: " << iter << "\n";

    return result_angles;
}

// ****************************************************************************
// ******************************  MISCELLANEOUS ******************************
// ****************************************************************************
void DoaEstimator::music_save_spectrum(double const grid_step) {

    // Data structure to export a csv
    std::vector<double> rows_names;
    std::vector<double> cols_names;
    std::vector<std::vector<double>> music_spectrum;
    static int file_index = 1;

    this->music_get_spectrum(rows_names, cols_names, music_spectrum, grid_step);

    // Print to file
    std::ofstream spectrum_csv;
    std::stringstream file_index_ss;
    file_index_ss << std::setw(2) << std::setfill('0') << file_index;
    spectrum_csv.open("data/music_spectrum_data/spectrum_" + file_index_ss.str() + ".csv");
    auto double_precision = std::numeric_limits<long double>::digits10 + 1;

    for (auto it = cols_names.begin(); it != cols_names.end(); it++) {
        spectrum_csv << std::setprecision(double_precision) << *it;
        if (std::next(it) != cols_names.end()) {
            spectrum_csv << ",";
        }
    }
    spectrum_csv << "\n";

    for (unsigned int i = 0; i < rows_names.size(); i++) {
        spectrum_csv << std::setprecision(double_precision) << rows_names[i] << ",";
        for (unsigned int j = 0; j < music_spectrum[0].size(); j++) {
            spectrum_csv << std::setprecision(double_precision) << music_spectrum[i][j];
            if ((j + 1) < music_spectrum[0].size()) {
                spectrum_csv << ",";
            }
        }
        spectrum_csv << "\n";
    }
    spectrum_csv.close();
    file_index++;

    return;
}

void DoaEstimator::music_get_spectrum(std::vector<double>& rows_names, std::vector<double>& cols_names,
                                      std::vector<std::vector<double>>& music_spectrum, const double grid_step) {

    int azimuth_max_iter = (int)(2 * M_PI / grid_step);
    int elevation_max_iter = (int)((M_PI / 2) / grid_step);

    for (int elevation_index = 0; elevation_index < elevation_max_iter; elevation_index++) {
        cols_names.push_back(elevation_index * grid_step);
    }
    for (int azimuth_index = 0; azimuth_index < azimuth_max_iter; azimuth_index++) {
        double azimuth = azimuth_index * grid_step;
        rows_names.push_back(azimuth);
        music_spectrum.push_back({});
        for (int elevation_index = 0; elevation_index < elevation_max_iter; elevation_index++) {
            double elevation = std::abs(elevation_index * grid_step);
            double result = this->estimate_music_result({azimuth, elevation});
            music_spectrum[azimuth_index].push_back(result);
        }
    }

    return;
}
