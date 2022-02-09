#include "libdoa.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

void DoaEstimator::load_samples(Eigen::Matrix<Eigen::dcomplex, n_antennas, n_samples>& in_samples,
                                Eigen::Matrix<Eigen::dcomplex, n_samples_ref, 1>& samples_reference,
                                double channel_frequency) {

    this->channel_frequency = channel_frequency;

    // Calculate phase shift due to sampling with samples_ref
    double phase_difference = 0;
    for (int i = 1; i < n_samples_ref; i++) {
        phase_difference += std::arg(samples_reference(i) * std::conj(samples_reference(i - 1)));
    }
    phase_difference /= (n_samples_ref - 1);

    // Load samples and compensate for phase difference calculated above
    for (int i = 0; i < in_samples.rows(); i++) {
        for (int j = 0; j < in_samples.cols(); j++) {
            samples(i, j) = in_samples(i, j) * std::polar<double>(1, i * phase_difference);
        }
    }
    return;
}

DoaAngles DoaEstimator::process_samples(DoaTechnique technique,
                                        MusicSearchOptim search_optmization,
                                        double grid_step,
                                        GradientOptimSpecs gradient_specs) {
    this->autocorrelation_matrix = (this->samples * this->samples.adjoint()) / n_samples;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(this->autocorrelation_matrix);
    this->noise_eigenvectors = eigensolver.eigenvectors().block(0, 0, n_antennas, (n_antennas - 1));
    this->signal_eigenvector = eigensolver.eigenvectors().block(0, (n_antennas - 1), n_antennas, 1);
    this->phase_constant = (2 * M_PI * this->antenna_gap_size * this->channel_frequency)
                           / this->speed_of_light;

    if (technique == DoaTechnique::music) {
        return this->process_music(search_optmization, grid_step, gradient_specs);
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

    angles.elevation = std::asin(phase_y / (this->phase_constant * std::sin(angles.azimuth)));

    return angles;
}

DoaAngles DoaEstimator::process_music(MusicSearchOptim search_optmization, double grid_step, GradientOptimSpecs gradient_specs) {
    DoaAngles angles;
    this->noise_eigenvectors_product = this->noise_eigenvectors * this->noise_eigenvectors.adjoint();
    switch (search_optmization) {
        case MusicSearchOptim::simple_grid:
            angles = music_simple_grid_search(grid_step);
            break;
        case MusicSearchOptim::linear_grid_gradient:
            angles = music_linear_grid_gradient_search(grid_step, gradient_specs);
            break;
        case MusicSearchOptim::coarse_grid_gradient:
            angles = music_coarse_grid_gradient_search(grid_step, gradient_specs);
            break;
        case MusicSearchOptim::music_mapping:
            angles = music_result_mapping(grid_step);
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

DoaAngles DoaEstimator::music_simple_grid_search(double grid_step) {
    static constexpr double azimuth_max = 2 * M_PI;
    static constexpr double elevation_max = M_PI / 2;
    DoaAngles result_angles = {0, 0};
    double maximum_result = 0;

    for (double azimuth = 0; azimuth < azimuth_max; azimuth += grid_step) {
        for (double elevation = 0; elevation < elevation_max; elevation += grid_step) {
            double result = this->estimate_music_result({azimuth, elevation});
            if (result > maximum_result) {
                maximum_result = result;
                result_angles = {azimuth, elevation};
            }
        }
    }
    return result_angles;
}

DoaAngles DoaEstimator::music_linear_grid_gradient_search(double grid_step, GradientOptimSpecs gradient_specs) {
    double azimuth_max = 2 * M_PI;
    double elevation_max = M_PI / 2;
    DoaAngles coarse_angles = {0, 0};
    double maximum_result = 0;
    double elevation, azimuth;

    elevation = elevation_max / 2;

    for (azimuth = 0; azimuth < azimuth_max; azimuth += grid_step) {
        double result = this->estimate_music_result({azimuth, elevation});
        if (result > maximum_result) {
            maximum_result = result;
            coarse_angles = {azimuth, elevation};
        }
    }

    maximum_result = 0;
    azimuth = coarse_angles.azimuth;
    for (elevation = 0; elevation < elevation_max; elevation += grid_step) {
        double result = this->estimate_music_result({azimuth, elevation});
        if (result > maximum_result) {
            maximum_result = result;
            coarse_angles = {azimuth, elevation};
        }
    }

    return music_gradient_search(coarse_angles, gradient_specs);
}

DoaAngles DoaEstimator::music_coarse_grid_gradient_search(double grid_step, GradientOptimSpecs gradient_specs) {
    static constexpr double azimuth_max = 2 * M_PI;
    static constexpr double elevation_max = M_PI / 2;
    DoaAngles coarse_angles = {0, 0};
    double maximum_result = 0;

    for (double azimuth = 0; azimuth < azimuth_max; azimuth += grid_step) {
        for (double elevation = 0; elevation < elevation_max; elevation += grid_step) {
            double result = this->estimate_music_result({azimuth, elevation});
            if (result > maximum_result) {
                maximum_result = result;
                coarse_angles = {azimuth, elevation};
            }
        }
    }

    return music_gradient_search(coarse_angles, gradient_specs);
}

// TODO Calculate derivative with: (f(x + h) - f(x - h)) / 2 * h. This will increase the calculation
// of estimate_music_result from 3 to 4 on each iteration. Check if it's worth it.
DoaAngles DoaEstimator::music_gradient_search(DoaAngles coarse_angles, GradientOptimSpecs gradient_specs) {
    double learning_rate = gradient_specs.learning_rate;
    double accuracy = gradient_specs.accuracy;
    double diff_step = gradient_specs.diff_step;
    bool continue_azimuth = true;
    bool continue_elevation = true;
    DoaAngles result_angles = coarse_angles;
    int counter = 0;
    do {
        counter++;
        double gradient_azimuth, gradient_elevation;
        double curr_azimuth = result_angles.azimuth;
        double music_result = estimate_music_result(result_angles);
        if (continue_azimuth) {
            gradient_azimuth = estimate_music_result({curr_azimuth + diff_step, result_angles.elevation});
            gradient_azimuth = (gradient_azimuth - music_result) / diff_step;
            double step = learning_rate * gradient_azimuth;
            result_angles.azimuth += step;
            continue_azimuth = std::abs(step) > accuracy;
        }
        if (continue_elevation) {
            gradient_elevation = estimate_music_result({curr_azimuth, result_angles.elevation + diff_step});
            gradient_elevation = (gradient_elevation - music_result) / diff_step;
            double step = learning_rate * gradient_elevation;
            result_angles.elevation += step;
            continue_elevation = std::abs(step) > accuracy;
        }
    } while ((continue_azimuth || continue_elevation) && counter < 1000);

    return result_angles;
}

DoaAngles DoaEstimator::music_result_mapping(double grid_step) {
    static constexpr double azimuth_max = 2 * M_PI;
    static constexpr double elevation_max = M_PI / 2;
    DoaAngles result_angles = {0, 0};
    double maximum_result = 0;

    // Data structure to export a csv
    std::vector<double> csv_rows_names;
    std::vector<double> csv_cols_names;
    std::vector<std::vector<double>> csv_values;
    int row_index = 0;

    for (double elevation = 0; elevation < elevation_max; elevation += grid_step) {
        csv_cols_names.push_back(elevation);
    }
    for (double azimuth = 0; azimuth < azimuth_max; azimuth += grid_step) {
        csv_rows_names.push_back(azimuth);
        csv_values.push_back({});
        for (double elevation = 0; elevation < elevation_max; elevation += grid_step) {
            double result = this->estimate_music_result({azimuth, elevation});
            if (result > maximum_result) {
                maximum_result = result;
                result_angles = {azimuth, elevation};
            }
            csv_values[row_index].push_back(result);
        }
        row_index++;
    }

    std::ofstream values_csv;
    values_csv.open("values.csv");
    auto double_precision = std::numeric_limits<long double>::digits10;

    for (auto it = csv_cols_names.begin(); it != csv_cols_names.end(); it++) {
        values_csv << std::fixed << std::setprecision(double_precision) << *it;
        if (std::next(it) != csv_cols_names.end()) {
            values_csv << ",";
        }
    }
    values_csv << "\n";

    for (unsigned int i = 0; i < csv_rows_names.size(); i++) {
        values_csv << std::fixed << std::setprecision(double_precision) << csv_rows_names[i] << ",";
        for (unsigned int j = 0; j < csv_values[0].size(); j++) {
            values_csv << std::fixed << std::setprecision(double_precision) << csv_values[i][j];
            if ((j + 1) < csv_values[0].size()) {
                values_csv << ",";
            }
        }
        values_csv << "\n";
    }
    values_csv.close();

    return result_angles;
}
