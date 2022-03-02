#include "estimator.h"

#include <chrono>
#include <cmath>
#include <iostream>

// Tests which is faster: for loop with a double increment or for loop with
// int increment with a multiplication after
void for_double_int_mult_test() {
    long int int_max = 1e7;
    double double_step = 1e-2;
    double double_max = int_max * double_step;
    volatile long int int_it = 0;
    volatile long int double_it = 0;

    std::cout << "Starting test double vs int for loop for 1e7 iteration\n";
    auto double_t1 = std::chrono::high_resolution_clock::now();
    for (double i = 0; i < double_max; i += double_step) {
        double_it++;
    }
    auto double_t2 = std::chrono::high_resolution_clock::now();

    auto int_t1 = std::chrono::high_resolution_clock::now();
    for (long int i = 0; i < int_max; i++) {
        volatile double temp = i * double_step;
        int_it++;
    }
    auto int_t2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<float, std::milli> runtime_int = int_t2 - int_t1;
    std::chrono::duration<float, std::milli> runtime_double = double_t2 - double_t1;

    std::cout << "Runtime Int: " << runtime_int.count() << "ms\n";
    std::cout << "Iterations: " << int_it << "\n";
    std::cout << "Runtime Double: " << runtime_double.count() << "ms\n";
    std::cout << "Iterations: " << double_it << "\n\n";
}

// Tests which is faster: for loop with a float increment or for loop with
// int increment with a multiplication after
void for_float_int_mult_test() {
    int int_max = 1e7;
    float float_step = 1e-2;
    float float_max = int_max * float_step;
    volatile long int int_it = 0;
    volatile long int float_it = 0;

    std::cout << "Starting test float vs int for loop for 1e7 iteration\n";
    auto float_t1 = std::chrono::high_resolution_clock::now();
    for (float i = 0.0; i < float_max; i += float_step) {
        float_it++;
    }
    auto float_t2 = std::chrono::high_resolution_clock::now();

    auto int_t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < int_max; i++) {
        volatile float temp = i * float_step;
        int_it++;
    }
    auto int_t2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<float, std::milli> runtime_int = int_t2 - int_t1;
    std::chrono::duration<float, std::milli> runtime_float = float_t2 - float_t1;

    std::cout << "Runtime Int: " << runtime_int.count() << "ms\n";
    std::cout << "Iterations: " << int_it << "\n";
    std::cout << "Runtime Float: " << runtime_float.count() << "ms\n";
    std::cout << "Iterations: " << float_it << "\n\n";
}

int main() {
    for_double_int_mult_test();
    for_float_int_mult_test();
}