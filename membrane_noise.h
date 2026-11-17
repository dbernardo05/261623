#pragma once

#include <random>

// The corrected model owns one generator and one distribution (including its
// cached normal variate) per OpenMP worker. Both belong in a restart snapshot.
struct MembraneNoise {
    std::mt19937 generator{std::random_device{}()};
    std::normal_distribution<double> distribution{0.0, 1.0};
};

inline MembraneNoise& membraneNoise() {
    static thread_local MembraneNoise noise;
    return noise;
}
