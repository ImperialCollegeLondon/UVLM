#pragma once

// This contains the necessary header files for using Eigen in UVLM

// #include <Eigen/Core>
// #include <Eigen/Dense>
// #include <Eigen/Geometry>
// No automatic parallelisation from Eigen, improves
// speed when the program is paralellised with openmp
// #define EIGEN_DONT_PARALLELIZE
#include <Eigen/Eigen>
