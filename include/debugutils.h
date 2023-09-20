#pragma once

#include <string>
#include <fstream>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <Eigen/Dense>

/**
 * @file debug_utils.h
 * @brief Contains utility functions for common debugging.
 */

/**
 * @brief Check if all elements in a matrix are finite.
 *
 * This function checks if all elements in a given Eigen matrix are finite.
 *
 * @param x The Eigen matrix to check for finiteness.
 * @return True if all elements are finite, otherwise false.
 */
template<typename Derived>
inline bool is_finite(const Eigen::MatrixBase<Derived>& x)
{
    return ((x - x).array() == (x - x).array()).all();
}

/**
 * @brief Check if all elements in a matrix are NaN (Not-a-Number).
 *
 * This function checks if all elements in a given Eigen matrix are NaN.
 *
 * @param x The Eigen matrix to check for NaN values.
 * @return True if all elements are NaN, otherwise false.
 */
template<typename Derived>
inline bool is_nan(const Eigen::MatrixBase<Derived>& x)
{
    return ((x.array() == x.array())).all();
}

/**
 * @brief Write an Eigen matrix to a CSV file.
 *
 * This function writes an Eigen matrix to a CSV file with the specified name.
 *
 * @param name The name of the CSV file to create.
 * @param matrix The Eigen matrix to be written to the file.
 */
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
void writeToCSVfile(std::string name, Eigen::MatrixXd matrix)
{
    std::ofstream file(name.c_str());
    file << matrix.format(CSVFormat);
}

/**
 * @brief Signal handler for capturing backtrace on segmentation fault.
 *
 * This function is a signal handler for capturing the backtrace when a segmentation fault occurs.
 *
 * @param sig The signal number.
 */
void handler(int sig) {
    void *array[10];
    size_t size;

    // Get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // Print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

/**
 * @brief Raise a segmentation fault (for debugging purposes).
 *
 * This function raises a segmentation fault, which can be used for debugging purposes
 * to capture the backtrace using the signal handler.
 */
void raise_segfault()
{
    int *foo = (int*)-1;        // Make a bad pointer
    printf("%d\n", *foo);       // Causes segfault
}

/**
 * @brief Export data to a CSV file.
 *
 * This function exports data to a CSV file with the specified file name.
 *
 * @param file_name The name of the CSV file to create.
 * @param data The data to be written to the file.
 */
template<typename t_file_name,
         typename t_data>
void export_data_to_csv_file
(
    const t_file_name& file_name,
    const t_data& data
)
{
    std::ofstream outfile;
    outfile.open(file_name);
    outfile << data << std::endl;
}
