/**
 * @file cpp_interface.h
 * @brief Header file containing definitions of interface functions (SHARPy <--> UVLM)
 */

#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "constants.h"
#include "geometry.h"
#include "steady.h"
#include "unsteady.h"
#include "struct_utils.h"
#include <iostream>

// Define DLLEXPORT as an extern "C" for C/C++ compatibility
#define DLLEXPORT extern "C"

// Uncomment the following lines if OpenMP is used
// #ifdef _OPENMP
// #include "omp.h"
// #endif

namespace UVLMlin {

    /**
     * @brief Call the derivative of the Biot-Savart panel method.
     *
     * This function calculates the derivative of the Biot-Savart panel method for a panel.
     * It computes the derivative of the induced velocity at a point due to a panel's circulation.
     *
     * @param p_DerP [out] - The calculated derivative of the induced velocity components (3x3 matrix).
     * @param p_DerVertices [in] - The coordinates of the panel's vertices (12 values - 4 vertices x, y, and z coordinates).
     * @param p_zetaP [in] - The coordinates of the point where the velocity is being calculated (3 values - x, y, and z coordinates).
     * @param p_ZetaPanel [in] - The coordinates of the panel's vertices (12 values - 4 vertices x, y, and z coordinates).
     * @param gamma [in] - The circulation strength of the panel.
     * @param vortex_radius [in/out] - The radius of the vortex core (modified during computation).
     */
    DLLEXPORT void call_der_biot_panel(double p_DerP[9],
        double p_DerVertices[36],
        double p_zetaP[3],
        double p_ZetaPanel[12],
        const double& gamma,
        double& vortex_radius);

    /**
     * @brief Call the Biot-Savart panel method.
     *
     * This function calculates the Biot-Savart panel method for a panel.
     * It computes the induced velocity at a point due to a panel's circulation.
     *
     * @param p_vel [out] - The calculated induced velocity components (3 values - x, y, and z components).
     * @param p_zetaP [in] - The coordinates of the point where the velocity is being calculated (3 values - x, y, and z coordinates).
     * @param p_ZetaPanel [in] - The coordinates of the panel's vertices (12 values - 4 vertices x, y, and z coordinates).
     * @param gamma [in] - The circulation strength of the panel.
     * @param vortex_radius [in/out] - The radius of the vortex core (modified during computation).
     */
    DLLEXPORT void call_biot_panel(double p_vel[3],
        double p_zetaP[3],
        double p_ZetaPanel[12],
        const double& gamma,
        double& vortex_radius);
}
