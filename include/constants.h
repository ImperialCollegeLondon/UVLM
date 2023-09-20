/**
 * @file constants.h
 * @brief Header file containing definitions of constants used in the UVLM program.
 */

#pragma once

#include "types.h"

namespace UVLM
{
    namespace Constants
    {
        /**
         * @brief Constants used in the UVLM program.
         */
        
        /// The number of dimensions in the program (always set to 3 for 3D space).
        const unsigned int NDIM = 3;

        /// Mathematical constant pi (π) with high precision.
        const UVLM::Types::Real PI = 3.1415926535897932384626433832795028841971;

        /// Four times the value of pi (4π).
        const UVLM::Types::Real PI4 = 4.0 * PI;

        /// The reciprocal of four times pi (1 / 4π).
        const UVLM::Types::Real INV_PI4 = 1.0 / PI4;

        /// Conversion factor for degrees to radians (π / 180.0).
        const UVLM::Types::Real DEGREES2RAD = PI / 180.0;

        /// A small positive value used as a numerical tolerance (10 times machine epsilon).
        const UVLM::Types::Real EPSILON = 10 * std::numeric_limits<UVLM::Types::Real>::epsilon();
    }
}
