#pragma once

#include "types.h"

namespace UVLM
{
    namespace Constants
    {
        // Constants in the program
        const unsigned int NDIM = 3;
        const UVLM::Types::Real PI = 3.1415926535897932384626433832795028841971;
        const UVLM::Types::Real PI4 = 4.0*PI;
        const UVLM::Types::Real INV_PI4 = 1.0/(4.0*PI);
        const UVLM::Types::Real DEGREES2RAD = PI/180.0;

        const UVLM::Types::Real EPSILON =
                10*std::numeric_limits<UVLM::Types::Real>::epsilon();
    }
}
