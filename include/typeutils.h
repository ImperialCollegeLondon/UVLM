#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "triads.h"

#include <cassert>

/**
 * @file typeutils.h
 * @brief Contains overloaded operators for vector-matrix operations.
 */

/**
 * @brief Subtraction operator for VecVecMatrixX types.
 *
 * Subtracts one VecVecMatrixX from another and returns the result.
 *
 * @param v1 The first VecVecMatrixX operand.
 * @param v2 The second VecVecMatrixX operand.
 * @return The result of the subtraction operation.
 */
UVLM::Types::VecVecMatrixX operator-
    (const UVLM::Types::VecVecMatrixX& v1,
     const UVLM::Types::VecVecMatrixX& v2)
{
    UVLM::Types::VecVecMatrixX vout;
    UVLM::Types::allocate_VecVecMat(vout, v1);

    UVLM::Triads::VecVecMatrix_difference(v1, v2, vout);
    return vout;
}

/**
 * @brief Subtraction operator for VecVecMapX and VecVecMatrixX types.
 *
 * Subtracts a VecVecMatrixX from a VecVecMapX and returns the result.
 *
 * @param v1 The VecVecMapX operand.
 * @param v2 The VecVecMatrixX operand.
 * @return The result of the subtraction operation.
 */
UVLM::Types::VecVecMatrixX operator-
    (const UVLM::Types::VecVecMapX& v1,
     const UVLM::Types::VecVecMatrixX& v2)
{
    UVLM::Types::VecVecMatrixX vout;
    UVLM::Types::allocate_VecVecMat(vout, v1);

    UVLM::Triads::VecVecMatrix_difference(v1, v2, vout);
    return vout;
}

/**
 * @brief Addition operator for VecVecMatrixX types.
 *
 * Adds two VecVecMatrixX and returns the result.
 *
 * @param v1 The first VecVecMatrixX operand.
 * @param v2 The second VecVecMatrixX operand.
 * @return The result of the addition operation.
 */
UVLM::Types::VecVecMatrixX operator+
    (const UVLM::Types::VecVecMatrixX& v1,
     const UVLM::Types::VecVecMatrixX& v2)
{
    UVLM::Types::VecVecMatrixX vout;
    UVLM::Types::allocate_VecVecMat(vout, v1);

    UVLM::Triads::VecVecMatrix_addition(v1, v2, vout);
    return vout;
}

/**
 * @brief Compound addition operator for VecVecMatrixX types.
 *
 * Adds a VecVecMatrixX to the current VecVecMatrixX and modifies the current object.
 *
 * @param own The VecVecMatrixX to be modified.
 * @param rhs The VecVecMatrixX to be added.
 * @return A reference to the modified VecVecMatrixX.
 */
UVLM::Types::VecVecMatrixX& operator+=
    (UVLM::Types::VecVecMatrixX& own,
     const UVLM::Types::VecVecMatrixX& rhs)
{
    unsigned int n_surf = own.size();
    for (unsigned int i_surf = 0; i_surf < n_surf; ++i_surf)
    {
        unsigned int n_dim = own[i_surf].size();
        for (unsigned int i_dim = 0; i_dim < n_dim; ++i_dim)
        {
            own[i_surf][i_dim] += rhs[i_surf][i_dim];
        }
    }
    return own;
}
