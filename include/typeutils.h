#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "triads.h"

#include <cassert>

// substraction operator
UVLM::Types::VecVecMatrixX operator-
        (const UVLM::Types::VecVecMatrixX& v1,
         const UVLM::Types::VecVecMatrixX& v2)
{
    UVLM::Types::VecVecMatrixX vout;
    UVLM::Types::allocate_VecVecMat(vout, v1);

    UVLM::Triads::VecVecMatrix_difference(v1, v2, vout);
    return vout;
}
UVLM::Types::VecVecMatrixX operator-
        (const UVLM::Types::VecVecMapX& v1,
         const UVLM::Types::VecVecMatrixX& v2)
{
    UVLM::Types::VecVecMatrixX vout;
    UVLM::Types::allocate_VecVecMat(vout, v1);

    UVLM::Triads::VecVecMatrix_difference(v1, v2, vout);
    return vout;
}


// addition operator
UVLM::Types::VecVecMatrixX operator+
        (const UVLM::Types::VecVecMatrixX& v1,
         const UVLM::Types::VecVecMatrixX& v2)
{
    UVLM::Types::VecVecMatrixX vout;
    UVLM::Types::allocate_VecVecMat(vout, v1);

    UVLM::Triads::VecVecMatrix_addition(v1, v2, vout);
    return vout;
}

UVLM::Types::VecVecMatrixX& operator+=
    (UVLM::Types::VecVecMatrixX& own,
     const UVLM::Types::VecVecMatrixX& rhs)
{
    unsigned int n_surf = own.size();
    for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
    {
        unsigned int n_dim = own[i_surf].size();
        for (unsigned int i_dim=0; i_dim<n_dim; ++i_dim)
        {
            own[i_surf][i_dim] += rhs[i_surf][i_dim];
        }
    }
    return own;
}
