#pragma once
#include "types.h"
#include "debugutils.h"

#include <iostream>

// DECLARATIONS
namespace UVLM
{
    namespace Symmetry
    {
        template <typename t_mat>
        inline void flip_sign_component_VecVecMat
        (
            t_mat& mat,
            const uint& component
        )
        {
            uint M, N;
            uint n_surf = mat.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                M = mat[i_surf][0].rows();
                N = mat[i_surf][0].cols();
                for (uint i=0; i<M; ++i)
                {
                    for (uint j=0; j<N; ++j)
                    {
                        mat[i_surf][component](i,j) *= -1.;
                    }
                }
            }    
        }

        template <typename t_zeta>
        void generate_symmetric_surface_grids
        (
            const t_zeta& zeta,
            UVLM::Types::VecVecMatrixX&  zeta_symmetry
        )
        {
            UVLM::Types::allocate_VecVecMat(zeta_symmetry, zeta);
            UVLM::Types::copy_VecVecMat(zeta, zeta_symmetry);
            flip_sign_component_VecVecMat(zeta_symmetry, 1);
        }
    }
}