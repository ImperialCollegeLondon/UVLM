#pragma once
#include "types.h"
#include "debugutils.h"

#include <iostream>

// DECLARATIONS
namespace UVLM
{
    namespace Symmetry
    {
        inline void flip_sign_component_VecVecMat
        (
            UVLM::Types::VecVecMatrixX& mat,
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
        inline void flip_sign_component_VecMat
        (
            UVLM::Types::VecMatrixX& mat
        )
        {
            uint M, N;
            uint n_surf = mat.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                M = mat[i_surf].rows();
                N = mat[i_surf].cols();
                for (uint i=0; i<M; ++i)
                {
                    for (uint j=0; j<N; ++j)
                    {
                        mat[i_surf](i,j) *= -1.;
                    }
                }
            }    
        }      
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
            const int& mirrored_component,
            UVLM::Types::VecVecMatrixX&  zeta_symmetry
        )
        {
            UVLM::Types::allocate_VecVecMat(zeta_symmetry, zeta);
            UVLM::Types::copy_VecVecMat(zeta, zeta_symmetry);
            flip_sign_component_VecVecMat(zeta_symmetry, mirrored_component);
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
            flip_sign_component_VecVecMat(zeta_symmetry, 2);
        }
        template <typename t_gamma>
        void generate_symmetric_gamma_grid
        (
            const t_gamma& gamma,
            UVLM::Types::VecMatrixX& gamma_symmetry
        )
        {
            // TODO: Check possibility of flipping zeta to have same gamma?
            UVLM::Types::allocate_VecMat(gamma_symmetry, gamma);
            UVLM::Types::copy_VecMat(gamma, gamma_symmetry);
            flip_sign_component_VecMat(gamma_symmetry);
        }
    }
}