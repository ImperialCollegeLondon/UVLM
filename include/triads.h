#pragma once
#include "types.h"

#include <iostream>

namespace UVLM
{
    namespace Triads
    {
        template <typename t_mat, typename t_out>
        inline void BilinearMap1d(const t_mat& mat,
                                  t_out& out)
        {
            const int n_rows = out.rows();
            const int n_cols = out.cols();

            for (int i=0; i<n_rows; ++i)
            {
                for (int j=0; j<n_cols; ++j)
                {
                    out(i,j) = mat.template block<2,2>(i,j).mean();
                }
            }
        }

        // transfer information from the collocation
        // points to the cornerpoints
        // mat has 1 less element per dimension than out
        template <typename t_mat,
                  typename t_out>
        inline void InvBilinearMap1d(const t_mat& mat,
                                     t_out& out)
        {
            const int n_rows = mat.rows();
            const int n_cols = mat.cols();

            for (int i=0; i<n_rows; ++i)
            {
                for (int j=0; j<n_cols; ++j)
                {
                    out.template block<2,2>(i,j) +=
                        0.25*mat(i,j);
                }
            }
        }

        template <typename t_1,
                  typename t_2,
                  typename t_out>
        void VecVecMatrix_difference(const t_1& mat1,
                                     const t_2& mat2,
                                     t_out& mat_out)
        {
            for (unsigned int i_surf=0; i_surf<mat1.size(); ++i_surf)
            {
                for (unsigned int i_dim=0; i_dim<mat1[i_surf].size(); ++i_dim)
                {
                    mat_out[i_surf][i_dim].noalias() = mat1[i_surf][i_dim] - mat2[i_surf][i_dim];
                }
            }
        }


        template <typename t_1,
                  typename t_2,
                  typename t_out>
        void VecVecMatrix_addition(const t_1& mat1,
                                    const t_2& mat2,
                                    t_out& mat_out)
        {
            for (unsigned int i_surf=0; i_surf<mat1.size(); ++i_surf)
            {
                for (unsigned int i_dim=0; i_dim<mat1[i_surf].size(); ++i_dim)
                {
                    mat_out[i_surf][i_dim].noalias() = mat1[i_surf][i_dim] + mat2[i_surf][i_dim];
                }
            }
        }
    }
}
