#pragma once

#include "types.h"

#include <iostream>

/**
 * @file triads.h
 * @brief Contains utility functions for matrix operations and transformations.
 */

namespace UVLM
{
    namespace Triads
    {
        /**
         * @brief Perform a bilinear mapping on a matrix.
         *
         * This function calculates the mean of each 2x2 block of a matrix and stores the result in another matrix.
         *
         * @tparam t_mat The input matrix type.
         * @tparam t_out The output matrix type.
         * @param mat The input matrix to be mapped.
         * @param out The output matrix to store the mapped values.
         */
        template <typename t_mat, typename t_out>
        inline void BilinearMap1d(const t_mat& mat, t_out& out)
        {
            const int n_rows = out.rows();
            const int n_cols = out.cols();

            for (int i = 0; i < n_rows; ++i)
            {
                for (int j = 0; j < n_cols; ++j)
                {
                    out(i, j) = mat.template block<2, 2>(i, j).mean();
                }
            }
        }

        /**
         * @brief Perform an inverse bilinear mapping on a matrix.
         *
         * This function transfers information from the collocation points to the corner points in the matrix.
         * The input matrix has one less element per dimension than the output matrix.
         *
         * @tparam t_mat The input matrix type.
         * @tparam t_out The output matrix type.
         * @param mat The input matrix to be mapped.
         * @param out The output matrix to store the mapped values.
         */
        template <typename t_mat, typename t_out>
        inline void InvBilinearMap1d(const t_mat& mat, t_out& out)
        {
            const int n_rows = mat.rows();
            const int n_cols = mat.cols();

            for (int i = 0; i < n_rows; ++i)
            {
                for (int j = 0; j < n_cols; ++j)
                {
                    out.template block<2, 2>(i, j) += 0.25 * mat(i, j);
                }
            }
        }

        /**
         * @brief Calculate the difference between two matrices element-wise.
         *
         * @tparam t_1 The first matrix type.
         * @tparam t_2 The second matrix type.
         * @tparam t_out The output matrix type.
         * @param mat1 The first input matrix.
         * @param mat2 The second input matrix.
         * @param mat_out The output matrix to store the element-wise difference.
         */
        template <typename t_1, typename t_2, typename t_out>
        void VecVecMatrix_difference(const t_1& mat1, const t_2& mat2, t_out& mat_out)
        {
            for (unsigned int i_surf = 0; i_surf < mat1.size(); ++i_surf)
            {
                for (unsigned int i_dim = 0; i_dim < mat1[i_surf].size(); ++i_dim)
                {
                    mat_out[i_surf][i_dim].noalias() = mat1[i_surf][i_dim] - mat2[i_surf][i_dim];
                }
            }
        }

        /**
         * @brief Calculate the element-wise addition of two matrices.
         *
         * @tparam t_1 The first matrix type.
         * @tparam t_2 The second matrix type.
         * @tparam t_out The output matrix type.
         * @param mat1 The first input matrix.
         * @param mat2 The second input matrix.
         * @param mat_out The output matrix to store the element-wise sum.
         */
        template <typename t_1, typename t_2, typename t_out>
        void VecVecMatrix_addition(const t_1& mat1, const t_2& mat2, t_out& mat_out)
        {
            for (unsigned int i_surf = 0; i_surf < mat1.size(); ++i_surf)
            {
                for (unsigned int i_dim = 0; i_dim < mat1[i_surf].size(); ++i_dim)
                {
                    mat_out[i_surf][i_dim].noalias() = mat1[i_surf][i_dim] + mat2[i_surf][i_dim];
                }
            }
        }

        /**
         * @brief Add two matrices element-wise and update the first matrix.
         *
         * @tparam t_1 The first matrix type.
         * @tparam t_2 The second matrix type.
         * @param mat_in_and_out The input matrix to be updated.
         * @param mat_in The second input matrix to add.
         */
        template <typename t_1, typename t_2>
        void VecVecMatrix_addition(t_1& mat_in_and_out, const t_2& mat_in)
        {
            for (uint i_surf = 0; i_surf < mat_in_and_out.size(); ++i_surf)
            {
                for (uint i_dim = 0; i_dim < mat_in_and_out[i_surf].size(); ++i_dim)
                {
                    mat_in_and_out[i_surf][i_dim] += mat_in[i_surf][i_dim];
                }
            }
        }

        /**
         * @brief Calculate the element-wise difference between two matrices and update the first matrix.
         *
         * @tparam t_1 The first matrix type.
         * @tparam t_2 The second matrix type.
         * @param mat_in_and_out The input matrix to be updated.
         * @param mat_in The second input matrix for subtraction.
         */
        template <typename t_1, typename t_2>
        void VecVecMatrix_difference(t_1& mat_in_and_out, const t_2& mat_in)
        {
            for (uint i_surf = 0; i_surf < mat_in_and_out.size(); ++i_surf)
            {
                for (uint i_dim = 0; i_dim < mat_in_and_out[i_surf].size(); ++i_dim)
                {
                    mat_in_and_out[i_surf][i_dim] -= mat_in[i_surf][i_dim];
                }
            }
        }
    }
}
