/**
 * @file mapping.h
 * @brief Header file containing functions and data structures related to mapping and transformations in the UVLM framework.
 */

#pragma once

#include "triads.h"
#include "types.h"

/**
 * @namespace UVLM
 * @brief Namespace for the UVLM (Unsteady Vortex Lattice Method) framework.
 */
namespace UVLM
{    
    /**
     * @namespace Mapping
     * @brief Namespace for functions and data structures related to mapping and transformations.
     */
    namespace Mapping
    {
       /**
         * @brief Matrix containing the mapping from corner index to matrix indices.
         *
         * This matrix defines the mapping of corner indices to matrix indices as follows:
         * vortex_indices = [0, 0
         *                   1, 0,
         *                   1, 1,
         *                   0, 1]
         *
         * With the corner numbering as:
         *
         *       N -->
         *   0---------3
         *   |         |
         *   |         |
         *   1---------2
         *
         * So, the first element (0) has the coordinates (0, 0).
         */
        const Eigen::Matrix<unsigned int, 4, 2>
                vortex_indices((Eigen::Matrix<unsigned int, 4, 2>()
                                        << 0,0,1,0,1,1,0,1).finished());


        /**
         * @brief Perform bilinear mapping on input and output vectors.
         *
         * This function performs bilinear mapping on the input and output vectors.
         *
         * @tparam t_in Type of the input vector.
         * @tparam t_out Type of the output vector.
         * @param in Input vector.
         * @param out Output vector.
         */
        template <typename t_in, typename t_out>
        void BilinearMapping(t_in& in,
                             t_out& out)
        {
            const unsigned int ndims = in.size();
            for (unsigned int idim=0; idim<ndims; ++idim)
            {
                UVLM::Triads::BilinearMap1d(in[idim],
                                            out[idim]);
            }
        }
        /**
         * @brief Map double matrices to a vector of map objects.
         *
         * This function maps double matrices to a vector of map objects.
         *
         * @param dimensions Dimensions of the matrices.
         * @param in Input matrices.
         * @param map Vector of map objects.
         * @param correction Correction value for dimensions.
         * @param n_dim Number of dimensions (default is UVLM::Constants::NDIM).
         */
        void map_VecVecMat(const UVLM::Types::VecDimensions& dimensions,
                           double** in,
                           UVLM::Types::VecVecMapX& map,
                           const int& correction=0,
                           const unsigned int& n_dim=UVLM::Constants::NDIM)
        {
            const unsigned int n_surf = dimensions.size();
            map.resize(n_surf);
            unsigned int counter = 0;
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                for (uint i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    map[i_surf].push_back(UVLM::Types::MapMatrixX (in[counter],
                                                                   dimensions[i_surf].first + correction,
                                                                   dimensions[i_surf].second + correction));
                    counter++;
                }
            }
        }
        /**
         * @brief Map double matrices to a vector of map objects.
         *
         * This function maps double matrices to a vector of map objects.
         *
         * @param dimensions Dimensions of the matrices.
         * @param in Input matrices.
         * @param map Vector of map objects.
         * @param correction Correction value for dimensions.
         */
        void map_VecMat(const UVLM::Types::VecDimensions& dimensions,
                        double** in,
                        UVLM::Types::VecMapX& map,
                        const int& correction=0)
        {
            const unsigned int n_surf = dimensions.size();
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                map.push_back(UVLM::Types::MapMatrixX (in[i_surf],
                                                       dimensions[i_surf].first + correction,
                                                       dimensions[i_surf].second + correction));
            }
        }
        /**
         * @brief Map double matrices to a vector of vector map objects.
         *
         * This function maps double matrices to a vector of vector map objects.
         *
         * @param dimensions Dimensions of the matrices.
         * @param in Input matrices.
         * @param map Vector of vector map objects.
         * @param correction Correction value for dimensions.
         */
        void map_VecVec1(const UVLM::Types::VecDimensions& dimensions,
                        double** in,
                        UVLM::Types::VecMapVX& map,
                        const int& correction=0)
        {
            // Generates a variable that will be indexed as map[i_surf](i_m)
            const unsigned int n_surf = dimensions.size();
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                map.push_back(UVLM::Types::MapVectorX (in[i_surf],
                                                       dimensions[i_surf].first + correction));
            }
        }
        /**
         * @brief Map a double array to a VectorX object.
         *
         * This function maps a double array to a VectorX object.
         *
         * @param N_rows Number of rows in the array.
         * @param in Input array.
         * @param out Output VectorX object.
         * @param correction Correction value for dimensions.
        */
        void map_VecX(const uint N_rows,
                        double* in,
                        UVLM::Types::VectorX& out,
                        const int& correction=0)
        {
            // Caution: Use map VecX only for small vectors like p_rbm_vel_g
            out.resize(N_rows);
            for (uint i_row=0; i_row < N_rows; ++i_row)
            {
                out[i_row] = in[i_row];   
            }
        }
        /**
         * @brief Transform dimensions from a double array to a vector of dimensions.
         *
         * This function transforms dimensions from a double array to a vector of dimensions.
         *
         * @param n_surf Number of surfaces.
         * @param dimensions_in Input dimensions as a double array.
         * @param dimensions Output vector of dimensions.
         */
        void transform_dimensions(unsigned int& n_surf,
                                  unsigned int** dimensions_in,
                                  UVLM::Types::VecDimensions& dimensions)
        {
            dimensions.resize(n_surf);
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                dimensions[i_surf].first = dimensions_in[i_surf][0];
                dimensions[i_surf].second= dimensions_in[i_surf][1];
            }
        }

    }
}
