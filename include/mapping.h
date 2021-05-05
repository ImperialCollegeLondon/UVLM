#pragma once

#include "triads.h"
#include "types.h"

namespace UVLM
{
    namespace Mapping
    {
        // this matrix contains the mapping from the corner index to the
        // matrix indices.
        // It contains:
        // vortex_indices = [0, 0
        //                   1, 0,
        //                   1, 1,
        //                   0, 1]
        // With the numbering as:
        //          N -->
        //      0---------3
        //   M  |         |
        //   |  |         |
        //   V  1---------2
        // so, the first element (0), has the coordinate (for example)
        // indices of: vortex_indices(0), that is, 0 and 0
        const Eigen::Matrix<unsigned int, 4, 2>
                vortex_indices((Eigen::Matrix<unsigned int, 4, 2>()
                                        << 0,0,1,0,1,1,0,1).finished());



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

        void map_VecVec1(const UVLM::Types::VecDimensions& dimensions,
                        double** in,
                        UVLM::Types::VecMapVX& map,,
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
        void map_VecX(const uint N_rows,
                        double* in,
                        UVLM::Types::VectorX& out,
                        //UVLM::Types::MapVectorX& map,
                        const int& correction=0)
        {
            // Caution: Use map VecX only for small vectors like p_rbm_vel_g
            out.resize(N_rows);
            for (uint i_row=0; i_row < N_rows; ++i_row)
            {
                out[i_row] = in[i_row];   
            }
        }

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
