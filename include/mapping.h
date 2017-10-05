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


    }
}
