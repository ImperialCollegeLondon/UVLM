#pragma once

#include "triads.h"
#include "types.h"

namespace UVLM
{
    namespace Mapping
    {
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
