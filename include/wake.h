#pragma once

#include "EigenInclude.h"
#include "types.h"

namespace UVLM
{
    namespace Wake
    {
        template <typename t_zeta,
                  typename t_zeta_star>
        void init_steady_wake
        (
            const t_zeta& zeta,
            t_zeta_star zeta_star
        )
        {
            // first implementation assumes the
            // wake is convected in the (1,0,0)
            // direction
            UVLM::Types::Vector3 dir_stream << 1, 0, 0
            
        }
    }
}
