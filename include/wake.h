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
            t_zeta_star zeta_star,
            const UVLM::Types::FlightConditions& flightconditions
        )
        {
            // wake convected in freestream direction
            UVLM::Types::Vector3 dir_stream(
                          flightconditions.uinf_direction);
            // UVLM::Types::Vector3 dir_stream;
            // dir_stream << 1, 0, 0;

            const uint n_surf = zeta.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const uint m_chordwise_panels = zeta[i_surf][0].rows() - 1;
                const uint n_spanwise_panels = zeta[i_surf][0].cols() - 1;
                for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                {
                    for (uint j=0; j<n_spanwise_panels; ++j)
                    {
                        auto b_block = zeta[i_surf][i_dim].template \
                                    block<2,2>(m_chordwise_panels-1, j);

                        // point 1 of the b is point 0 of the wake
                        // point 2 of the b is point 3 of the wake
                        zeta_star[i_surf][i_dim](0, j) = b_block(1, 0);
                        zeta_star[i_surf][i_dim](0, j+1) = b_block(1, 1);
                        zeta_star[i_surf][i_dim](1, j) = b_block(1, 0) + dir_stream(i_dim);
                        zeta_star[i_surf][i_dim](1, j+1) = b_block(1, 1) + dir_stream(i_dim);
                    }
                }
            }
        }


        template <typename t_gamma,
                  typename t_gamma_star>
        void steady_wake_circulation
        (
            const t_gamma& gamma,
            t_gamma_star& gamma_star
        )
        {
            const uint n_surf = gamma.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                gamma_star[i_surf] = gamma[i_surf].template bottomRows<1>();
            }
        }
    }
}
