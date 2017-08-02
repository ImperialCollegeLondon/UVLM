#pragma once

#include "EigenInclude.h"
#include "types.h"

namespace UVLM
{
    namespace Wake
    {
        namespace Horseshoe
        {
            template <typename t_zeta,
                      typename t_zeta_star>
            void init
            (
                const t_zeta& zeta,
                t_zeta_star zeta_star,
                const UVLM::Types::FlightConditions& flightconditions
            )
            {
                // wake convected in freestream direction
                UVLM::Types::Vector3 dir_stream(
                              flightconditions.uinf_direction);

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

            template <typename t_zeta_star,
                      typename t_gamma_star>
            void to_wake
            // TODO TO FINISH!!!! note: maybe it is finished
            (
                t_zeta_star& zeta_star,
                t_gamma_star& gamma_star,
                const double& delta_x
            )
            {
                UVLM::Types::Vector3 dir_stream;
                dir_stream << zeta_star[0][0](1, 0),
                              zeta_star[0][1](1, 0),
                              zeta_star[0][2](1, 0);
                dir_stream.normalize();
                UVLM::Types::Vector3 delta_x_vec = dir_stream*delta_x;

                const uint n_surf = zeta_star.size();
                for (uint i_surf=0; i_surf<n_surf; ++i_surf)
                {
                    const uint n_spanwise_panels =
                        zeta_star[i_surf][0].cols() - 1;
                    const uint mstar =
                        zeta_star[i_surf][0].rows() - 1;

                    for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                    {
                        for (uint j=0; j<n_spanwise_panels + 1; ++j)
                        {
                            for (uint i=1; i<mstar + 1; ++i)
                            {
                                zeta_star[i_surf][i_dim](i, j) =
                                    zeta_star[i_surf][i_dim](i - 1, j)\
                                        + delta_x_vec(i_dim);
                            }
                        }
                    }
                    std::cout << "   " << std::endl;
                    std::cout << gamma_star[i_surf].rows() << std::endl;
                    std::cout << gamma_star[i_surf].cols() << std::endl;
                    for (uint j=0; j<n_spanwise_panels; ++j)
                    {
                        for (uint i=1; i<mstar; ++i)
                        {
                            gamma_star[i_surf](i, j) = \
                                gamma_star[i_surf](i-1, j);
                        }
                    }
                }
            }

            template <typename t_gamma,
                      typename t_gamma_star>
            void circulation_transfer
            (
                const t_gamma& gamma,
                t_gamma_star& gamma_star
            )
            {
                const uint n_surf = gamma.size();
                for (uint i_surf=0; i_surf<n_surf; ++i_surf)
                {
                    gamma_star[i_surf].template topRows<1>() = gamma[i_surf].template bottomRows<1>();
                }
            }
        }
    }
}
