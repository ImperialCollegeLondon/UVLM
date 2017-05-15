#pragma once

#include "EigenInclude.h"
#include "types.h"

namespace UVLM
{
    namespace PostProc
    {
        template <typename t_zeta,
                  typename t_gamma,
                  typename t_uext,
                  typename t_forces>
        void calculate_static_forces
        (
            const t_zeta& zeta,
            const t_gamma& gamma,
            const t_uext& uext,
            t_forces&  forces,
            const UVLM::Types::FlightConditions& flightconditions
        )
        {
            // not bothered with effciency.
            // if it is so critical, it could be improved
            const uint n_surf = zeta.size();
            const UVLM::Types::Real q_inf =  0.5
                                            *flightconditions.rho
                                            *flightconditions.uinf
                                            *flightconditions.uinf
                                            *flightconditions.c_ref;

            UVLM::Types::Vector3 dl;
            UVLM::Types::Vector3 v;
            UVLM::Types::Vector3 f;
            uint start;
            uint end;
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const uint M = gamma[i_surf].rows();
                const uint N = gamma[i_surf].cols();

                for (uint i_M=0; i_M<M; ++i_M)
                {
                    for (uint i_N=0; i_N<N; ++i_N)
                    {
                        dl.setZero();
                        UVLM::Types::Vector3 v1;
                        UVLM::Types::Vector3 v2;
                        const unsigned int n_segment = 4;
                        for (unsigned int i_segment=0; i_segment<n_segment; ++i_segment)
                        {
                            unsigned int start = i_segment;
                            unsigned int end = (start + 1)%n_segment;
                            uint i_start = i_M + UVLM::Mapping::vortex_indices(start, 0);
                            uint j_start = i_N + UVLM::Mapping::vortex_indices(start, 1);
                            uint i_end = i_M + UVLM::Mapping::vortex_indices(end, 0);
                            uint j_end = i_N + UVLM::Mapping::vortex_indices(end, 1);

                            v1 << zeta[i_surf][0](i_start, j_start),
                                  zeta[i_surf][1](i_start, j_start),
                                  zeta[i_surf][2](i_start, j_start);
                            v2 << zeta[i_surf][0](i_end, j_end),
                                  zeta[i_surf][1](i_end, j_end),
                                  zeta[i_surf][2](i_end, j_end);
                            dl = v2-v1;

                            v << 0.5*(uext[i_surf][0](i_start, j_start) +
                                      uext[i_surf][0](i_end, j_end)),
                                 0.5*(uext[i_surf][1](i_start, j_start) +
                                      uext[i_surf][1](i_end, j_end)),
                                 0.5*(uext[i_surf][2](i_start, j_start) +
                                      uext[i_surf][2](i_end, j_end));

                            f = flightconditions.rho*gamma[i_surf](i_M, i_N)*v.cross(dl)/q_inf;

                            // transfer forces to matrix
                            // there are no moments
                            for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                            {
                                forces[i_surf][i_dim](i_start, j_start) +=
                                    0.5*f(i_dim);
                                forces[i_surf][i_dim](i_end, j_end) +=
                                    0.5*f(i_dim);
                            }
                        }
                    }
                }
            }
        }
    }
}
