#pragma once

#include "EigenInclude.h"
#include "types.h"

namespace UVLM
{
    namespace PostProc
    {
        template <typename t_zeta,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_uext,
                  typename t_forces>
        void calculate_static_forces
        (
            const t_zeta& zeta,
            const t_zeta_star& zeta_star,
            const t_gamma& gamma,
            const t_gamma_star& gamma_star,
            const t_uext& uext,
            t_forces&  forces,
            const UVLM::Types::VMopts options,
            const UVLM::Types::FlightConditions& flightconditions
        )
        {
            // first calculate all the velocities at the corner points
            UVLM::Types::VecVecMatrixX velocities;
            UVLM::Types::allocate_VecVecMat(velocities, zeta);
            // free stream contribution
            UVLM::Types::copy_VecVecMat(uext, velocities);

            // not bothered with effciency.
            // if it is so critical, it could be improved
            const uint n_surf = zeta.size();
            UVLM::Types::Vector3 dl;
            UVLM::Types::Vector3 v;
            UVLM::Types::Vector3 f;
            UVLM::Types::Vector3 v_ind;
            UVLM::Types::Vector3 rp;
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
                        UVLM::Types::Vector3 r1;
                        UVLM::Types::Vector3 r2;
                        const unsigned int n_segment = 4;
                        for (unsigned int i_segment=0; i_segment<n_segment; ++i_segment)
                        {
                            if ((i_segment == 1) && (i_M == M - 1))
                            {
                                // trailing edge
                                continue;
                            }
                            unsigned int start = i_segment;
                            unsigned int end = (start + 1)%n_segment;
                            uint i_start = i_M + UVLM::Mapping::vortex_indices(start, 0);
                            uint j_start = i_N + UVLM::Mapping::vortex_indices(start, 1);
                            uint i_end = i_M + UVLM::Mapping::vortex_indices(end, 0);
                            uint j_end = i_N + UVLM::Mapping::vortex_indices(end, 1);


                            r1 << zeta[i_surf][0](i_start, j_start),
                                  zeta[i_surf][1](i_start, j_start),
                                  zeta[i_surf][2](i_start, j_start);
                            r2 << zeta[i_surf][0](i_end, j_end),
                                  zeta[i_surf][1](i_end, j_end),
                                  zeta[i_surf][2](i_end, j_end);

                            // position of the center point of the vortex filament
                            rp = 0.5*(r1 + r2);

                            // induced vel by vortices at vp
                            v_ind.setZero();
                            for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
                            {
                                UVLM::Types::VecMatrixX temp_uout;
                                UVLM::Types::allocate_VecMat(temp_uout,
                                                             zeta[ii_surf],
                                                             -1);
                                UVLM::BiotSavart::surface_with_horseshoe
                                (
                                    zeta[ii_surf],
                                    zeta_star[ii_surf],
                                    gamma[ii_surf],
                                    gamma_star[ii_surf],
                                    rp,
                                    temp_uout,
                                    options.ImageMethod
                                );
                                v_ind(0) += temp_uout[0].sum();
                                v_ind(1) += temp_uout[1].sum();
                                v_ind(2) += temp_uout[2].sum();
                            }

                            dl = r2-r1;

                            v << 0.5*(uext[i_surf][0](i_start, j_start) +
                                      uext[i_surf][0](i_end, j_end)),
                                 0.5*(uext[i_surf][1](i_start, j_start) +
                                      uext[i_surf][1](i_end, j_end)),
                                 0.5*(uext[i_surf][2](i_start, j_start) +
                                      uext[i_surf][2](i_end, j_end));

                            v = (v + v_ind).eval();

                            f = flightconditions.rho*gamma[i_surf](i_M, i_N)*v.cross(dl);

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
