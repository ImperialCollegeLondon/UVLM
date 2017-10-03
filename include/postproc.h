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
            // Set forces to 0
            UVLM::Types::initialise_VecVecMat(forces);

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
                                UVLM::BiotSavart::surface_with_steady_wake
                                (
                                    zeta[ii_surf],
                                    zeta_star[ii_surf],
                                    gamma[ii_surf],
                                    gamma_star[ii_surf],
                                    rp,
                                    options.horseshoe,
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

        // Forces is not set to 0, forces are added
        // NOTE: uext has to include all the velocities contribution
        template <typename t_zeta,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_previous_gamma,
                  typename t_uext,
                  typename t_forces>
        void calculate_dynamic_forces
        (
            const t_zeta& zeta,
            const t_zeta_star& zeta_star,
            const t_gamma& gamma,
            const t_gamma_star& gamma_star,
            const t_previous_gamma& previous_gamma,
            const t_uext& uext,
            t_forces&  forces,
            const UVLM::Types::VMopts options,
            const UVLM::Types::FlightConditions& flightconditions
        )
        {
            // TODO not finished

            const UVLM::Types::Real dt = options.dt;
            const uint n_surf = zeta.size();

            UVLM::Types::VecVecMatrixX unsteady_force;
            UVLM::Types::allocate_VecVecMat(unsteady_force, forces, -1);

            // calculate gamma_dot
            UVLM::Types::VecMatrixX gamma_dot;
            UVLM::Types::allocate_VecMat(gamma_dot, gamma);
            // simple finite differences:
            // f' = (f+ - f-)/dt
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const uint n_rows = gamma[i_surf].rows();
                const uint n_cols = gamma[i_surf].cols();

                for (uint i=0; i<n_rows; ++i)
                {
                    for (uint j=0; j<n_cols; ++j)
                    {
                        gamma_dot[i_surf](i,j) =
                            (gamma[i_surf](i,j)
                             -
                             previous_gamma[i_surf](i,j))/dt;
                    }
                }
            }

            // transfer forces to vortex corners

        }

    }
}
