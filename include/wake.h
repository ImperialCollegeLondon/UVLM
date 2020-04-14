#pragma once

#include "EigenInclude.h"
#include "types.h"
#include <math.h>
// #include "unsteady.h"
// #include "steady.h"

namespace UVLM
{
    namespace Wake
    {
        namespace General
        {
            template <typename t_mat>
            void displace_VecVecMat
            (
                t_mat& mat
            )
            {
                const uint n_surf = mat.size();
                for (uint i_surf=0; i_surf<n_surf; ++i_surf)
                {
                    const uint n_rows = mat[i_surf][0].rows();
                    const uint n_dim = mat[i_surf].size();
                    for (uint i_dim=0; i_dim<n_dim; ++i_dim)
                    {
                        for (uint i_row=n_rows - 1; i_row>0; --i_row)
                        {
                            mat[i_surf][i_dim].row(i_row) =
                                mat[i_surf][i_dim].row(i_row - 1);
                        }
                        mat[i_surf][i_dim].template topRows<1>().setZero();
                    }
                }
            }

            template <typename t_mat>
            void displace_VecMat
            (
                t_mat& mat
            )
            {
                const uint n_surf = mat.size();
                for (uint i_surf=0; i_surf<n_surf; ++i_surf)
                {
                const uint n_rows = mat[i_surf].rows();
                    for (uint i_row=n_rows - 1; i_row>0; --i_row)
                    {
                    mat[i_surf].row(i_row) =
                        mat[i_surf].row(i_row - 1);
                    }
                mat[i_surf].template topRows<1>().setZero();
                }
            }
        }

        namespace Discretised
        {
            template <typename t_zeta_star,
                      typename t_zeta>
            void generate_new_row
            (
                t_zeta_star& zeta_star,
                const t_zeta& zeta
            )
            {
                const uint n_surf = zeta_star.size();
                for (uint i_surf=0; i_surf<n_surf; ++i_surf)
                {
                    for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                    {
                        zeta_star[i_surf][i_dim].template topRows<1>() =
                            zeta[i_surf][i_dim].template bottomRows<1>();
                    }
                }
            }

            /*******************************************************************
            Given a velocity field at the vertices of the grid, convect it
            following:
            x^{n+1} = x^{n} + u_ind \cdot dt
            *******************************************************************/
            template <typename t_zeta_star,
                      typename t_u_ind>
            void convect
            (
                t_zeta_star& zeta_star,
                const t_u_ind& u_ind,
                const double& delta_t
            )
            {
                const uint n_surf = zeta_star.size();
                for (uint i_surf=0; i_surf<n_surf; ++i_surf)
                {
                    for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                    {
                        const uint n_M = zeta_star[i_surf][i_dim].rows();
                        const uint n_N = zeta_star[i_surf][i_dim].cols();
                        for (uint i_M=0; i_M<n_M; ++i_M)
                        {
                            for (uint j_N=0; j_N<n_N; ++j_N)
                            {
                                zeta_star[i_surf][i_dim](i_M, j_N) += (
                                    u_ind[i_surf][i_dim](i_M, j_N)*delta_t
                                );
                            }
                        }
                    }
                }
            }


            template <typename t_zeta,
                      typename t_zeta_star,
                      typename t_gamma,
                      typename t_gamma_star,
                      typename t_uext_total_col>
            void circulation_transfer
            (
                const t_zeta& zeta,
                const t_zeta_star& zeta_star,
                const t_gamma& gamma,
                t_gamma_star& gamma_star,
                const t_uext_total_col& uext_total_col,
                double dt=0.
            )
            {
                uint n_cols, M;
                double cfl;
                UVLM::Types::Vector3 vel, dist;

                const uint n_surf = gamma.size();
                for (uint i_surf=0; i_surf<n_surf; ++i_surf)
                {
                    n_cols = gamma[i_surf].cols();
                    M = gamma[i_surf].rows();
                    for (uint i_n=0; i_n<n_cols; ++i_n)
                    {
                        dist << 0.25*(zeta_star[i_surf][0](1, i_n) + zeta_star[i_surf][0](1, i_n+1)
                                        - zeta[i_surf][0](M-1, i_n) - zeta[i_surf][0](M-1, i_n+1)),
                                0.25*(zeta_star[i_surf][1](1, i_n) + zeta_star[i_surf][1](1, i_n+1)
                                        - zeta[i_surf][1](M-1, i_n) - zeta[i_surf][1](M-1, i_n+1)),
                                0.25*(zeta_star[i_surf][2](1, i_n) + zeta_star[i_surf][2](1, i_n+1)
                                        - zeta[i_surf][2](M-1, i_n) - zeta[i_surf][2](M-1, i_n+1));
                        // delta = 0.25*pow((zeta_star[i_surf][0](1, i_n) + zeta_star[i_surf][0](1, i_n+1)
                        //                 - zeta[i_surf][0](-2, i_n) - zeta[i_surf][0](-2, i_n+1))*
                        //              (zeta_star[i_surf][1](1, i_n) + zeta_star[i_surf][1](1, i_n+1)
                        //                 - zeta[i_surf][1](-2, i_n) - zeta[i_surf][1](-2, i_n+1))*
                        //              (zeta_star[i_surf][2](1, i_n) + zeta_star[i_surf][2](1, i_n+1)
                        //                 - zeta[i_surf][2](-2, i_n) - zeta[i_surf][2](-2, i_n+1)), 1./3.);
                        vel << uext_total_col[i_surf][0](M-1, i_n),
                               uext_total_col[i_surf][1](M-1, i_n),
                               uext_total_col[i_surf][2](M-1, i_n);
                        cfl = dt*vel.norm()/dist.norm();
                        // if(cfl > 1.){
                        //     std::cout << "WARNING: CFL=" << cfl << " > 1 at isurf:" << i_surf << " in:" << i_n << std::endl;
                        // }
                        // std::cout << "dist" << dist << std::endl;
                        // std::cout << "vel" << vel << std::endl;
                        // std::cout << "cfl" << cfl << std::endl;
                        gamma_star[i_surf](0, i_n) = (1. - cfl)*gamma_star[i_surf](1, i_n) +
                                                       cfl*gamma[i_surf](M-1, i_n);
                        // The wake has already been convected so I should use gamma_star[i_surf](1, i_n)
                        // but in the new convection scheme I am copying back the old value to the first cell because otherwise its lost
                    }
                }
            }

        }
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
            void to_discretised
            (
                t_zeta_star& zeta_star,
                t_gamma_star& gamma_star,
                const double& delta_x
            )
            {
                UVLM::Types::Vector3 dir_stream;
                dir_stream << zeta_star[0][0](1, 0) - zeta_star[0][0](0, 0),
                              zeta_star[0][1](1, 0) - zeta_star[0][1](0, 0),
                              zeta_star[0][2](1, 0) - zeta_star[0][2](0, 0);
                dir_stream.normalize();
                UVLM::Types::Vector3 delta_x_vec = dir_stream*delta_x;

                const uint n_surf = zeta_star.size();
                for (uint i_surf=0; i_surf<n_surf; ++i_surf)
                {
                    const uint n_spanwise_panels =
                        zeta_star[i_surf][0].cols() - 1;
                    const uint mstar =
                        zeta_star[i_surf][0].rows() - 1;
                    // Now the wake shape generation will be dealt with from SHARPy
                    // It might be acceptable to keep it for horseshoe cases
                    // for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                    // {
                    //     for (uint j=0; j<n_spanwise_panels + 1; ++j)
                    //     {
                    //         for (uint i=1; i<mstar + 1; ++i)
                    //         {
                    //             zeta_star[i_surf][i_dim](i, j) =
                    //                 zeta_star[i_surf][i_dim](i - 1, j)\
                    //                     + delta_x_vec(i_dim);
                    //         }
                    //     }
                    // }
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

            // It can be used for unsteady too with in_n_rows = 1
            // (only one row to be copied from the trailing edge)
            template <typename t_gamma,
                      typename t_gamma_star>
            void circulation_transfer
            (
                const t_gamma& gamma,
                t_gamma_star& gamma_star,
                const int in_n_rows = -1
            )
            {

                const uint n_surf = gamma.size();
                for (uint i_surf=0; i_surf<n_surf; ++i_surf)
                {
                    uint n_rows = in_n_rows;
                    if (n_rows == -1)
                    {
                        n_rows = gamma_star[i_surf].rows();
                    }
                    for (uint i_m=0; i_m<n_rows; ++i_m)
                    {
                        gamma_star[i_surf].row(i_m) = gamma[i_surf].template bottomRows<1>();
                    }
                }
            }
        }// End namespace horseshoe
    }
}
