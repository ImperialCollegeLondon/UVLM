#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "geometry.h"
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
                      typename t_gamma_star>
            void circulation_transfer
            (
                const t_zeta& zeta,
                const t_zeta_star& zeta_star,
                const t_gamma& gamma,
                t_gamma_star& gamma_star,
                const int in_n_rows = -1
            )
            {
                uint N, M, Mstar, n_rows;
                UVLM::Types::Real panel_area, te_panel_area;

                const uint n_surf = gamma.size();
                for (uint i_surf=0; i_surf<n_surf; ++i_surf)
                {
                    N = gamma_star[i_surf].cols();
                    Mstar = gamma_star[i_surf].rows();
                    M = gamma[i_surf].rows();
                    if (in_n_rows == -1){
                        n_rows = Mstar;
                    } else {
                        n_rows = in_n_rows;
                    }

                    for (uint i_n=0; i_n<N; ++i_n)
                    {
                        te_panel_area = UVLM::Geometry::panel_area
                            (
                                zeta[i_surf][0].template block<2,2>(M - 1, i_n),
                                zeta[i_surf][1].template block<2,2>(M - 1, i_n),
                                zeta[i_surf][2].template block<2,2>(M - 1, i_n)
                            );
                        for (uint i_m=0; i_m<n_rows; ++i_m)
                        {
                            panel_area = UVLM::Geometry::panel_area
                                (
                                    zeta[i_surf][0].template block<2,2>(i_m, i_n),
                                    zeta[i_surf][1].template block<2,2>(i_m, i_n),
                                    zeta[i_surf][2].template block<2,2>(i_m, i_n)
                                );
                            gamma_star[i_surf](i_m, i_n) = gamma[i_surf](M - 1, i_n)*panel_area/te_panel_area;
                        }
                    }
                }
            }

            template <typename t_zeta_star,
                      typename t_gamma_star,
                      typename t_extra_gamma_star,
                      typename t_extra_zeta_star,
                      typename t_dist_to_orig,
                      typename t_wake_conv_vel>
            void cfl_n1
            (
                const UVLM::Types::UVMopts& options,
                t_zeta_star& zeta_star,
                t_gamma_star& gamma_star,
                t_extra_gamma_star& extra_gamma_star,
                t_extra_zeta_star& extra_zeta_star,
                const t_dist_to_orig& dist_to_orig,
                t_wake_conv_vel& wake_conv_vel,
                double dt=0.
            )
            {
                unsigned int M, N;
                double cfl, dist;
                UVLM::Types::VectorX dist_to_orig_conv, coord0, coord1, coord2;
                UVLM::Types::VectorX new_coord0, new_coord1, new_coord2;
                UVLM::Types::Vector3 point;
                UVLM::Types::Real total_dist;

                // Allocate zeta_star_conv
                const uint n_surf = gamma_star.size();
                // UVLM::Types::VecVecMatrixX zeta_star_conv;
                // zeta_star_conv.resize(n_surf);
                // for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
                // {
                //     for (unsigned int i_dim=0; i_dim<3; ++i_dim)
                //     {
                //         zeta_star_conv[i_surf].push_back(UVLM::Types::MatrixX(gamma_star[i_surf].rows() + 1,
                //                                                                gamma_star[i_surf].cols() + 1));
                //     }
                // }

                // backup the wake shape
                // UVLM::Types::copy_VecVecMat(zeta_star, zeta_star_conv);

                // Loop through the surfaces
                for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
                {
                    // Retrieve array sizes
                    M = gamma_star[i_surf].rows();
                    N = gamma_star[i_surf].cols();
                    // Allocate vectors
                    dist_to_orig_conv.resize(M + 2);
                    dist_to_orig_conv.setZero();

                    // Loop through streamline vortices
                    for (unsigned int i_n=0; i_n<N + 1; ++i_n)
                    {
                        // Recompute the geometry
                        if ((options.convection_scheme == 2) or (options.convection_scheme == 3))
                        {
                            // Compute the distance of each wake vertice to the first point of the
                            dist_to_orig_conv(0) = 0.;
                            for (unsigned int i_m=1; i_m<M+1; ++i_m)
                            {
                                point << zeta_star[i_surf][0](i_m, i_n) - zeta_star[i_surf][0](i_m - 1, i_n),
                                         zeta_star[i_surf][1](i_m, i_n) - zeta_star[i_surf][1](i_m - 1, i_n),
                                         zeta_star[i_surf][2](i_m, i_n) - zeta_star[i_surf][2](i_m - 1, i_n);

                                dist_to_orig_conv(i_m) = point.norm() + dist_to_orig_conv(i_m - 1);
                            }
                            // Compute the last point
                            point << extra_zeta_star[i_surf][0](0, i_n) - zeta_star[i_surf][0](M, i_n),
                                     extra_zeta_star[i_surf][1](0, i_n) - zeta_star[i_surf][1](M, i_n),
                                     extra_zeta_star[i_surf][2](0, i_n) - zeta_star[i_surf][2](M, i_n);
                            dist_to_orig_conv(M + 1) = point.norm() + dist_to_orig_conv(M);
                            total_dist = dist_to_orig_conv(M+1) + 0.;
                            // Set maximum to one
                            for (unsigned int i_m=0; i_m<M+2; ++i_m)
                            {
                                dist_to_orig_conv(i_m) /= dist_to_orig_conv(M + 1);
                            }

                            // Change of coordinates
                            if (options.interp_coords == 0)
                            {
                                // Cartesian coordinates
                                for (unsigned int i_m=0; i_m<M+1; ++i_m)
                                {
                                    coord0(i_m) = zeta_star[i_surf][0](i_m, i_n);
                                    coord1(i_m) = zeta_star[i_surf][1](i_m, i_n);
                                    coord2(i_m) = zeta_star[i_surf][2](i_m, i_n);
                                }
                                coord0(M + 1) = extra_zeta_star[i_surf][0](0, i_n);
                                coord1(M + 1) = extra_zeta_star[i_surf][1](0, i_n);
                                coord2(M + 1) = extra_zeta_star[i_surf][2](0, i_n);

                            }else if (options.interp_coords == 1)
                            {
                                // Cylindrical coordinates along the z axis
                                for (unsigned int i_m=0; i_m<M+1; ++i_m)
                                {
                                    coord0(i_m) = std::sqrt(zeta_star[i_surf][0](i_m, i_n)*zeta_star[i_surf][0](i_m, i_n) +
                                                      zeta_star[i_surf][1](i_m, i_n)*zeta_star[i_surf][1](i_m, i_n));
                                    coord1(i_m) = atan2(zeta_star[i_surf][1](i_m, i_n), zeta_star[i_surf][0](i_m, i_n));
                                    coord2(i_m) = zeta_star[i_surf][2](i_m, i_n);
                                }
                                coord0(M + 1) = std::sqrt(extra_zeta_star[i_surf][0](0, i_n)*extra_zeta_star[i_surf][0](0, i_n) +
                                                  extra_zeta_star[i_surf][1](0, i_n)*extra_zeta_star[i_surf][1](0, i_n));
                                coord1(M + 1) = atan2(extra_zeta_star[i_surf][1](0, i_n), extra_zeta_star[i_surf][0](0, i_n));
                                coord2(M + 1) = extra_zeta_star[i_surf][2](0, i_n);
                            } else
                            {
                                std::cerr << "interp_coords == "
                                          << options.interp_coords
                                          << " is not supported by the UVLM solver. \n"
                                          << "Supported options are from [0->1]"
                                          << std::endl;
                            }

                            // Filter the values
                            if (options.filter_method == 0)
                            {
                                // No filter
                                // UVLM::Filter::
                                continue;
                            } else if (options.filter_method == 1)
                            {
                                // Splines
                                UVLM::Filters::splines(M + 1,
                                                       dist_to_orig_conv,
                                                       coord0, coord1, coord2);

                            } else
                            {
                                std::cerr << "filter_method == "
                                          << options.filter_method
                                          << " is not supported by the UVLM solver. \n"
                                          << "Supported options are from [0->1]"
                                          << std::endl;
                            }

                            // Redefine the location of the vertices
                            if (options.interp_method == 0)
                            {
                                // Linear interpolation
                                UVLM::Interpolation::linear(M + 1,
                                                            dist_to_orig[i_surf].col(i_n), dist_to_orig_conv,
                                                            coord0, coord1, coord2,
                                                            new_coord0, new_coord1, new_coord2);
                            } else if (options.interp_method == 1)
                            {
                                // Parabolic interpolation
                                UVLM::Interpolation::parabolic(M + 1,
                                                            dist_to_orig[i_surf].col(i_n), dist_to_orig_conv,
                                                            coord0, coord1, coord2,
                                                            new_coord0, new_coord1, new_coord2);
                            } else if (options.interp_method == 2)
                            {
                                // Splines interpolation
                                UVLM::Interpolation::splines(M + 1,
                                                            dist_to_orig[i_surf].col(i_n), dist_to_orig_conv,
                                                            coord0, coord1, coord2,
                                                            new_coord0, new_coord1, new_coord2);
                            } else
                            {
                                std::cerr << "interp_method == "
                                          << options.interp_method
                                          << " is not supported by the UVLM solver. \n"
                                          << "Supported options are from [0->1]"
                                          << std::endl;
                            }

                            // Change the coordinates back
                            if (options.interp_coords == 0)
                            {
                                // Cartesian coordinates
                                for (unsigned int i_m=0; i_m<M+2; ++i_m)
                                {
                                    zeta_star[i_surf][0](i_m, i_n) = new_coord0(i_m);
                                    zeta_star[i_surf][1](i_m, i_n) = new_coord1(i_m);
                                    zeta_star[i_surf][2](i_m, i_n) = new_coord2(i_m);
                                }
                            }else if (options.interp_coords == 1)
                            {
                                // Cylindrical coordinates
                                for (unsigned int i_m=0; i_m<M+2; ++i_m)
                                {
                                    zeta_star[i_surf][0](i_m, i_n) = new_coord0(i_m)*std::cos(new_coord1(i_m));
                                    zeta_star[i_surf][1](i_m, i_n) = new_coord0(i_m)*std::sin(new_coord1(i_m));
                                    zeta_star[i_surf][2](i_m, i_n) = new_coord2(i_m);
                                }
                            }

                        } // end if convection schemes

                        // Convect circulation and velocities
                        if (i_n < N)
                        {
                            // The first value (i_m=0) has already been adequately assigned in the convection step
                            for (unsigned int i_m=1; i_m<M-1; ++i_m)
                            {
                                dist = 0.25*(dist_to_orig[i_surf](i_m + 1, i_n) - dist_to_orig[i_surf](i_m - 1, i_n) +
                                             dist_to_orig[i_surf](i_m + 1, i_n + 1) - dist_to_orig[i_surf](i_m - 1, i_n + 1));
                                cfl = dt*wake_conv_vel[i_surf](i_m, i_n)/dist/total_dist;
                                gamma_star[i_surf](i_m, i_n) = (1. - cfl)*gamma_star[i_surf](i_m + 1, i_n) +
                                                                  cfl*gamma_star[i_surf](i_m, i_n);
                                wake_conv_vel[i_surf](i_m, i_n) = (1. - cfl)*wake_conv_vel[i_surf](i_m + 1, i_n) +
                                                                  cfl*wake_conv_vel[i_surf](i_m, i_n);
                            }
                            dist = 0.25*(dist_to_orig[i_surf](M, i_n) - dist_to_orig[i_surf](M - 1, i_n) +
                                         dist_to_orig[i_surf](M, i_n + 1) - dist_to_orig[i_surf](M - 1, i_n + 1));
                            cfl = dt*wake_conv_vel[i_surf](M - 1, i_n)/dist/total_dist;
                            gamma_star[i_surf](M - 1, i_n) = (1. - cfl)*extra_gamma_star[i_surf](0, i_n) +
                                                                  cfl*gamma_star[i_surf](M - 1, i_n);
                            // wake_conv_vel[i_surf](M - 1, i_n) = wake_conv_vel[i_surf](M - 2, i_n);
                        }
                      } // end of N loop
                  } // end of surfaces loop
            } // end cfl_n1 function
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
