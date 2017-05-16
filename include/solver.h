#pragma once

#include "types.h"
#include "triads.h"
#include "constants.h"
#include "mapping.h"
#include "geometry.h"
#include "biotsavart.h"
#include "matrix.h"
#include "wake.h"
#include "postproc.h"

#include <iostream>

// DECLARATIONS
namespace UVLM
{
    namespace Solver
    {
        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_uext,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_normals,
                  typename t_forces>
        void solve
        (
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_normals& normals,
            t_forces& forces,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );

        template <typename t_in,
                  typename t_out>
        void generate_colocationMesh
        (
            t_in& vortex_mesh,
            t_out& collocation_mesh
        );
    }
}

// SOURCE CODE
/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_in,
          typename t_out>
void UVLM::Solver::generate_colocationMesh
(
    t_in& vortex_mesh,
    t_out& collocation_mesh
)
{
    // Size of surfaces contained in a vector of tuples
    UVLM::Types::VecDimensions dimensions;
    dimensions.resize(vortex_mesh.size());
    for (unsigned int i_surf=0; i_surf<dimensions.size(); ++i_surf)
    {
        dimensions[i_surf] = UVLM::Types::IntPair(
                                            vortex_mesh[i_surf][0].rows(),
                                            vortex_mesh[i_surf][0].cols());
    }

    if (collocation_mesh.empty())
    {
        UVLM::Types::allocate_VecVecMat(collocation_mesh,
                                        UVLM::Constants::NDIM,
                                        dimensions,
                                        -1);
    }
    for (unsigned int i_surf=0; i_surf<dimensions.size(); ++i_surf)
    {
        UVLM::Mapping::BilinearMapping(vortex_mesh[i_surf],
                                       collocation_mesh[i_surf]);
    }
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta,
          typename t_zeta_dot,
          typename t_uext,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_normals,
          typename t_forces>
void UVLM::Solver::solve
(
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_normals& normals,
    t_forces& forces,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // Generate collocation points info
    //  Declaration
    UVLM::Types::VecVecMatrixX zeta_col;
    UVLM::Types::VecVecMatrixX zeta_dot_col;
    UVLM::Types::VecVecMatrixX uext_col;
    UVLM::Types::VecVecMatrixX zeta_star_col;

    //  Allocation and mapping
    UVLM::Solver::generate_colocationMesh(zeta, zeta_col);
    UVLM::Solver::generate_colocationMesh(zeta_dot, zeta_dot_col);
    UVLM::Solver::generate_colocationMesh(uext, uext_col);
    UVLM::Solver::generate_colocationMesh(zeta_star, zeta_star_col);

    // panel normals
    UVLM::Geometry::generate_surfaceNormal(zeta, normals);

    // wake generation for steady solutions
    if (options.Steady)
    {
        UVLM::Wake::init_steady_wake(zeta, zeta_star, flightconditions);
    }


    // RHS generation
    UVLM::Types::VectorX rhs;
    unsigned int Ktotal;
    UVLM::Matrix::RHS(zeta_col,
                      zeta_star,
                      uext_col,
                      zeta_dot_col,
                      gamma_star,
                      normals,
                      options,
                      rhs,
                      Ktotal);

    // AIC generation
    UVLM::Types::MatrixX aic = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
    UVLM::Matrix::AIC(Ktotal,
                      zeta,
                      zeta_col,
                      zeta_star,
                      zeta_star_col,
                      uext_col,
                      normals,
                      options,
                      aic);

    // std::cout << rhs << std::endl;

    UVLM::Types::VectorX gamma_flat;
    gamma_flat = aic.partialPivLu().solve(rhs);

    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(gamma_flat,
                                    gamma,
                                    zeta_col,
                                    zeta_star_col,
                                    options);

    // copy gamma from trailing edge to wake if steady solution
    UVLM::Wake::steady_wake_circulation(gamma,
                                        gamma_star);

    // static forces
    UVLM::PostProc::calculate_static_forces
    (
        zeta,
        zeta_star,
        gamma,
        gamma_star,
        uext,
        forces,
        options,
        flightconditions
    );


    // // // temporary output of velocities at nodes
    // uint n_surf = gamma.size();
    // UVLM::Types::Vector3 dl;
    // UVLM::Types::Vector3 v;
    // UVLM::Types::Vector3 f;
    // UVLM::Types::Vector3 v_ind;
    // UVLM::Types::Vector3 vp;
    // uint start;
    // uint end;
    // for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    // {
    //     const uint M = gamma[i_surf].rows();
    //     const uint N = gamma[i_surf].cols();
    //
    //     for (uint i_M=0; i_M<M; ++i_M)
    //     {
    //         for (uint i_N=0; i_N<N; ++i_N)
    //         {
    //             UVLM::Types::Vector3 v1;
    //             const unsigned int n_segment = 4;
    //             for (unsigned int i_segment=0; i_segment<n_segment; ++i_segment)
    //             {
    //                 unsigned int start = i_segment;
    //                 unsigned int end = (start + 1)%n_segment;
    //                 uint i_start = i_M + UVLM::Mapping::vortex_indices(start, 0);
    //                 uint j_start = i_N + UVLM::Mapping::vortex_indices(start, 1);
    //                 uint i_end = i_M + UVLM::Mapping::vortex_indices(end, 0);
    //                 uint j_end = i_N + UVLM::Mapping::vortex_indices(end, 1);
    //
    //                 v1 << zeta[i_surf][0](i_start, j_start),
    //                       zeta[i_surf][1](i_start, j_start),
    //                       zeta[i_surf][2](i_start, j_start);
    //
    //                 // induced vel by vortices at v1
    //                 v_ind.setZero();
    //                 forces[i_surf][0](i_start, j_start) = 0.0;
    //                 forces[i_surf][1](i_start, j_start) = 0.0;
    //                 forces[i_surf][2](i_start, j_start) = 0.0;
    //                 for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
    //                 {
    //                     UVLM::Types::VecMatrixX temp_uout;
    //                     UVLM::Types::allocate_VecMat(temp_uout,
    //                                                  zeta[ii_surf],
    //                                                  -1);
    //                     UVLM::BiotSavart::surface_with_horseshoe
    //                     (
    //                         zeta[ii_surf],
    //                         zeta_star[ii_surf],
    //                         gamma[ii_surf],
    //                         gamma_star[ii_surf],
    //                         v1,
    //                         temp_uout,
    //                         options.ImageMethod
    //                     );
    //                     forces[i_surf][0](i_start, j_start) += temp_uout[0].sum();//+uext[i_surf][0](i_start, j_start);
    //                     forces[i_surf][1](i_start, j_start) += temp_uout[1].sum();//+uext[i_surf][1](i_start, j_start);
    //                     forces[i_surf][2](i_start, j_start) += temp_uout[2].sum();//+uext[i_surf][2](i_start, j_start);
    //                 }
    //             }
    //         }
    //     }
    // }
}
