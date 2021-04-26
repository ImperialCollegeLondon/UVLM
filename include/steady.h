#pragma once

#include "types.h"
#include "constants.h"
#include "triads.h"
#include "mapping.h"
#include "geometry.h"
#include "biotsavart.h"
#include "matrix.h"
#include "wake.h"
#include "postproc.h"
#include "linear_solver.h"

#include "debugutils.h"

#include <iostream>

// DECLARATIONS
namespace UVLM
{
    namespace Steady
    {
        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_uext,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_forces,
                  typename t_rbm_vel_g,
                  typename t_centre_rot_g>
        void solver
        (
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_forces& forces,
            t_rbm_vel_g& rbm_vel_g,
            t_centre_rot_g& centre_rot_g,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );

        template <typename t_zeta,
                  typename t_zeta_col,
                  typename t_uext_col,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_normals>
        void solve_horseshoe
        (
            t_zeta& zeta,
            t_zeta_col& zeta_col,
            t_uext_col& uext_col,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_normals& normals,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );
        template <typename t_zeta,
                  typename t_zeta_col,
                  typename t_uext_col,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_normals>
        void solve_discretised
        (
            t_zeta& zeta,
            t_zeta_col& zeta_col,
            t_uext_col& uext_col,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_normals& normals,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );
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
          typename t_forces,
          typename t_rbm_vel_g,
          typename t_centre_rot_g>
void UVLM::Steady::solver
(
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_forces& forces,
    t_rbm_vel_g& rbm_vel_g,
    t_centre_rot_g& centre_rot_g,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // Generate collocation points info
    //  Declaration
    UVLM::Types::VecVecMatrixX zeta_col;
    UVLM::Types::VecVecMatrixX uext_col;
    UVLM::Types::VecVecMatrixX uext_total;
    UVLM::Types::VecVecMatrixX uext_total_col;

    //  Allocation and mapping
    UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
    UVLM::Geometry::generate_colocationMesh(uext, uext_col);

    UVLM::Types::allocate_VecVecMat(uext_total, uext);
    UVLM::Types::copy_VecVecMat(uext, uext_total);
    UVLM::Types::allocate_VecVecMat(uext_total_col, uext, -1);

    // panel normals
    UVLM::Types::VecVecMatrixX normals;
    UVLM::Types::allocate_VecVecMat(normals, zeta_col);
    UVLM::Geometry::generate_surfaceNormal(zeta, normals);

    UVLM::Types::Vector3 u_steady;
    u_steady << uext[0][0](0,0),
                uext[0][1](0,0),
                uext[0][2](0,0);
    double delta_x = u_steady.norm()*options.dt;

    // total stream velocity
    UVLM::Unsteady::Utils::compute_resultant_grid_velocity
    (
        zeta,
        zeta_dot,
        uext,
        rbm_vel_g,
        centre_rot_g,
        uext_total
    );

    UVLM::Geometry::generate_colocationMesh(uext_total, uext_total_col);

    // if options.horseshoe, it is finished.
    if (options.horseshoe)
    {
        // solve the steady horseshoe problem
        UVLM::Steady::solve_horseshoe
        (
            zeta,
            zeta_col,
            uext_total_col,
            zeta_star,
            gamma,
            gamma_star,
            normals,
            options,
            flightconditions
        );

        UVLM::PostProc::calculate_static_forces
        (
            zeta,
            zeta_star,
            gamma,
            gamma_star,
            uext_total,
            forces,
            options,
            flightconditions
        );
        UVLM::Wake::Horseshoe::to_discretised(zeta_star,
                                              gamma_star,
                                              delta_x);
        return;
    }


    // create Wake
    // UVLM::Wake::Horseshoe::init(zeta, zeta_star, flightconditions);
    // UVLM::Wake::Horseshoe::to_discretised(zeta_star,
    //                                       gamma_star,
    //                                       delta_x);

    UVLM::Steady::solve_discretised
    (
        zeta,
        zeta_col,
        uext_total_col,
        zeta_star,
        gamma,
        gamma_star,
        normals,
        options,
        flightconditions
    );

    double zeta_star_norm_first = 0.0;
    double zeta_star_norm_previous = 0.0;
    double zeta_star_norm = 0.0;
    unsigned int N;

    (void)zeta_star_norm_previous; // used to silence warning

    UVLM::Types::VecVecMatrixX zeta_star_previous;
    if (options.n_rollup != 0)
    {
        zeta_star_norm_first = UVLM::Types::norm_VecVec_mat(zeta_star);
        zeta_star_norm_previous = zeta_star_norm_first;
        zeta_star_norm = 0.0;

        UVLM::Types::allocate_VecVecMat(zeta_star_previous, zeta_star);
        UVLM::Types::copy_VecVecMat(zeta_star, zeta_star_previous);
    }

    // ROLLUP LOOP--------------------------------------------------------
    for (uint i_rollup=0; i_rollup<options.n_rollup; ++i_rollup)
    {
        // determine convection velocity u_ind
        UVLM::Types::VecVecMatrixX u_ind;
        UVLM::Types::allocate_VecVecMat(u_ind,
                                        zeta_star);
        // induced velocity by vortex rings
        UVLM::BiotSavart::total_induced_velocity_on_wake(
            zeta,
            zeta_star,
            gamma,
            gamma_star,
            u_ind,
            options.ImageMethod,
            options.vortex_radius_wake_ind);
        // Do not move the vertices in the TE
        for (uint i_surf=0; i_surf<zeta_star.size(); ++i_surf)
        {
            N = zeta_star[i_surf][0].cols();
            for (uint i_n=0; i_n<N; ++i_n)
            {
                for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                {
                    u_ind[i_surf][i_dim](0, i_n) = 0.;
                }
            }
        }

        // convect based on u_ind for all the grid.
        UVLM::Wake::Discretised::convect(zeta_star,
                                         u_ind,
                                         options.dt);

        // generate AIC again
        if (i_rollup%options.rollup_aic_refresh == 0)
        {
            UVLM::Steady::solve_discretised
            (
                zeta,
                zeta_col,
                uext_total_col,
                zeta_star,
                gamma,
                gamma_star,
                normals,
                options,
                flightconditions
            );
        }

        // convergence check -------------------
        zeta_star_norm = UVLM::Types::norm_VecVec_mat(zeta_star);
        if (i_rollup != 0)
        {
            // double eps = std::abs((zeta_star_norm - zeta_star_norm_previous)
            //                       /zeta_star_norm_first);
            double eps = std::abs(UVLM::Types::norm_VecVec_mat(zeta_star - zeta_star_previous))/zeta_star_norm_first;
            std::cout << "    UVLM: Rollup iteration: " << i_rollup << ". Error: " << eps << std::endl;
            if (eps < options.rollup_tolerance)
            {
                break;
            }
            zeta_star_norm_previous = zeta_star_norm;
            UVLM::Types::copy_VecVecMat(zeta_star, zeta_star_previous);
        }
    }

    UVLM::PostProc::calculate_static_forces_unsteady
    (
        zeta,
        zeta_dot,
        zeta_star,
        gamma,
        gamma_star,
        uext,
        rbm_vel_g,
        centre_rot_g,
        forces,
        options,
        flightconditions
    );
}



/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta,
          typename t_zeta_col,
          typename t_uext_col,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_normals>
void UVLM::Steady::solve_horseshoe
(
    t_zeta& zeta,
    t_zeta_col& zeta_col,
    t_uext_col& uext_col,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_normals& normals,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // wake generation for horseshoe initialisation
    UVLM::Wake::Horseshoe::init(zeta,
                                zeta_star,
                                flightconditions);

    const uint n_surf = options.NumSurfaces;
    // size of rhs
    uint ii = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        uint M = uext_col[i_surf][0].rows();
        uint N = uext_col[i_surf][0].cols();

        ii += M*N;
    }

    const uint Ktotal = ii;
    // RHS generation
    UVLM::Types::VectorX rhs;
    UVLM::Matrix::RHS(zeta_col,
                      zeta_star,
                      uext_col,
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
                      uext_col,
                      normals,
                      options,
                      true,
                      aic);

    UVLM::Types::VectorX gamma_flat;
    UVLM::Matrix::deconstruct_gamma(gamma,
                                    gamma_flat,
                                    zeta_col);

    UVLM::LinearSolver::solve_system
    (
        aic,
        rhs,
        options,
        gamma_flat
    );

    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(gamma_flat,
                                    gamma,
                                    zeta_col);

    // copy gamma from trailing edge to wake if steady solution
    UVLM::Wake::Horseshoe::circulation_transfer(gamma,
                                                gamma_star);

}



/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta,
          typename t_zeta_col,
          typename t_uext_col,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_normals>
void UVLM::Steady::solve_discretised
(
    t_zeta& zeta,
    t_zeta_col& zeta_col,
    t_uext_col& uext_col,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_normals& normals,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    const uint n_surf = options.NumSurfaces;
    // size of rhs
    uint ii = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        uint M = uext_col[i_surf][0].rows();
        uint N = uext_col[i_surf][0].cols();

        ii += M*N;
    }
    const uint Ktotal = ii;

    UVLM::Types::VectorX rhs;
    UVLM::Types::MatrixX aic = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
    // RHS generation
    UVLM::Matrix::RHS(zeta_col,
                      zeta_star,
                      uext_col,
                      gamma_star,
                      normals,
                      options,
                      rhs,
                      Ktotal);

    // AIC generation
    UVLM::Matrix::AIC(Ktotal,
                      zeta,
                      zeta_col,
                      zeta_star,
                      uext_col,
                      normals,
                      options,
                      false,
                      aic);

    // linear system solution
    UVLM::Types::VectorX gamma_flat;
    UVLM::Matrix::deconstruct_gamma(gamma,
                                    gamma_flat,
                                    zeta_col);


    UVLM::LinearSolver::solve_system
    (
        aic,
        rhs,
        options,
        gamma_flat
    );

    // gamma flat to gamma
    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(gamma_flat,
                                    gamma,
                                    zeta_col);

    // copy gamma from trailing edge to wake
    int in_n_rows = -1;
    // if (!options.Steady) {in_n_rows = 1;}
    // UVLM::Wake::Horseshoe::circulation_transfer(gamma,
    //                                             gamma_star,
    //                                             in_n_rows);
    if (options.Steady) {
        UVLM::Wake::Horseshoe::circulation_transfer(gamma,
                                                gamma_star,
                                                in_n_rows);
    }
}
