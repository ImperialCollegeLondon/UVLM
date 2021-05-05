#pragma once

#include "types.h"
#include "triads.h"
#include "unsteady_utils.h"
#include "constants.h"
#include "mapping.h"
#include "geometry.h"
#include "biotsavart.h"
#include "matrix.h"
#include "postproc.h"
#include "steady.h"
#include "wake.h"

#include <iostream>

// DECLARATIONS
namespace UVLM
{
    namespace Unsteady
    {
        template <typename t_lifting_surfaces>
        void solver
        (
            const uint& i_iter,
            t_lifting_surfaces& lifting_surfaces_unsteady,
            const UVLM::Types::UVMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );

        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_zeta_star,
                  typename t_zeta_star_dot,
                  typename t_uext,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_rbm_velocity,
                  typename t_forces>
        void initialise
        (
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_zeta_star& zeta_star,
            t_zeta_star_dot& zeta_star_dot,
            t_uext& uext,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_rbm_velocity& rbm_velocity,
            t_forces& forces,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );
    }
}

// Initialisation routine for unsteady UVLM routines.
// It runs VLM in different configurations and provides
// a steady force calculation if required.
// Input conditions should be the same as for the first time step.
// In order to avoid force steps at initialisation for coupled simulations, the
// wake modelling should be as close as possible to the model
// used in the unsteady part.
template <typename t_zeta,
          typename t_zeta_dot,
          typename t_zeta_star,
          typename t_zeta_star_dot,
          typename t_uext,
          typename t_gamma,
          typename t_gamma_star,
          typename t_rbm_velocity,
          typename t_forces>
void UVLM::Unsteady::initialise
(
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_zeta_star& zeta_star,
    t_zeta_star_dot& zeta_star_dot,
    t_uext& uext,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_rbm_velocity& rbm_velocity,
    t_forces& forces,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // incident velocity taking into account RBM, zeta_dot:
    // UVLM::Types::VecVecMatrixX uext_resultant;
    // UVLM::Types::allocate_VecVecMat(uext_resultant,
    //                                 uext);
    // UVLM::Unsteady::Utils::compute_resultant_grid_velocity
    // (
    //     zeta,
    //     zeta_dot,
    //     uext,
    //     rbm_velocity,
    //     uext_resultant
    // );

    // call steady solver
    UVLM::Steady::solver
    (
        zeta,
        zeta_dot,
        uext,
        zeta_star,
        gamma,
        gamma_star,
        forces,
        rbm_velocity,
        options,
        flightconditions
    );
    return;
}


template <typename t_lifting_surfaces>
void UVLM::Unsteady::solver
(
    const uint& i_iter,
    t_lifting_surfaces& lifting_surfaces_unsteady,
    const UVLM::Types::UVMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // SOLVE------------------------------------------
    const uint n_surf = options.NumSurfaces;
    const double dt = options.dt;
    // Generate collocation points info
    //  Declaration
    lifting_surfaces_unsteady.get_surface_parameters();

    UVLM::Unsteady::Utils::compute_resultant_grid_velocity_solid_vel
    (
        lifting_surfaces_unsteady.zeta,
        lifting_surfaces_unsteady.zeta_dot,
        lifting_surfaces_unsteady.u_ext,
        lifting_surfaces_unsteady.rbm_vel_g,
        lifting_surfaces_unsteady.uext_total,
        lifting_surfaces_unsteady.solid_vel
    );

    //  Allocation and mapping
    // Same in steady    
    UVLM::Geometry::generate_colocationMesh(lifting_surfaces_unsteady.uext_total,
                                            lifting_surfaces_unsteady.uext_total_col);

    // Unsteady specific
    UVLM::Types::VecVecMatrixX uext_star_total;
    UVLM::Types::allocate_VecVecMat(uext_star_total, lifting_surfaces_unsteady.uext_star);

    // for what is extra gamma star?
    UVLM::Types::VecMatrixX extra_gamma_star;
    UVLM::Types::VecVecMatrixX extra_zeta_star;
    extra_zeta_star.resize(n_surf);
    for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
    {
        extra_gamma_star.push_back(UVLM::Types::MatrixX());
        extra_gamma_star[i_surf].setZero(1, lifting_surfaces_unsteady.gamma_star[i_surf].cols());
        for (unsigned int i_dim=0; i_dim<3; ++i_dim)
        {
            extra_zeta_star[i_surf].push_back(UVLM::Types::MatrixX(1,
                                                                   lifting_surfaces_unsteady.gamma_star[i_surf].cols() + 1));
        }
    }

    UVLM::Types::VMopts steady_options = UVLM::Types::UVMopts2VMopts(options);

    // different than in steady
    if (options.convect_wake)
    {
        UVLM::Unsteady::Utils::convect_unsteady_wake
        (
            options,
            lifting_surfaces_unsteady.zeta,
            lifting_surfaces_unsteady.zeta_star,
            lifting_surfaces_unsteady.gamma,
            lifting_surfaces_unsteady.gamma_star,
            lifting_surfaces_unsteady.u_ext,
            lifting_surfaces_unsteady.uext_star,
            uext_star_total,
            lifting_surfaces_unsteady.rbm_vel_g,
            extra_gamma_star,
            extra_zeta_star
        );
    }

    UVLM::Wake::Discretised::circulation_transfer(lifting_surfaces_unsteady.zeta,
                                                  lifting_surfaces_unsteady.zeta_star,
                                                  lifting_surfaces_unsteady.gamma,
                                                  lifting_surfaces_unsteady.gamma_star,
                                                  lifting_surfaces_unsteady.uext_total_col,
                                                  dt);

    if (!options.cfl1)
    {
        UVLM::Wake::Discretised::cfl_n1(options,
                                        lifting_surfaces_unsteady.zeta_star,
                                        lifting_surfaces_unsteady.gamma_star,
                                        extra_gamma_star,
                                        extra_zeta_star,
                                        lifting_surfaces_unsteady.dist_to_orig,
                                        uext_star_total,
                                        lifting_surfaces_unsteady.solid_vel,
                                        dt);
    }

    // we can use UVLM::Steady::solve_discretised if uext_col
    // is the total velocity including non-steady contributions.
    // TODO: Check if struct for unsteady surfaces can be transformed to steady to save mem and time?
    UVLM::Steady::solve_discretised
    (
        lifting_surfaces_unsteady,
        steady_options,
        flightconditions
    );

    if (options.quasi_steady)
    {
        UVLM::Wake::Horseshoe::circulation_transfer(lifting_surfaces_unsteady.gamma,
                                                    lifting_surfaces_unsteady.gamma_star,
                                                    -1);
    }

    // forces calculation
    // set forces to 0 just in case
    UVLM::Types::initialise_VecVecMat(lifting_surfaces_unsteady.forces);
    // static:
    UVLM::PostProc::calculate_static_forces_unsteady
    (
        lifting_surfaces_unsteady.zeta,
        lifting_surfaces_unsteady.zeta_dot,
        lifting_surfaces_unsteady.zeta_star,
        lifting_surfaces_unsteady.gamma,
        lifting_surfaces_unsteady.gamma_star,
        lifting_surfaces_unsteady.u_ext,
        lifting_surfaces_unsteady.rbm_vel_g,
        lifting_surfaces_unsteady.forces,
        steady_options,
        flightconditions
    );
    // TODO: Check if delete?
    UVLM::Types::initialise_VecVecMat(lifting_surfaces_unsteady.dynamic_forces);

    }

    // total stream velocity
    UVLM::Unsteady::Utils::compute_resultant_grid_velocity_solid_vel
    (
        zeta,
        zeta_dot,
        uext,
        rbm_velocity,
        uext_total,
        solid_vel
    );

    UVLM::Types::VMopts steady_options = UVLM::Types::UVMopts2VMopts(options);

    //  Allocation and mapping
    UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
    UVLM::Geometry::generate_colocationMesh(uext_total, uext_total_col);

    UVLM::Types::VecVecMatrixX uext_star_total;
    UVLM::Types::allocate_VecVecMat(uext_star_total, uext_star);

    // panel normals
    UVLM::Geometry::generate_surfaceNormal(zeta, normals);

    if (options.convect_wake)
    {
        UVLM::Unsteady::Utils::convect_unsteady_wake
        (
            options,
            zeta,
            zeta_star,
            gamma,
            gamma_star,
            uext,
            uext_star,
            uext_star_total,
            rbm_velocity,
            extra_gamma_star,
            extra_zeta_star
        );
    }

    UVLM::Wake::Discretised::circulation_transfer(zeta,
                                                  zeta_star,
                                                  gamma,
                                                  gamma_star,
                                                  uext_total_col,
                                                  dt);

    if (!options.cfl1)
    {
        UVLM::Wake::Discretised::cfl_n1(options,
                                        zeta_star,
                                        gamma_star,
                                        extra_gamma_star,
                                        extra_zeta_star,
                                        dist_to_orig,
                                        uext_star_total,
                                        solid_vel,
                                        dt);
    }

    // we can use UVLM::Steady::solve_discretised if uext_col
    // is the total velocity including non-steady contributions.
    UVLM::Steady::solve_discretised
    (
        zeta,
        zeta_col,
        uext_total_col,
        zeta_star,
        gamma,
        gamma_star,
        normals,
        steady_options,
        flightconditions
    );

    if (options.quasi_steady)
    {
        UVLM::Wake::Horseshoe::circulation_transfer(gamma,
                                                    gamma_star,
                                                    -1);
    }

    // forces calculation
    // set forces to 0 just in case
    UVLM::Types::initialise_VecVecMat(forces);
    // static:
    UVLM::PostProc::calculate_static_forces_unsteady
    (
        zeta,
        zeta_dot,
        zeta_star,
        gamma,
        gamma_star,
        uext,
        rbm_velocity,
        forces,
        steady_options,
        flightconditions
    );
    // dynamic::
    // if (i_iter > 0)
    // {
        UVLM::Types::initialise_VecVecMat(dynamic_forces);
        // std::cout << "Max dynamic forces:" << std::endl;
        // std::cout << dynamic_forces[0][2].maxCoeff() << std::endl;
        // calculate dynamic forces
        //std::cerr << "No unsteady forces will be calculated with the old routine!" << std::endl;
        // UVLM::PostProc::calculate_dynamic_forces
        // (
        //     zeta,
        //     zeta_star,
        //     zeta_col,
        //     gamma,
        //     gamma_star,
        //     previous_gamma,
        //     uext_total,
        //     normals,
        //     dynamic_forces,
        //     options,
        //     flightconditions
        // );
        // std::cout << dynamic_forces[0][2].maxCoeff() << std::endl;
    // }
}
