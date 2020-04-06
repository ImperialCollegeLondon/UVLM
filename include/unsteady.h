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
        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_uext,
                  typename t_uext_star,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_normals,
                //   typename t_previous_gamma,
                  typename t_rbm_velocity,
                  typename t_forces>
        void solver
        (
            const uint& i_iter,
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_uext_star& uext_star,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_normals& normals,
            // const t_previous_gamma& previous_gamma,
            t_rbm_velocity& rbm_velocity,
            t_forces& forces,
            t_forces& dynamic_forces,
            const UVLM::Types::UVMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );

        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_uext,
                  typename t_uext_star,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_normals,
                //   typename t_previous_gamma,
                  typename t_rbm_velocity,
                  typename t_forces,
                  typename t_uindw>
        void solver_uindw
        (
            const uint& i_iter,
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_uext_star& uext_star,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_normals& normals,
            // const t_previous_gamma& previous_gamma,
            t_rbm_velocity& rbm_velocity,
            t_forces& forces,
            t_forces& dynamic_forces,
            t_uindw& uindw,
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

        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_uext,
                  typename t_uext_star,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_normals,
                //   typename t_previous_gamma,
                  typename t_rbm_velocity,
                  typename t_forces>
        void shw_solver
        (
            const uint& i_iter,
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_uext_star& uext_star,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_normals& normals,
            // const t_previous_gamma& previous_gamma,
            t_rbm_velocity& rbm_velocity,
            t_forces& forces,
            t_forces& dynamic_forces,
            const UVLM::Types::UVMopts& options,
            const UVLM::Types::FlightConditions& flightconditions,
            const UVLM::Types::SHWOptions& shwoptions
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
    UVLM::Types::VecVecMatrixX uext_resultant;
    UVLM::Types::allocate_VecVecMat(uext_resultant,
                                    uext);
    UVLM::Unsteady::Utils::compute_resultant_grid_velocity
    (
        zeta,
        zeta_dot,
        uext,
        rbm_velocity,
        uext_resultant
    );

    // call steady solver
    UVLM::Steady::solver
    (
        zeta,
        zeta_dot,
        uext_resultant,
        zeta_star,
        gamma,
        gamma_star,
        forces,
        options,
        flightconditions
    );
    return;
}


template <typename t_zeta,
          typename t_zeta_dot,
          typename t_uext,
          typename t_uext_star,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_normals,
        //   typename t_previous_gamma,
          typename t_rbm_velocity,
          typename t_forces>
void UVLM::Unsteady::solver
(
    const uint& i_iter,
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_uext_star& uext_star,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_normals& normals,
    // const t_previous_gamma& previous_gamma,
    t_rbm_velocity& rbm_velocity,
    t_forces& forces,
    t_forces& dynamic_forces,
    const UVLM::Types::UVMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // SOLVE------------------------------------------
    const uint n_surf = options.NumSurfaces;
    const double dt = options.dt;
    // Generate collocation points info
    //  Declaration
    UVLM::Types::VecVecMatrixX zeta_col;
    UVLM::Types::VecVecMatrixX uext_total;
    UVLM::Types::allocate_VecVecMat(uext_total, uext);
    UVLM::Types::copy_VecVecMat(uext, uext_total);

    UVLM::Types::VecVecMatrixX uext_total_col;
    UVLM::Types::allocate_VecVecMat(uext_total_col, uext, -1);

    // total stream velocity
    UVLM::Unsteady::Utils::compute_resultant_grid_velocity
    (
        zeta,
        zeta_dot,
        uext,
        rbm_velocity,
        uext_total
    );

    UVLM::Types::VMopts steady_options = UVLM::Types::UVMopts2VMopts(options);

    //  Allocation and mapping
    UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
    UVLM::Geometry::generate_colocationMesh(uext_total, uext_total_col);

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
            rbm_velocity
        );
    }

    UVLM::Wake::Discretised::circulation_transfer(zeta,
                                                  zeta_star,
                                                  gamma,
                                                  gamma_star,
                                                  uext_total_col,
                                                  dt);

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

template <typename t_zeta,
          typename t_zeta_dot,
          typename t_uext,
          typename t_uext_star,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_normals,
          typename t_rbm_velocity,
          typename t_forces,
          typename t_uindw>
void UVLM::Unsteady::solver_uindw
(
    const uint& i_iter,
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_uext_star& uext_star,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_normals& normals,
    // const t_previous_gamma& previous_gamma,
    t_rbm_velocity& rbm_velocity,
    t_forces& forces,
    t_forces& dynamic_forces,
    t_uindw& uindw,
    const UVLM::Types::UVMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // SOLVE------------------------------------------
    const uint n_surf = options.NumSurfaces;
    const double dt = options.dt;
    // Generate collocation points info
    //  Declaration
    UVLM::Types::VecVecMatrixX zeta_col;
    UVLM::Types::VecVecMatrixX uext_total;
    UVLM::Types::allocate_VecVecMat(uext_total, uext);
    UVLM::Types::copy_VecVecMat(uext, uext_total);

    UVLM::Types::VecVecMatrixX uext_total_col;
    UVLM::Types::allocate_VecVecMat(uext_total_col, uext, -1);

    // total stream velocity
    UVLM::Unsteady::Utils::compute_resultant_grid_velocity
    (
        zeta,
        zeta_dot,
        uext,
        rbm_velocity,
        uext_total
    );

    UVLM::Types::VMopts steady_options = UVLM::Types::UVMopts2VMopts(options);

    //  Allocation and mapping
    UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
    UVLM::Geometry::generate_colocationMesh(uext_total, uext_total_col);

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
            rbm_velocity
        );
    }

    UVLM::Wake::Discretised::circulation_transfer(zeta,
                                                  zeta_star,
                                                  gamma,
                                                  gamma_star,
                                                  uext_total_col,
                                                  dt);

    // we can use UVLM::Steady::solve_discretised if uext_col
    // is the total velocity including non-steady contributions.
    UVLM::Steady::solve_discretised_uindw
    (
        zeta,
        zeta_col,
        uext_total_col,
        zeta_star,
        gamma,
        gamma_star,
        normals,
        uindw,
        steady_options,
        flightconditions
    );

    UVLM::Types::initialise_VecVecMat(forces);
    UVLM::Types::initialise_VecVecMat(dynamic_forces);

    // // forces calculation
    // // set forces to 0 just in case
    // UVLM::Types::initialise_VecVecMat(forces);
    // // static:
    // UVLM::PostProc::calculate_static_forces_unsteady
    // (
    //     zeta,
    //     zeta_dot,
    //     zeta_star,
    //     gamma,
    //     gamma_star,
    //     uext,
    //     rbm_velocity,
    //     forces,
    //     steady_options,
    //     flightconditions
    // );
    // // dynamic::
    // // if (i_iter > 0)
    // // {
    //     UVLM::Types::initialise_VecVecMat(dynamic_forces);
    //     // std::cout << "Max dynamic forces:" << std::endl;
    //     // std::cout << dynamic_forces[0][2].maxCoeff() << std::endl;
    //     // calculate dynamic forces
    //     //std::cerr << "No unsteady forces will be calculated with the old routine!" << std::endl;
    //     // UVLM::PostProc::calculate_dynamic_forces
    //     // (
    //     //     zeta,
    //     //     zeta_star,
    //     //     zeta_col,
    //     //     gamma,
    //     //     gamma_star,
    //     //     previous_gamma,
    //     //     uext_total,
    //     //     normals,
    //     //     dynamic_forces,
    //     //     options,
    //     //     flightconditions
    //     // );
    //     // std::cout << dynamic_forces[0][2].maxCoeff() << std::endl;
    // // }
}
template <typename t_zeta,
          typename t_zeta_dot,
          typename t_uext,
          typename t_uext_star,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_normals,
        //   typename t_previous_gamma,
          typename t_rbm_velocity,
          typename t_forces>
void UVLM::Unsteady::shw_solver
(
    const uint& i_iter,
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_uext_star& uext_star,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_normals& normals,
    // const t_previous_gamma& previous_gamma,
    t_rbm_velocity& rbm_velocity,
    t_forces& forces,
    t_forces& dynamic_forces,
    const UVLM::Types::UVMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    const UVLM::Types::SHWOptions& shwoptions
)
{
    // SOLVE------------------------------------------
    const uint n_surf = options.NumSurfaces;
    // Generate collocation points info
    //  Declaration
    UVLM::Types::VecVecMatrixX zeta_col;
    UVLM::Types::VecVecMatrixX uext_total;
    UVLM::Types::allocate_VecVecMat(uext_total, uext);
    UVLM::Types::copy_VecVecMat(uext, uext_total);

    UVLM::Types::VecVecMatrixX uext_total_col;
    UVLM::Types::allocate_VecVecMat(uext_total_col, uext, -1);

    // total stream velocity
    UVLM::Unsteady::Utils::compute_resultant_grid_velocity
    (
        zeta,
        zeta_dot,
        uext,
        rbm_velocity,
        uext_total
    );

    UVLM::Types::VMopts steady_options = UVLM::Types::UVMopts2VMopts(options);
    steady_options.Steady = true;

    //  Allocation and mapping
    UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
    UVLM::Geometry::generate_colocationMesh(uext_total, uext_total_col);

    // panel normals
    UVLM::Geometry::generate_surfaceNormal(zeta, normals);

    // In principle, it is not going to be needed
    // // std::cout << options.convect_wake << std::endl;
    // if (options.convect_wake)
    // {
    //     // std::cout << "Convecting wake" << std::endl;
    //     UVLM::Unsteady::Utils::convect_unsteady_wake
    //     (
    //         options,
    //         zeta,
    //         zeta_star,
    //         gamma,
    //         gamma_star,
    //         uext,
    //         uext_star,
    //         rbm_velocity
    //     );
    // }

    // I am going to replace the call to UVLM::Steady::solve_discretised by
    // the whole code. This should be integrated at some point
    // BEGINNING

    // Define the two first rows of the wake (the second one will be overwritten afterwards)
    // UVLM::Wake::Horseshoe::init(zeta, zeta_star, flightconditions);
    // UVLM::Wake::SHW::to_discretised(zeta_star,
    //                                 gamma_star,
    //                                 flightconditions,
    //                                 shwoptions);

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

    // END

    // forces calculation
    // set forces to 0 just in case
    UVLM::Types::initialise_VecVecMat(forces);
    // UVLM::Types::
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
