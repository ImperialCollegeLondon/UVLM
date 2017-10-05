#pragma once

#include "types.h"
#include "triads.h"
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
                  typename t_previous_gamma,
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
            const t_previous_gamma& previous_gamma,
            t_rbm_velocity& rbm_velocity,
            t_forces& forces,
            t_forces& dynamic_forces,
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
                  typename t_rbm_velocity,
                  typename t_uext_out>
        void compute_resultant_grid_velocity
        (
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_rbm_velocity& rbm_velocity,
            t_uext_out& uext_out
        );
        template <typename t_zeta,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_uext,
                  typename t_uext_star,
                  typename t_rbm_velocity>
        void convect_unsteady_wake
        (
            const UVLM::Types::UVMopts& options,
            const t_zeta& zeta,
            t_zeta_star& zeta_star,
            const t_gamma& gamma,
            t_gamma_star& gamma_star,
            const t_uext& uext,
            const t_uext_star& uext_star,
            const t_rbm_velocity& rbm_velocity
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
    UVLM::Unsteady::compute_resultant_grid_velocity
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
          typename t_previous_gamma,
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
    const t_previous_gamma& previous_gamma,
    t_rbm_velocity& rbm_velocity,
    t_forces& forces,
    t_forces& dynamic_forces,
    const UVLM::Types::UVMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
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
    UVLM::Unsteady::compute_resultant_grid_velocity
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
    UVLM::Types::VecVecMatrixX normals;
    UVLM::Types::allocate_VecVecMat(normals, zeta_col);
    UVLM::Geometry::generate_surfaceNormal(zeta, normals);

    convect_unsteady_wake
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

    // // wake convection
    // // the UVMopts flag convection_scheme determines how the
    // // wake is convected.
    // // convection_scheme == 0 => prescribed and fixed
    // // convection_scheme == 1 => prescribed following deformations of the wing
    // // convection_scheme == 2 => free, convection based on u_ext
    // // convection_scheme == 3 => free, convection based on u_ext and induced velocities.
    // if (options.convection_scheme == 0)
    // {
    //     UVLM::Wake::General::displace_VecMat(gamma_star);
    // } else if (options.convection_scheme == 1)
    // {
    //     std::cerr << "convection_scheme == "
    //               << options.convection_scheme
    //               << " is not yet implemented in the UVLM solver"
    //               << std::endl;
    // } else if (options.convection_scheme == 2)
    // {
    //     UVLM::Types::VecVecMatrixX uext_star_total;
    //     UVLM::Types::allocate_VecVecMat(uext_star_total, uext_star);
    //     UVLM::Types::VecVecMatrixX zeros;
    //     UVLM::Types::allocate_VecVecMat(zeros, uext_star);
    //     // total stream velocity
    //     UVLM::Unsteady::compute_resultant_grid_velocity
    //     (
    //         zeta_star,
    //         zeros,
    //         uext_star,
    //         rbm_velocity,
    //         uext_star_total
    //     );
    //     // convection with uext + delta u (perturbation)
    //     // (no u_induced)
    //     UVLM::Wake::Discretised::convect(zeta_star,
    //                                      uext_star_total,
    //                                      options.dt);
    //     // displace both zeta and gamma
    //     UVLM::Wake::General::displace_VecMat(gamma_star);
    //     UVLM::Wake::General::displace_VecVecMat(zeta_star);
    //
    //     // copy last row of zeta into zeta_star
    //     UVLM::Wake::Discretised::generate_new_row
    //     (
    //         zeta_star,
    //         zeta
    //     );
    // } else if (options.convection_scheme == 3)
    // {
    //     UVLM::Types::VecVecMatrixX uext_star_total;
    //     UVLM::Types::allocate_VecVecMat(uext_star_total, uext_star);
    //     UVLM::Types::VecVecMatrixX zeros;
    //     UVLM::Types::allocate_VecVecMat(zeros, uext_star);
    //     // total stream velocity
    //     UVLM::Unsteady::compute_resultant_grid_velocity
    //     (
    //         zeta_star,
    //         zeros,
    //         uext_star,
    //         rbm_velocity,
    //         uext_star_total
    //     );
    //     // convection with uext + delta u (perturbation) + u_ind
    //     UVLM::Types::VecVecMatrixX u_convection;
    //     UVLM::Types::allocate_VecVecMat
    //     (
    //         u_convection,
    //         uext_star_total
    //     );
    //     // induced velocity by vortex rings
    //     UVLM::BiotSavart::total_induced_velocity_on_wake
    //     (
    //         zeta,
    //         zeta_star,
    //         gamma,
    //         gamma_star,
    //         u_convection
    //     );
    //     // u_convection = u_convection + uext_star;
    //     for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    //     {
    //         for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
    //         {
    //             u_convection[i_surf][i_dim] = u_convection[i_surf][i_dim] + uext_star_total[i_surf][i_dim];
    //         }
    //     }
    //
    //     UVLM::Wake::Discretised::convect(zeta_star,
    //                                      u_convection,
    //                                      options.dt);
    //     // displace both zeta and gamma
    //     UVLM::Wake::General::displace_VecMat(gamma_star);
    //     UVLM::Wake::General::displace_VecVecMat(zeta_star);
    //
    //     // copy last row of zeta into zeta_star
    //     UVLM::Wake::Discretised::generate_new_row
    //     (
    //         zeta_star,
    //         zeta
    //     );
    // } else
    // {
    //     std::cerr << "convection_scheme == "
    //               << options.convection_scheme
    //               << " is not supported by the UVLM solver. \n"
    //               << "Supported options are from [0->3]"
    //               << std::endl;
    // }


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
    // UVLM::Types::
    // static:
    UVLM::PostProc::calculate_static_forces
    (
        zeta,
        zeta_star,
        gamma,
        gamma_star,
        uext,
        forces,
        steady_options,
        flightconditions
    );
    // dynamic::
    if (i_iter > 0)
    {
        // calculate dynamic forces
        UVLM::PostProc::calculate_dynamic_forces
        (
            zeta,
            zeta_star,
            zeta_col,
            gamma,
            gamma_star,
            previous_gamma,
            uext_total,
            normals,
            dynamic_forces,
            options,
            flightconditions
        );
    }
}


template <typename t_zeta,
          typename t_zeta_dot,
          typename t_uext,
          typename t_rbm_velocity,
          typename t_uext_out>
void UVLM::Unsteady::compute_resultant_grid_velocity
(
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_rbm_velocity& rbm_velocity,
    t_uext_out& uext_out
)
{
    const uint n_surf = zeta.size();
    UVLM::Types::initialise_VecVecMat(uext_out);
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        const uint n_col = zeta[i_surf][0].cols();
        const uint n_row = zeta[i_surf][0].rows();
        for (uint i_col=0; i_col<n_col; ++i_col)
        {
            for (uint i_row=0; i_row<n_row; ++i_row)
            {
                UVLM::Types::Vector3 zeta_temp;
                zeta_temp << zeta[i_surf][0](i_row, i_col),
                             zeta[i_surf][1](i_row, i_col),
                             zeta[i_surf][2](i_row, i_col);
                UVLM::Types::Vector3 w_cross_zeta;
                w_cross_zeta =
                    rbm_velocity.template block<3,1> (3, 0).cross(zeta_temp);
                for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                {
                    uext_out[i_surf][i_dim](i_row, i_col) =
                                              uext[i_surf][i_dim](i_row, i_col)
                                            - w_cross_zeta(i_dim)
                                            - zeta_dot[i_surf][i_dim](i_row, i_col)
                                            - rbm_velocity(i_dim);
                }
            }
        }
    }
}

// wake convection
// the UVMopts flag convection_scheme determines how the
// wake is convected.
// convection_scheme == 0 => prescribed and fixed
// convection_scheme == 1 => prescribed following deformations of the wing
// convection_scheme == 2 => free, convection based on u_ext
// convection_scheme == 3 => free, convection based on u_ext and induced velocities.
template <typename t_zeta,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_uext,
          typename t_uext_star,
          typename t_rbm_velocity>
void UVLM::Unsteady::convect_unsteady_wake
(
    const UVLM::Types::UVMopts& options,
    const t_zeta& zeta,
    t_zeta_star& zeta_star,
    const t_gamma& gamma,
    t_gamma_star& gamma_star,
    const t_uext& uext,
    const t_uext_star& uext_star,
    const t_rbm_velocity& rbm_velocity
)
{
    const uint n_surf = options.NumSurfaces;

    if (options.convection_scheme == 0)
    {
        UVLM::Wake::General::displace_VecMat(gamma_star);
    } else if (options.convection_scheme == 1)
    {
        std::cerr << "convection_scheme == "
                  << options.convection_scheme
                  << " is not yet implemented in the UVLM solver"
                  << std::endl;
    } else if (options.convection_scheme == 2)
    {
        UVLM::Types::VecVecMatrixX uext_star_total;
        UVLM::Types::allocate_VecVecMat(uext_star_total, uext_star);
        UVLM::Types::VecVecMatrixX zeros;
        UVLM::Types::allocate_VecVecMat(zeros, uext_star);
        // total stream velocity
        UVLM::Unsteady::compute_resultant_grid_velocity
        (
            zeta_star,
            zeros,
            uext_star,
            rbm_velocity,
            uext_star_total
        );
        // convection with uext + delta u (perturbation)
        // (no u_induced)
        UVLM::Wake::Discretised::convect(zeta_star,
                                         uext_star_total,
                                         options.dt);
        // displace both zeta and gamma
        UVLM::Wake::General::displace_VecMat(gamma_star);
        UVLM::Wake::General::displace_VecVecMat(zeta_star);

        // copy last row of zeta into zeta_star
        UVLM::Wake::Discretised::generate_new_row
        (
            zeta_star,
            zeta
        );
    } else if (options.convection_scheme == 3)
    {
        UVLM::Types::VecVecMatrixX uext_star_total;
        UVLM::Types::allocate_VecVecMat(uext_star_total, uext_star);
        UVLM::Types::VecVecMatrixX zeros;
        UVLM::Types::allocate_VecVecMat(zeros, uext_star);
        // total stream velocity
        UVLM::Unsteady::compute_resultant_grid_velocity
        (
            zeta_star,
            zeros,
            uext_star,
            rbm_velocity,
            uext_star_total
        );
        // convection with uext + delta u (perturbation) + u_ind
        UVLM::Types::VecVecMatrixX u_convection;
        UVLM::Types::allocate_VecVecMat
        (
            u_convection,
            uext_star_total
        );
        // induced velocity by vortex rings
        UVLM::BiotSavart::total_induced_velocity_on_wake
        (
            zeta,
            zeta_star,
            gamma,
            gamma_star,
            u_convection
        );
        // u_convection = u_convection + uext_star;
        for (uint i_surf=0; i_surf<n_surf; ++i_surf)
        {
            for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
            {
                u_convection[i_surf][i_dim] = u_convection[i_surf][i_dim] + uext_star_total[i_surf][i_dim];
            }
        }

        UVLM::Wake::Discretised::convect(zeta_star,
                                         u_convection,
                                         options.dt);
        // displace both zeta and gamma
        UVLM::Wake::General::displace_VecMat(gamma_star);
        UVLM::Wake::General::displace_VecVecMat(zeta_star);

        // copy last row of zeta into zeta_star
        UVLM::Wake::Discretised::generate_new_row
        (
            zeta_star,
            zeta
        );
    } else
    {
        std::cerr << "convection_scheme == "
                  << options.convection_scheme
                  << " is not supported by the UVLM solver. \n"
                  << "Supported options are from [0->3]"
                  << std::endl;
    }
    return;
}
