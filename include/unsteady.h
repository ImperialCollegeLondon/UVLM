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
#include "steady.h"

#include <iostream>

// DECLARATIONS
namespace UVLM
{
    namespace Unsteady
    {
        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_uext,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_rbm_velocity,
                  typename t_forces>
        void solver
        (
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_rbm_velocity& rbm_velocity,
            t_forces& forces,
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
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_rbm_velocity,
          typename t_forces>
void UVLM::Unsteady::solver
(
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_rbm_velocity& rbm_velocity,
    t_forces& forces,
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


    // wake convection
    // the UVMopts flag convection_scheme determines how the
    // wake is convected.
    // convection_scheme == 0 => prescribed and fixed
    // convection_scheme == 1 => prescribed following deformations of the wing
    // convection_scheme == 2 => free, convection based on u_ext
    // convection_scheme == 3 => free, convection based on u_ext and induced velocities.
    if (options.convection_scheme == 0)
    {
        UVLM::Wake::General::displace_VecMat(gamma_star);
        // set first row's gamma_star to 0
        for (uint i_surf=0; i_surf<n_surf; ++i_surf)
        {
            gamma_star[i_surf].template row(0).setZero();
        }
    } else
    {
        std::cerr << "convection_scheme == "
                  << options.convection_scheme
                  << " is not yet supported by the UVLM solver"
                  << std::endl;
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

    // forces calculation
    // UVLM::Types::
    // static:

    // dynamic::
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
