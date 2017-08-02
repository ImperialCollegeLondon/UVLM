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
    namespace Steady
    {
        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_uext,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_forces>
        void solver
        (
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_forces& forces,
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

    }
}

// SOURCE CODE

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta,
          typename t_zeta_dot,
          typename t_uext,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_forces>
void UVLM::Steady::solver
(
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_forces& forces,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // Generate collocation points info
    //  Declaration
    UVLM::Types::VecVecMatrixX zeta_col;
    // UVLM::Types::VecVecMatrixX zeta_dot_col;
    UVLM::Types::VecVecMatrixX uext_col;
    // UVLM::Types::VecVecMatrixX zeta_star_col;

    //  Allocation and mapping
    UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
    // UVLM::Geometry::generate_colocationMesh(zeta_dot, zeta_dot_col);
    UVLM::Geometry::generate_colocationMesh(uext, uext_col);
    // UVLM::Geometry::generate_colocationMesh(zeta_star, zeta_star_col);

    // panel normals
    UVLM::Types::VecVecMatrixX normals;
    UVLM::Types::allocate_VecVecMat(normals, zeta_col);
    UVLM::Geometry::generate_surfaceNormal(zeta, normals);

    // solve the steady horseshoe problem
    solve_horseshoe
    (
        zeta,
        zeta_col,
        uext_col,
        zeta_star,
        gamma,
        gamma_star,
        normals,
        options,
        flightconditions
    );

    // if options.horseshoe, it is finished.
    if (options.horseshoe)
    {
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
        return;
    }

    // if not, the wake has to be transformed into a normal, non-horseshoe
    // one:
    double delta_x = 1.0;

    UVLM::Wake::Horseshoe::to_wake(zeta_star,
                                   gamma_star,
                                   delta_x);
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

    // RHS generation
    UVLM::Types::VectorX rhs;
    unsigned int Ktotal;
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
                      aic);

    UVLM::Types::VectorX gamma_flat;
    gamma_flat = aic.partialPivLu().solve(rhs);

    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(gamma_flat,
                                    gamma,
                                    zeta_col,
                                    zeta_star,
                                    options);

    // copy gamma from trailing edge to wake if steady solution
    UVLM::Wake::Horseshoe::circulation_transfer(gamma,
                                                gamma_star);

}
