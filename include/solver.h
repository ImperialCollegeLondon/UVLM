#pragma once

#include "types.h"
#include "triads.h"
#include "constants.h"
#include "mapping.h"
#include "geometry.h"
#include "biotsavart.h"
#include "matrix.h"

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
            const UVLM::Types::VMopts& options
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
    const UVLM::Types::VMopts& options
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
    std::cout << "gamma-------------------------" << std::endl;
    std::cout << gamma_flat<< std::endl;

    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(gamma_flat,
                                    gamma,
                                    zeta_col,
                                    zeta_star_col,
                                    options);
    // std::cout << gamma_flat << std::endl;

}
