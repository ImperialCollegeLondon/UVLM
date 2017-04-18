#pragma once

#include "types.h"
#include "triads.h"
#include "constants.h"
#include "mapping.h"
#include "geometry.h"
#include "biotsavart.h"

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
                  typename t_gamma_star>
        void solve
        (
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            UVLM::Types::VMopts& VMOPTS
        );

        template <typename t_in,
                  typename t_out>
        void generate_colocationMesh
        (
            t_in& vortex_mesh,
            t_out& collocation_mesh
        );

        template <typename t_uinc,
                  typename t_normal>
        void generate_RHS
        (
            const t_uinc& uinc,
            const t_normal& normal,
            UVLM::Types::VectorX& rhs,
            unsigned int& Ktotal
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
    for (int i_surf=0; i_surf<dimensions.size(); ++i_surf)
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
          typename t_gamma_star>
void UVLM::Solver::solve
(
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    UVLM::Types::VMopts& VMOPTS
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

    // panel normal
    UVLM::Types::VecVecMatrixX normal;
    UVLM::Types::allocate_VecVecMat(normal, zeta_col);
    UVLM::Geometry::generate_surfaceNormal(zeta, normal);

    // normal wash
    UVLM::Types::VecVecMatrixX uinc;
    uinc = uext_col - zeta_dot_col;

    // contribution of the wake to the incident velocity at the bound panels
    if (!VMOPTS.Steady)
    {
        UVLM::Types::VecVecMatrixX uout;
        UVLM::Types::allocate_VecVecMat(uout, zeta_col);
        UVLM::BiotSavart::multimultisurface
        (
            zeta_star,
            gamma_star,
            zeta_col,
            uout
        );
        uinc += uout;
    }

    // RHS generation
    UVLM::Types::VectorX rhs;
    unsigned int Ktotal;
    UVLM::Solver::generate_RHS(uinc, normal, rhs, Ktotal);

    // AIC generation
    UVLM::Types::MatrixX aic = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
    UVLM::Solver::generate_AIC()




}


/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_uinc,
          typename t_normal>
void UVLM::Solver::generate_RHS
(
    const t_uinc& uinc,
    const t_normal& normal,
    UVLM::Types::VectorX& rhs,
    unsigned int& Ktotal
)
{
    unsigned int n_surf = uinc.size();

    // size of rhs
    unsigned int ii = 0;
    for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
    {
        unsigned int M = uinc[i_surf][0].rows();
        unsigned int N = uinc[i_surf][0].cols();

        ii += M*N;
    }
    Ktotal = ii;
    rhs.setZero(Ktotal);

    // filling up RHS
    ii = -1;
    for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
    {
        unsigned int M = uinc[i_surf][0].rows();
        unsigned int N = uinc[i_surf][0].cols();

        for (unsigned int i=0; i<M; ++i)
        {
            for (unsigned int j=0; j<N; ++j)
            {
                // dot product of uinc and normal
                rhs(++ii) =
                (
                    uinc[i_surf][0](i,j) + normal[i_surf][0](i,j) +
                    uinc[i_surf][1](i,j) + normal[i_surf][1](i,j) +
                    uinc[i_surf][2](i,j) + normal[i_surf][2](i,j)
                );
            }
        }
    }
}
