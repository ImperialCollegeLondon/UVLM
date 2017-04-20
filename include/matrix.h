#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "biotsavart.h"

namespace UVLM
{
    namespace Matrix
    {
        // DECLARATIONS
        template <typename t_zeta,
                  typename t_zeta_star,
                  typename t_normals,
                  typename t_aic>
        void AIC
        (
            t_zeta& zeta,
            t_zeta_star& zeta_star,
            t_normals& normals,
            UVLM::Types::VMopts& options,
            t_aic& aic
        );


        template <typename t_zeta_col,
                  typename t_zeta_star,
                  typename t_uext_col,
                  typename t_zeta_dot_col>
        void RHS
        (
            const t_zeta_col& zeta_col,
            const t_zeta_star& zeta_star,
            const t_uext_col& uext_col,
            const t_zeta_dot_col& zeta_dot_col,
            const UVLM::Types::VMopts& options,
            UVLM::Types::VectorX& rhs,
            unsigned int& Ktotal
        );
    }
}
// SOURCE CODE
/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta,
          typename t_zeta_star,
          typename t_normals,
          typename t_aic>
void UVLM::Matrix::AIC
(
    t_zeta& zeta,
    t_zeta_star& zeta_star,
    t_normals& normals,
    UVLM::Types::VMopts& options,
    t_aic& aic
)
{
    const unsigned int Ktotal = aic.rows();
    const unsigned int n_surf = options.NumSurfaces;


}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta_col,
          typename t_zeta_star,
          typename t_uext_col,
          typename t_zeta_dot_col>
void UVLM::Matrix::RHS
(
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    const t_uext_col& uext_col,
    const t_zeta_dot_col& zeta_dot_col,
    const UVLM::Types::VMopts& options,
    UVLM::Types::VectorX& rhs,
    unsigned int& Ktotal
)
{
    const unsigned int n_surf = options.NumSurfaces;

    // normal wash
    UVLM::Types::VecVecMatrixX uinc;
    uinc = uext_col - zeta_dot_col;

    // contribution of the wake to the incident velocity at the bound panels
    if (!options.Steady)
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
