#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "biotsavart.h"

#include <fstream>

namespace UVLM
{
    namespace Matrix
    {
        // DECLARATIONS
        template <typename t_zeta,
                  typename t_zeta_col,
                  typename t_zeta_star,
                //   typename t_zeta_star_col,
                  typename t_uext_col,
                  typename t_normals,
                  typename t_aic>
        void AIC
        (
            const uint& Ktotal,
            const t_zeta& zeta,
            const t_zeta_col& zeta_col,
            const t_zeta_star& zeta_star,
            // const t_zeta_star_col& zeta_star_col,
            const t_uext_col& uext_col,
            const t_normals& normals,
            const UVLM::Types::VMopts& options,
            const bool horseshoe,
            t_aic& aic
        );


        template <typename t_zeta_col,
                  typename t_zeta_star,
                  typename t_uext_col,
                //   typename t_zeta_dot_col,
                  typename t_gamma_star,
                  typename t_normal>
        void RHS
        (
            const t_zeta_col& zeta_col,
            const t_zeta_star& zeta_star,
            const t_uext_col& uext_col,
            // const t_zeta_dot_col& zeta_dot_col,
            const t_gamma_star& gamma_star,
            const t_normal& normal,
            const UVLM::Types::VMopts& options,
            UVLM::Types::VectorX& rhs,
            uint& Ktotal
        );


        void generate_assembly_offset
        (
            const UVLM::Types::VecDimensions& dimensions,
            const UVLM::Types::VecDimensions& dimensions_star,
            UVLM::Types::MatrixX& offset,
            const bool steady
        );


        template <typename t_gamma,
                  typename t_zeta_col,
                  typename t_zeta_star>
        void reconstruct_gamma
        (
            const UVLM::Types::VectorX& gamma_flat,
            t_gamma& gamma,
            const t_zeta_col& zeta_col,
            const t_zeta_star& zeta_star,
            const UVLM::Types::VMopts& options
        );
    }
}
// SOURCE CODE
/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta,
          typename t_zeta_col,
          typename t_zeta_star,
        //   typename t_zeta_star_col,
          typename t_uext_col,
          typename t_normals,
          typename t_aic>
void UVLM::Matrix::AIC
(
    const uint& Ktotal,
    const t_zeta& zeta,
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    // const t_zeta_star_col& zeta_star_col,
    const t_uext_col& uext_col,
    const t_normals& normals,
    const UVLM::Types::VMopts& options,
    const bool horseshoe,
    t_aic& aic
)
{
    const uint n_surf = options.NumSurfaces;
    UVLM::Types::VecDimensions dimensions;
    UVLM::Types::generate_dimensions(zeta_col, dimensions);
    UVLM::Types::VecDimensions dimensions_star;
    UVLM::Types::generate_dimensions(zeta_star, dimensions_star, -1);

    // fill up AIC
    uint i_offset = 0;
    for (uint icol_surf=0; icol_surf<n_surf; ++icol_surf)
    {
        uint k_surf = dimensions[icol_surf].first*
                              dimensions[icol_surf].second;

        uint ii_offset = 0;
        for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
        {
            // SURFACE - SURFACE coeffs
            uint kk_surf = dimensions[ii_surf].first*
                                   dimensions[ii_surf].second;
            UVLM::Types::MatrixX dummy_gamma;
            UVLM::Types::MatrixX dummy_gamma_star;
            UVLM::Types::Block block = aic.block(i_offset, ii_offset, k_surf, kk_surf);
            if (options.Steady)
            {
                // steady wake coefficients
                dummy_gamma.setOnes(dimensions[ii_surf].first,
                                    dimensions[ii_surf].second);
                dummy_gamma_star.setOnes(dimensions_star[ii_surf].first,
                                         dimensions_star[ii_surf].second);
                UVLM::BiotSavart::multisurface_steady_wake
                (
                    zeta[ii_surf],
                    zeta_star[ii_surf],
                    dummy_gamma,
                    dummy_gamma_star,
                    zeta_col[icol_surf],
                    horseshoe,
                    block,
                    options.ImageMethod,
                    normals[icol_surf]
                );
            } else // unsteady case/free wake/no rollup
            {
                std::cerr << "Not implemented in matrix.h, line "
                          << __LINE__
                          << std::endl;
            }

            ii_offset += kk_surf;
        }
        i_offset += k_surf;

    }
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta_col,
          typename t_zeta_star,
          typename t_uext_col,
        //   typename t_zeta_dot_col,
          typename t_gamma_star,
          typename t_normal>
void UVLM::Matrix::RHS
(
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    const t_uext_col& uext_col,
    // const t_zeta_dot_col& zeta_dot_col,
    const t_gamma_star& gamma_star,
    const t_normal& normal,
    const UVLM::Types::VMopts& options,
    UVLM::Types::VectorX& rhs,
    uint& Ktotal
)
{
    const uint n_surf = options.NumSurfaces;

    // normal wash
    UVLM::Types::VecVecMatrixX uinc;
    // uinc = uext_col - zeta_dot_col;
    uinc = uext_col;
    // for dynamic simulations, add the influence of
    // RBM, deformation velocty, and wake induced velocity (except first row).



    if (!options.Steady)
    {
        // RBM

        // zeta_dot due to deformation
        
        // contribution of the wake to the incident velocity at the bound panels
        UVLM::Types::VecVecMatrixX uout;
        UVLM::Types::allocate_VecVecMat(uout, zeta_col);
        // TODO
        std::cerr << "Not implemented, matrix.h " << __LINE__ << std::endl;
        // UVLM::BiotSavart::multimultisurface
        // (
        //     zeta_star,
        //     gamma_star,
        //     zeta_col,
        //     uout
        // );
        uinc += uout;
    }


    // size of rhs
    uint ii = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        uint M = uinc[i_surf][0].rows();
        uint N = uinc[i_surf][0].cols();

        ii += M*N;
    }
    Ktotal = ii;
    rhs.setZero(Ktotal);

    // filling up RHS
    ii = -1;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        uint M = uinc[i_surf][0].rows();
        uint N = uinc[i_surf][0].cols();

        for (uint i=0; i<M; ++i)
        {
            for (uint j=0; j<N; ++j)
            {
                ++ii;
                // dot product of uinc and normal
                rhs(ii) =
                -(
                    uinc[i_surf][0](i,j) * normal[i_surf][0](i,j) +
                    uinc[i_surf][1](i,j) * normal[i_surf][1](i,j) +
                    uinc[i_surf][2](i,j) * normal[i_surf][2](i,j)
                );
            }
        }
    }
}


void UVLM::Matrix::generate_assembly_offset
(
    const UVLM::Types::VecDimensions& dimensions,
    const UVLM::Types::VecDimensions& dimensions_star,
    UVLM::Types::MatrixX& offset,
    const bool steady
)
{
    uint n_surf = dimensions.size();
    offset.setZero(n_surf, 2);

    uint counter = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        offset(i_surf, 0) = counter;
        counter += dimensions[i_surf].first*
                   dimensions[i_surf].second;

        offset(i_surf, 1) = counter;
        if (!steady)
        {
            counter += dimensions_star[i_surf].first*
                       dimensions_star[i_surf].second;
        }
    }
}

template <typename t_gamma,
          typename t_zeta_col,
          typename t_zeta_star>
void UVLM::Matrix::reconstruct_gamma
(
    const UVLM::Types::VectorX& gamma_flat,
    t_gamma& gamma,
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    const UVLM::Types::VMopts& options
)
{
    if (!options.Steady)
    {
        std::cerr << "Not implemented in matrix.h, line=" << __LINE__ << std::endl;
    }
    const uint n_surf = zeta_col.size();
    UVLM::Types::VecDimensions dimensions;
    UVLM::Types::generate_dimensions(zeta_col, dimensions);
    UVLM::Types::VecDimensions dimensions_star;
    UVLM::Types::generate_dimensions(zeta_star, dimensions_star, -1);

    uint i_flat = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        for (uint i=0; i<dimensions[i_surf].first; ++i)
        {
            for (uint j=0; j<dimensions[i_surf].second; ++j)
            {
                gamma[i_surf](i, j) = gamma_flat(i_flat++);
            }
        }
        // TODO unsteady and wake contribution
    }

}
