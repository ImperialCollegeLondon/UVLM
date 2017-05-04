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
                  typename t_zeta_star_col,
                  typename t_uext_col,
                  typename t_normals,
                  typename t_aic>
        void AIC
        (
            const unsigned int& Ktotal,
            const t_zeta& zeta,
            const t_zeta_col& zeta_col,
            const t_zeta_star& zeta_star,
            const t_zeta_star_col& zeta_star_col,
            const t_uext_col& uext_col,
            const t_normals& normals,
            const UVLM::Types::VMopts& options,
            t_aic& aic
        );


        template <typename t_zeta_col,
                  typename t_zeta_star,
                  typename t_uext_col,
                  typename t_zeta_dot_col,
                  typename t_gamma_star,
                  typename t_normal>
        void RHS
        (
            const t_zeta_col& zeta_col,
            const t_zeta_star& zeta_star,
            const t_uext_col& uext_col,
            const t_zeta_dot_col& zeta_dot_col,
            const t_gamma_star& gamma_star,
            const t_normal& normal,
            const UVLM::Types::VMopts& options,
            UVLM::Types::VectorX& rhs,
            unsigned int& Ktotal
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
                  typename t_zeta_star_col>
        void reconstruct_gamma
        (
            const UVLM::Types::VectorX& gamma_flat,
            t_gamma& gamma,
            const t_zeta_col& zeta_col,
            const t_zeta_star_col& zeta_star_col,
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
          typename t_zeta_star_col,
          typename t_uext_col,
          typename t_normals,
          typename t_aic>
void UVLM::Matrix::AIC
(
    const unsigned int& Ktotal,
    const t_zeta& zeta,
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    const t_zeta_star_col& zeta_star_col,
    const t_uext_col& uext_col,
    const t_normals& normals,
    const UVLM::Types::VMopts& options,
    t_aic& aic
)
{
    const unsigned int n_surf = options.NumSurfaces;
    UVLM::Types::VecDimensions dimensions;
    UVLM::Types::generate_dimensions(zeta_col, dimensions);
    UVLM::Types::VecDimensions dimensions_star;
    UVLM::Types::generate_dimensions(zeta_star_col, dimensions_star);

    // fill up AIC
    unsigned int i_offset = 0;
    for (unsigned int icol_surf=0; icol_surf<n_surf; ++icol_surf)
    {
        unsigned int k_surf = dimensions[icol_surf].first*
                              dimensions[icol_surf].second;

        unsigned int ii_offset = 0;
        // SURFACE - SURFACE coeffs
        for (unsigned int ii_surf=0; ii_surf<n_surf; ++ii_surf)
        {
            unsigned int kk_surf = dimensions[ii_surf].first*
                                   dimensions[ii_surf].second;
            UVLM::Types::MatrixX dummy_gamma;
            dummy_gamma.setOnes(dimensions[ii_surf].first,
                                dimensions[ii_surf].second);

            UVLM::Types::Block block = aic.block(i_offset, ii_offset, k_surf, kk_surf);
            UVLM::BiotSavart::multisurface
            (
                zeta[ii_surf],
                dummy_gamma,
                zeta_col[icol_surf],
                block,
                true,
                options.ImageMethod,
                normals[icol_surf]
            );
            ii_offset += kk_surf;
        } // end of surface-surface coeffs
        i_offset += k_surf;
        // SURFACE-WAKE coeffs
        // TODO
    }
    
    std::ofstream file;
    file.open("file.dat");
    file << aic << std::endl;
    file.close();


}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta_col,
          typename t_zeta_star,
          typename t_uext_col,
          typename t_zeta_dot_col,
          typename t_gamma_star,
          typename t_normal>
void UVLM::Matrix::RHS
(
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    const t_uext_col& uext_col,
    const t_zeta_dot_col& zeta_dot_col,
    const t_gamma_star& gamma_star,
    const t_normal& normal,
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
    unsigned int n_surf = dimensions.size();
    offset.setZero(n_surf, 2);

    unsigned int counter = 0;
    for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
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
          typename t_zeta_star_col>
void UVLM::Matrix::reconstruct_gamma
(
    const UVLM::Types::VectorX& gamma_flat,
    t_gamma& gamma,
    const t_zeta_col& zeta_col,
    const t_zeta_star_col& zeta_star_col,
    const UVLM::Types::VMopts& options
)
{
    if (!options.Steady)
    {
        std::cerr << "Not implemented in matrix.h, line=" << __LINE__ << std::endl;
    }
    const unsigned int n_surf = zeta_col.size();
    UVLM::Types::VecDimensions dimensions;
    UVLM::Types::generate_dimensions(zeta_col, dimensions);
    UVLM::Types::VecDimensions dimensions_star;
    UVLM::Types::generate_dimensions(zeta_star_col, dimensions_star);

    unsigned int i_flat = 0;
    for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
    {
        for (unsigned int i=0; i<dimensions[i_surf].first; ++i)
        {
            for (unsigned int j=0; j<dimensions[i_surf].second; ++j)
            {
                gamma[i_surf](i, j) = gamma_flat(i_flat++);
            }
        }
        // TODO unsteady and wake contribution
    }

}
