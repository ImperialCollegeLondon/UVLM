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
            const uint& Ktotal
        );


        void generate_assembly_offset
        (
            const UVLM::Types::VecDimensions& dimensions,
            const UVLM::Types::VecDimensions& dimensions_star,
            UVLM::Types::MatrixX& offset,
            const bool steady
        );


        template <typename t_gamma,
                  typename t_zeta_col>
        void reconstruct_gamma
        (
            const UVLM::Types::VectorX& gamma_flat,
            t_gamma& gamma,
            const t_zeta_col& zeta_col
        );


        template <typename t_gamma,
                  typename t_zeta_col>
        void deconstruct_gamma
        (
            const t_gamma& gamma,
            UVLM::Types::VectorX& gamma_flat,
            const t_zeta_col& zeta_col
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
    UVLM::Types::generate_dimensions(zeta, dimensions, - 1);

    UVLM::Types::VecDimensions dimensions_star;
    UVLM::Types::generate_dimensions(zeta_star, dimensions_star, -1);

    // build the offsets beforehand
    // (parallel variation)
    std::vector<uint> offset;
    uint i_offset = 0;
    for (uint icol_surf=0; icol_surf<n_surf; ++icol_surf)
    {
        offset.push_back(i_offset);
        uint k_surf = dimensions[icol_surf].first*
                      dimensions[icol_surf].second;
        i_offset += k_surf;
    }

    // fill up AIC
    for (uint icol_surf=0; icol_surf<n_surf; ++icol_surf)
    {
        uint k_surf = dimensions[icol_surf].first*
                      dimensions[icol_surf].second;

        // uint ii_offset = 0;
        for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
        {
            uint kk_surf = dimensions[ii_surf].first*
                           dimensions[ii_surf].second;
            UVLM::Types::MatrixX dummy_gamma;
            UVLM::Types::MatrixX dummy_gamma_star;
            UVLM::Types::Block block = aic.block(offset[icol_surf], offset[ii_surf], k_surf, kk_surf);
            // steady wake coefficients
            dummy_gamma.setOnes(dimensions[ii_surf].first,
                                dimensions[ii_surf].second);
            if (options.Steady)
            {
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
            } else // unsteady case
            {
                dummy_gamma_star.setOnes(1,
                                         dimensions_star[ii_surf].second);
                UVLM::BiotSavart::multisurface_unsteady_wake
                (
                    zeta[ii_surf],
                    zeta_star[ii_surf],
                    dummy_gamma,
                    dummy_gamma_star,
                    zeta_col[icol_surf],
                    block,
                    options.ImageMethod,
                    normals[icol_surf],
                    1
                );
            }
        }
    }
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta_col,
          typename t_zeta_star,
          typename t_uext_col,
          typename t_gamma_star,
          typename t_normal>
void UVLM::Matrix::RHS
(
    const t_zeta_col& zeta_col,
    const t_zeta_star& zeta_star,
    const t_uext_col& uinc_col,
    const t_gamma_star& gamma_star,
    const t_normal& normal,
    const UVLM::Types::VMopts& options,
    UVLM::Types::VectorX& rhs,
    const uint& Ktotal
)
{
    const uint n_surf = options.NumSurfaces;

    // make a copy of uinc in order to add wake effects
    UVLM::Types::VecVecMatrixX u_col;
    UVLM::Types::allocate_VecVecMat(u_col, uinc_col);
    UVLM::Types::copy_VecVecMat(uinc_col, u_col);

    rhs.setZero(Ktotal);

    // filling up RHS
    int ii = -1;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        uint M = uinc_col[i_surf][0].rows();
        uint N = uinc_col[i_surf][0].cols();


        for (uint i=0; i<M; ++i)
        {
            for (uint j=0; j<N; ++j)
            {
                if (!options.Steady)
                {
                    // we have to add the wake effect on the induced velocity.
                    UVLM::Types::Vector3 collocation_coords;
                    collocation_coords << zeta_col[i_surf][0](i,j),
                                          zeta_col[i_surf][1](i,j),
                                          zeta_col[i_surf][2](i,j);
                    for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
                    {
                        UVLM::Types::VecMatrixX induced_vel;
                        UVLM::Types::allocate_VecMat(induced_vel, zeta_star[ii_surf], -1);

                        UVLM::BiotSavart::surface
                        (
                            zeta_star[ii_surf],
                            gamma_star[ii_surf],
                            collocation_coords,
                            induced_vel,
                            1
                        );
                        u_col[i_surf][0](i, j) += induced_vel[0].sum();
                        u_col[i_surf][1](i, j) += induced_vel[1].sum();
                        u_col[i_surf][2](i, j) += induced_vel[2].sum();
                    }
                }
                // dot product of uinc and panel normal
                rhs(++ii) =
                -(
                    u_col[i_surf][0](i,j)*normal[i_surf][0](i,j) +
                    u_col[i_surf][1](i,j)*normal[i_surf][1](i,j) +
                    u_col[i_surf][2](i,j)*normal[i_surf][2](i,j)
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
    std::cerr << "Not sure this is correct, dont use (or debug)!!!" << std::endl;
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
          typename t_zeta_col>
void UVLM::Matrix::reconstruct_gamma
(
    const UVLM::Types::VectorX& gamma_flat,
    t_gamma& gamma,
    const t_zeta_col& zeta_col
)
{
    const uint n_surf = gamma.size();
    UVLM::Types::VecDimensions dimensions;
    UVLM::Types::generate_dimensions(zeta_col, dimensions);

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
    }

}

template <typename t_gamma,
          typename t_zeta_col>
void UVLM::Matrix::deconstruct_gamma
(
    const t_gamma& gamma,
    UVLM::Types::VectorX& gamma_flat,
    const t_zeta_col& zeta_col
)
{
    const uint n_surf = gamma.size();
    UVLM::Types::VecDimensions dimensions;
    UVLM::Types::generate_dimensions(zeta_col, dimensions);

    uint n_total = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        n_total += dimensions[i_surf].first*dimensions[i_surf].second;
    }
    gamma_flat.resize(n_total);

    uint i_flat = 0;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        for (uint i=0; i<dimensions[i_surf].first; ++i)
        {
            for (uint j=0; j<dimensions[i_surf].second; ++j)
            {
                gamma_flat(i_flat++) = gamma[i_surf](i, j);
            }
        }
    }

}
