#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "mapping.h"

#include <limits>

namespace UVLM
{
    namespace BiotSavart
    {
        // DECLARATIONS
        template <typename t_zeta,
                  typename t_gamma,
                  typename t_tsurface,
                  typename t_uout,
                  typename t_normals>
        void multisurface
        (
            const t_zeta&       zeta,
            const t_gamma&      gamma,
            const t_tsurface&   target_surface,
            t_uout&             uout,
            const bool&         reduction = false,
            const bool&         image_method = false,
            const t_normals&    normal = NULL,
            const UVLM::Types::Real vortex_radius = 1e-5
        );

        template <typename t_zeta,
                  typename t_gamma,
                  typename t_ttriad,
                  typename t_uout>
        void surface
        (
            const t_zeta&       zeta,
            const t_gamma&      gamma,
            const t_ttriad&     target_triad,
            t_uout&             uout,
            unsigned int        Mstart = 0,
            unsigned int        Nstart = 0,
            unsigned int        Mend   = -1,
            unsigned int        Nend   = -1,
            const bool&         image_method = false,
            const UVLM::Types::Real vortex_radius = 1e-8
        );

        template <typename t_triad,
                  typename t_block>
        void vortex_ring
        (
            const t_triad& target_triad,
            const t_block& x,
            const t_block& y,
            const t_block& z,
            const UVLM::Types::Real& gamma,
            UVLM::Types::Vector3& uind,
            const UVLM::Types::Real vortex_radius = 1e-8
        );

        template <typename t_triad>
        void segment
        (
            const t_triad& target_triad,
            const UVLM::Types::Vector3& v1,
            const UVLM::Types::Vector3& v2,
            const UVLM::Types::Real& gamma,
            UVLM::Types::Vector3& uind,
            const UVLM::Types::Real vortex_radius = 1e-8
        );
    }
}



// SOURCE CODE
template <typename t_triad>
void UVLM::BiotSavart::segment
        (
            const t_triad& rp,
            const UVLM::Types::Vector3& v1,
            const UVLM::Types::Vector3& v2,
            const UVLM::Types::Real& gamma,
            UVLM::Types::Vector3& uind,
            const UVLM::Types::Real vortex_radius
        )
{
    UVLM::Types::Vector3 r0 = v2 - v1;
    UVLM::Types::Vector3 r1 = rp - v1;
    UVLM::Types::Vector3 r2 = rp - v2;
    UVLM::Types::Vector3 r1_cross_r2 = r1.cross(r2);
    UVLM::Types::Real r1_cross_r2_mod_sq = r1_cross_r2.squaredNorm();

    UVLM::Types::Real r1_mod = r1.norm();
    UVLM::Types::Real r2_mod = r2.norm();

    if (r1_mod < vortex_radius ||
        r2_mod < vortex_radius ||
        r1_cross_r2_mod_sq < vortex_radius)
    {
        std::cout << "rp = " << rp.transpose() << std::endl;
        std::cout << "v1 = " << v1.transpose() << std::endl;
        std::cout << "v2 = " << v2.transpose() << std::endl;
        std::cout << "r1 = " << r1.transpose() << std::endl;
        std::cout << "r2 = " << r2.transpose() << std::endl;
        std::cout << "r1_mod = " << r1_mod << std::endl;
        std::cout << "r2_mod = " << r2_mod << std::endl;
        std::cout << "r1_cross_r2_mod_sq = " << r1_cross_r2_mod_sq << std::endl;
        std::cerr << "In core" << std::endl;
        return;
    }

    UVLM::Types::Real r0_dot_r1 = r0.dot(r1);
    UVLM::Types::Real r0_dot_r2 = r0.dot(r2);

    UVLM::Types::Real K;
    K = (gamma/(UVLM::Constants::PI4*r1_cross_r2_mod_sq))*
        (r0_dot_r1/r1_mod - r0_dot_r2/r2_mod);

    uind += K*r1_cross_r2;
}


template <typename t_triad,
          typename t_block>
void UVLM::BiotSavart::vortex_ring
(
    const t_triad& target_triad,
    const t_block& x,
    const t_block& y,
    const t_block& z,
    const UVLM::Types::Real& gamma,
    UVLM::Types::Vector3& uind,
    const UVLM::Types::Real vortex_radius
)
{
    if (std::abs(gamma) < UVLM::Constants::EPSILON)
    {
        std::cerr << "Gamma almost 0 " << std::endl;
        return;
    }

    UVLM::Types::Vector3 v1;
    UVLM::Types::Vector3 v2;
    const unsigned int n_segment = 4;
    for (unsigned int i_segment=0; i_segment<n_segment; ++i_segment)
    {
        unsigned int start = i_segment;
        unsigned int end = (start + 1)%n_segment;

        v1 << x(UVLM::Mapping::vortex_indices(start, 0),
                UVLM::Mapping::vortex_indices(start, 1)),
              y(UVLM::Mapping::vortex_indices(start, 0),
                UVLM::Mapping::vortex_indices(start, 1)),
              z(UVLM::Mapping::vortex_indices(start, 0),
                UVLM::Mapping::vortex_indices(start, 1));
        v2 << x(UVLM::Mapping::vortex_indices(end, 0),
                UVLM::Mapping::vortex_indices(end, 1)),
              y(UVLM::Mapping::vortex_indices(end, 0),
                UVLM::Mapping::vortex_indices(end, 1)),
              z(UVLM::Mapping::vortex_indices(end, 0),
                UVLM::Mapping::vortex_indices(end, 1));

        UVLM::BiotSavart::segment(target_triad,
                                  v1,
                                  v2,
                                  gamma,
                                  uind);
    }
}



template <typename t_zeta,
          typename t_gamma,
          typename t_ttriad,
          typename t_uout>
void UVLM::BiotSavart::surface
(
    const t_zeta&       zeta,
    const t_gamma&      gamma,
    const t_ttriad&     target_triad,
    t_uout&             uout,
    unsigned int        Mstart,
    unsigned int        Nstart,
    unsigned int        Mend,
    unsigned int        Nend,
    const bool&         image_method,
    const UVLM::Types::Real vortex_radius
)
{
    // If Mend or Nend are == -1, their values are taken as the surface M and N
    if (Mend == -1) {Mend = gamma.rows();}
    if (Nend == -1) {Nend = gamma.cols();}

    UVLM::Types::Vector3 temp_uout;
    for (unsigned int i=Mstart; i<Mend; ++i)
    {
        for (unsigned int j=Nstart; j<Nend; ++j)
        {
            temp_uout.setZero();
            UVLM::BiotSavart::vortex_ring(target_triad,
                                          zeta[0].template block<2, 2>(i,j),
                                          zeta[1].template block<2, 2>(i,j),
                                          zeta[2].template block<2, 2>(i,j),
                                          gamma(i,j),
                                          temp_uout);
            uout[0](i, j) = temp_uout(0);
            uout[1](i, j) = temp_uout(1);
            uout[2](i, j) = temp_uout(2);
        }
    }
}



template <typename t_zeta,
          typename t_gamma,
          typename t_tsurface,
          typename t_uout,
          typename t_normals>
void UVLM::BiotSavart::multisurface
(
    const t_zeta&       zeta,
    const t_gamma&      gamma,
    const t_tsurface&   target_surface,
    t_uout&             uout,
    const bool&         reduction,
    const bool&         image_method,
    const t_normals&    normal,
    const UVLM::Types::Real vortex_radius
)
{
    const unsigned int rows_collocation = target_surface[0].rows();
    const unsigned int cols_collocation = target_surface[0].cols();

    UVLM::Types::VecMatrixX temp_uout;
    // UVLM::Types::allocate_VecMat(temp_uout, target_surface);
    UVLM::Types::allocate_VecMat(temp_uout, zeta, -1);

    if (!reduction)
    {
        std::cerr << "Not implemented, biotsavart.h, line "
                  << __LINE__
                  << std::endl;
    }
    int collocation_counter = -1;
    int surface_counter;
    for (unsigned int i_col=0; i_col<rows_collocation; ++i_col)
    {
        for (unsigned int j_col=0; j_col<cols_collocation; ++j_col)
        {
            ++collocation_counter;
            UVLM::Types::initialise_VecMat(temp_uout, 0.0);
            UVLM::Types::Vector3 target_triad;
            target_triad << target_surface[0](i_col, j_col),
                            target_surface[1](i_col, j_col),
                            target_surface[2](i_col, j_col);
            UVLM::BiotSavart::surface(zeta,
                                      gamma,
                                      target_triad,
                                      temp_uout);
            unsigned int surf_rows = gamma.rows();
            unsigned int surf_cols = gamma.cols();

            surface_counter = -1;
            for (unsigned int i_surf=0; i_surf<surf_rows; ++i_surf)
            {
                for (unsigned int j_surf=0; j_surf<surf_cols; ++j_surf)
                {
                    ++surface_counter;
                    uout(collocation_counter, surface_counter) =
                        temp_uout[0](i_surf, j_surf)*normal[0](i_col, j_col) +
                        temp_uout[1](i_surf, j_surf)*normal[1](i_col, j_col) +
                        temp_uout[2](i_surf, j_surf)*normal[2](i_col, j_col);
                }
            }
        }
    }
}
