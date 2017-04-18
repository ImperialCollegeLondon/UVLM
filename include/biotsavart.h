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
                  typename t_uout>
        void multimultisurface
        (
            const t_zeta&       zeta,
            const t_gamma&      gamma,
            const t_tsurface&   target_surface,
            t_uout&             uout,
            const bool&         image_method = false,
            const UVLM::Types::Real vortex_radius = 1e-5
        );

        template <typename t_zeta,
                  typename t_gamma,
                  typename t_ttriad,
                  typename t_uout>
        void multisurface
        (
            const t_zeta&       zeta,
            const t_gamma&      gamma,
            const t_ttriad&     target_triad,
            t_uout&       uout,
            const bool&         image_method = false,
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
            const UVLM::Types::Real vortex_radius = 1e-5
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
            const UVLM::Types::Real vortex_radius = 1e-5
        );

        template <typename t_triad>
        void segment
        (
            const t_triad& target_triad,
            const UVLM::Types::Vector3& v1,
            const UVLM::Types::Vector3& v2,
            const UVLM::Types::Real& gamma,
            UVLM::Types::Vector3& uind,
            const UVLM::Types::Real vortex_radius = 1e-5
        );
    }
}



// SOURCE CODE
template <typename t_triad>
void UVLM::BiotSavart::segment(const t_triad& target_triad,
             const UVLM::Types::Vector3& v1,
             const UVLM::Types::Vector3& v2,
             const UVLM::Types::Real& gamma,
             UVLM::Types::Vector3& uind,
             const UVLM::Types::Real vortex_radius
         )
{
    UVLM::Types::Vector3 r1 = target_triad - v1;
    UVLM::Types::Vector3 r2 = target_triad - v2;

    // Vortex core
    if ((r1.norm() < vortex_radius) || (r2.norm() < vortex_radius))
    {
        return;
    }

    UVLM::Types::Vector3 rx = r1.cross(r2);
    UVLM::Types::Real K = gamma/(UVLM::Constants::PI4*rx.squaredNorm())*
                                (target_triad.dot(r1)/r1.norm() -
                                 target_triad.dot(r2)/r2.norm());
    uind += K*rx;
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
        // return;
    }

    const unsigned int n_segment = 4;
    for (unsigned int i_segment=0; i_segment<n_segment; ++i_segment)
    {
        UVLM::Types::Vector3 v1;
        UVLM::Types::Vector3 v2;

        v1 << x(UVLM::Mapping::vortex_indices(0, 0),
                UVLM::Mapping::vortex_indices(0, 1)),
              y(UVLM::Mapping::vortex_indices(0, 0),
                UVLM::Mapping::vortex_indices(0, 1)),
              z(UVLM::Mapping::vortex_indices(0, 0),
                UVLM::Mapping::vortex_indices(0, 1));
        v2 << x(UVLM::Mapping::vortex_indices(1, 0),
                UVLM::Mapping::vortex_indices(1, 1)),
              y(UVLM::Mapping::vortex_indices(1, 0),
                UVLM::Mapping::vortex_indices(1, 1)),
              z(UVLM::Mapping::vortex_indices(1, 0),
                UVLM::Mapping::vortex_indices(1, 1));

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

    for (unsigned int i=Mstart; i<Mend; ++i)
    {
        for (unsigned int j=Nstart; j<Nend; ++j)
        {
            UVLM::BiotSavart::vortex_ring(target_triad,
                                          zeta[0].template block<2, 2>(i,j),
                                          zeta[1].template block<2, 2>(i,j),
                                          zeta[2].template block<2, 2>(i,j),
                                          gamma(i,j),
                                          uout);
        }
    }
}



template <typename t_zeta,
          typename t_gamma,
          typename t_tsurface,
          typename t_uout>
void UVLM::BiotSavart::multisurface
(
    const t_zeta&       zeta,
    const t_gamma&      gamma,
    const t_tsurface&   target_surface,
    t_uout&             uout,
    const bool&         image_method,
    const UVLM::Types::Real vortex_radius
)
{
    unsigned int n_surf = zeta.size();

    unsigned int i_col = target_surface[0].rows();
    unsigned int j_col = target_surface[0].cols();
    UVLM::Types::Vector3 u_temp;
    for (unsigned int i=0; i<i_col; ++i)
    {
        for (unsigned int j=0; j<j_col; ++j)
        {
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                u_temp.setZero();
                UVLM::Types::Vector3 triad;
                triad << target_surface[0](i,j),
                         target_surface[1](i,j),
                         target_surface[2](i,j);
                UVLM::BiotSavart::surface(zeta[i_surf],
                                          gamma[i_surf],
                                          triad,
                                          u_temp);
                uout[0](i, j) = u_temp(0);
                uout[1](i, j) = u_temp(1);
                uout[2](i, j) = u_temp(2);
            }
        }
    }
}



template <typename t_zeta,
          typename t_gamma,
          typename t_tsurface,
          typename t_uout>
void UVLM::BiotSavart::multimultisurface
(
    const t_zeta&       zeta,
    const t_gamma&      gamma,
    const t_tsurface&   target_surface,
    t_uout&             uout,
    const bool&         image_method,
    const UVLM::Types::Real vortex_radius
)
{
    unsigned int n_surf = target_surface.size();

    // #pragma omp parallel for
    for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
    {
        UVLM::BiotSavart::multisurface(zeta,
                                       gamma,
                                       target_surface[i_surf],
                                       uout[i_surf]);
    }
}
