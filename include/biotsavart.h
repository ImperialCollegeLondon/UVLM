#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "mapping.h"

#include <limits>
#include <math.h>
#include <cmath>

#define VORTEX_RADIUS 1e-2

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
            const UVLM::Types::IntPair& dimensions,
            const bool&         image_method = false,
            const t_normals&    normal = NULL,
            const bool&         horseshoe = false,
            const UVLM::Types::Real vortex_radius = VORTEX_RADIUS
        );

        template <typename t_zeta,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_tsurface,
                  typename t_uout,
                  typename t_normals>
        void multisurface_steady_wake
        (
            const t_zeta&       zeta,
            const t_zeta_star&  zeta_star,
            const t_gamma&      gamma,
            const t_gamma_star& gamma_star,
            const t_tsurface&   target_surface,
            const bool&         horseshoe,
            t_uout&             uout,
            const bool&         image_method = false,
            const t_normals&    normal = NULL,
            const UVLM::Types::Real vortex_radius = VORTEX_RADIUS
        );

        template <typename t_zeta,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_tsurface,
                  typename t_uout,
                  typename t_normals>
        void multisurface_unsteady_wake
        (
            const t_zeta&       zeta,
            const t_zeta_star&  zeta_star,
            const t_gamma&      gamma,
            const t_gamma_star& gamma_star,
            const t_tsurface&   target_surface,
            t_uout&             uout,
            const bool&         image_method,
            const t_normals&    normal,
            const int&          n_rows = -1
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
            const UVLM::Types::Real vortex_radius = VORTEX_RADIUS
        );

        template <typename t_zeta,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_ttriad,
                  typename t_uout>
        void surface_with_steady_wake
        (
            const t_zeta&       zeta,
            const t_zeta_star&  zeta_star,
            const t_gamma&      gamma,
            const t_gamma_star& gamma_star,
            const t_ttriad&     target_triad,
            const bool&         horseshoe,
            t_uout&             uout,
            const bool&         image_method = false,
            const UVLM::Types::Real vortex_radius = VORTEX_RADIUS
        );

        template <typename t_zeta,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_ttriad,
                  typename t_uout>
        void surface_with_unsteady_wake
        (
            const t_zeta&       zeta,
            const t_zeta_star&  zeta_star,
            const t_gamma&      gamma,
            const t_gamma_star& gamma_star,
            const t_ttriad&     target_triad,
            t_uout&             uout,
            const bool&         image_method,
            const int&          n_rows = -1 // default val = -1
        );

        template <typename t_triad,
                  typename t_block>
        void vortex_ring
        (
            const t_triad& target_triad,
            const t_block& x,
            const t_block& y,
            const t_block& z,
            const UVLM::Types::Real& gamma_star,
            UVLM::Types::Vector3& uind,
            const UVLM::Types::Real vortex_radius = VORTEX_RADIUS
        );

        template <typename t_triad,
                  typename t_block>
        void horseshoe
        (
            const t_triad& target_triad,
            const t_block& x,
            const t_block& y,
            const t_block& z,
            const UVLM::Types::Real& gamma,
            UVLM::Types::Vector3& uind,
            const UVLM::Types::Real vortex_radius = VORTEX_RADIUS
        );

        template <typename t_triad>
        void segment
        (
            const t_triad& target_triad,
            const UVLM::Types::Vector3& v1,
            const UVLM::Types::Vector3& v2,
            const UVLM::Types::Real& gamma,
            UVLM::Types::Vector3& uind,
            const UVLM::Types::Real vortex_radius = VORTEX_RADIUS
        );



        template <typename t_zeta,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_uout>
        void total_induced_velocity_on_wake
        (
            const t_zeta&       zeta,
            const t_zeta_star&  zeta_star,
            const t_gamma&      gamma,
            const t_gamma_star& gamma_star,
            t_uout&             uout,
            const bool&         image_method = false,
            const UVLM::Types::Real vortex_radius = VORTEX_RADIUS
        );



        template <typename t_zeta,
                  typename t_gamma,
                  typename t_zeta_col,
                  typename t_u_ind>
        void whole_surface_on_surface
        (
            const t_zeta& zeta,
            const t_gamma& gamma,
            const t_zeta_col& zeta_col,
            t_u_ind& u_ind,
            const bool image_method = false
        );


        template <typename t_zeta,
                  typename t_gamma,
                  typename t_ttriad,
                  typename t_uout>
        void whole_surface
        (
            const t_zeta&       zeta,
            const t_gamma&      gamma,
            const t_ttriad&     target_triad,
            t_uout&             uout,
            unsigned int        Mstart = 0,
            unsigned int        Nstart = 0,
            unsigned int        Mend = -1,
            unsigned int        Nend = -1,
            const bool&         image_method = false,
            const UVLM::Types::Real vortex_radius = VORTEX_RADIUS
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
        r1_cross_r2_mod_sq < vortex_radius*vortex_radius)
    {
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
void UVLM::BiotSavart::horseshoe
(
    const t_triad& target_triad,
    const t_block& x,
    const t_block& y,
    const t_block& z,
    const UVLM::Types::Real& gamma_star,
    UVLM::Types::Vector3& uind,
    const UVLM::Types::Real vortex_radius
)
{
    if (std::abs(gamma_star) < UVLM::Constants::EPSILON)
    {
        // std::cerr << "Gamma almost 0 " << std::endl;
        return;
    }
    // three segments.
    //
    //     0___________3
    //      |         |
    //      |         |
    //      |         |
    //      |         |
    //      |         |
    //      |         |
    //      |         |
    //     1|         |2
    //
    // segments 0-1 and 2-3 are represented as length 1, but they are effectively
    // infinite
    // segment 3-0 is considered as a normal one

    UVLM::Types::Vector3 v1;
    UVLM::Types::Vector3 v2;
    // segment 3-0
    v1 << x(UVLM::Mapping::vortex_indices(3, 0),
            UVLM::Mapping::vortex_indices(3, 1)),
          y(UVLM::Mapping::vortex_indices(3, 0),
            UVLM::Mapping::vortex_indices(3, 1)),
          z(UVLM::Mapping::vortex_indices(3, 0),
            UVLM::Mapping::vortex_indices(3, 1));
    v2 << x(UVLM::Mapping::vortex_indices(0, 0),
            UVLM::Mapping::vortex_indices(0, 1)),
          y(UVLM::Mapping::vortex_indices(0, 0),
            UVLM::Mapping::vortex_indices(0, 1)),
          z(UVLM::Mapping::vortex_indices(0, 0),
            UVLM::Mapping::vortex_indices(0, 1));
    UVLM::BiotSavart::segment(target_triad,
                              v1,
                              v2,
                              gamma_star,
                              uind);

    // segment 0-1
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

    // here the segment will be considered as 1----->2 and the
    // point 2 is in infinity, so beta2=pi
    UVLM::Types::Vector3 r0 = v2 - v1;
    UVLM::Types::Vector3 r1 = target_triad - v1;
    UVLM::Types::Vector3 r2 = target_triad - v2;
    UVLM::Types::Vector3 r1_cross_r2 = r1.cross(r2);
    UVLM::Types::Real dist = (r1_cross_r2).norm()/r0.norm();
    UVLM::Types::Real beta1;
    UVLM::Types::Real beta2;
    UVLM::Types::Vector3 u_radial;
    if (!((r1.norm() < vortex_radius) ||
         (r2.norm() < vortex_radius) ||
         (r1_cross_r2.norm() < vortex_radius)))
    {
        beta1 = r0.dot(r1)/(r0.norm()*r1.norm());
        beta2 = UVLM::Constants::PI;

        u_radial = (r1_cross_r2)/(r1_cross_r2).norm();

        uind += gamma_star/(UVLM::Constants::PI4*dist)*(beta1 + 1.0)*
                u_radial;
    }

    // segment 2-3
    v1 << x(UVLM::Mapping::vortex_indices(2, 0),
            UVLM::Mapping::vortex_indices(2, 1)),
          y(UVLM::Mapping::vortex_indices(2, 0),
            UVLM::Mapping::vortex_indices(2, 1)),
          z(UVLM::Mapping::vortex_indices(2, 0),
            UVLM::Mapping::vortex_indices(2, 1));
    v2 << x(UVLM::Mapping::vortex_indices(3, 0),
            UVLM::Mapping::vortex_indices(3, 1)),
          y(UVLM::Mapping::vortex_indices(3, 0),
            UVLM::Mapping::vortex_indices(3, 1)),
          z(UVLM::Mapping::vortex_indices(3, 0),
            UVLM::Mapping::vortex_indices(3, 1));

    // here the segment will be considered as 1----->2 and the
    // point 1 is in infinity, so beta1=0
    r0 = v2 - v1;
    r1 = target_triad - v1;
    r2 = target_triad - v2;
    r1_cross_r2 = r1.cross(r2);
    dist = (r1_cross_r2).norm()/r0.norm();
    if (!((r1.norm() < vortex_radius) ||
         (r2.norm() < vortex_radius) ||
         (r1_cross_r2.norm() < vortex_radius)))
    {
        beta2 = r0.dot(r2)/(r0.norm()*r2.norm());
        dist = (r1.cross(r2)).norm()/r0.norm();

        u_radial = (r1_cross_r2)/(r1_cross_r2).norm();

        uind += gamma_star/(UVLM::Constants::PI4*dist)*(1.0 - beta2)*
                u_radial;
    }
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
            uout[0](i, j) += temp_uout(0);
            uout[1](i, j) += temp_uout(1);
            uout[2](i, j) += temp_uout(2);
        }
    }
}


template <typename t_zeta,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_ttriad,
          typename t_uout>
void UVLM::BiotSavart::surface_with_steady_wake
(
    const t_zeta&       zeta,
    const t_zeta_star&  zeta_star,
    const t_gamma&      gamma,
    const t_gamma_star& gamma_star,
    const t_ttriad&     target_triad,
    const bool&         horseshoe,
    t_uout&             uout,
    const bool&         image_method,
    const UVLM::Types::Real vortex_radius
)
{
    const uint Mstart = 0;
    const uint Nstart = 0;
    const uint Mend = gamma.rows();
    const uint Nend = gamma.cols();

    UVLM::Types::Vector3 temp_uout;
    const uint ii = 0;
    for (unsigned int i=Mstart; i<Mend; ++i)
    {
        for (unsigned int j=Nstart; j<Nend; ++j)
        {
            temp_uout.setZero();
            UVLM::BiotSavart::vortex_ring(target_triad,
                                          zeta[0].template block<2,2>(i,j),
                                          zeta[1].template block<2,2>(i,j),
                                          zeta[2].template block<2,2>(i,j),
                                          gamma(i,j),
                                          temp_uout);
            uout[0](i, j) = temp_uout(0);
            uout[1](i, j) = temp_uout(1);
            uout[2](i, j) = temp_uout(2);
            if (horseshoe)
            {
                if (i == Mend - 1)
                {
                    temp_uout.setZero();
                    UVLM::BiotSavart::horseshoe(target_triad,
                                                zeta_star[0].template block<2,2>(ii,j),
                                                zeta_star[1].template block<2,2>(ii,j),
                                                zeta_star[2].template block<2,2>(ii,j),
                                                gamma_star(ii,j),
                                                temp_uout);
                    uout[0](i, j) += temp_uout(0);
                    uout[1](i, j) += temp_uout(1);
                    uout[2](i, j) += temp_uout(2);
                }
            }
        }
    }
    if (!horseshoe)
    {
        const uint mstar = gamma_star.rows();
        const uint i = Mend - 1;
        for (uint j=Nstart; j<Nend; ++j)
        {
            for (uint i_star=0; i_star<mstar; ++i_star)
            {
                temp_uout.setZero();
                UVLM::BiotSavart::vortex_ring(target_triad,
                                              zeta_star[0].template block<2,2>(i_star, j),
                                              zeta_star[1].template block<2,2>(i_star, j),
                                              zeta_star[2].template block<2,2>(i_star, j),
                                              gamma_star(i_star, j),
                                              temp_uout);
                uout[0](i, j) += temp_uout(0);
                uout[1](i, j) += temp_uout(1);
                uout[2](i, j) += temp_uout(2);
            }
        }
    }
}

template <typename t_zeta,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_ttriad,
          typename t_uout>
void UVLM::BiotSavart::surface_with_unsteady_wake
(
    const t_zeta&       zeta,
    const t_zeta_star&  zeta_star,
    const t_gamma&      gamma,
    const t_gamma_star& gamma_star,
    const t_ttriad&     target_triad,
    t_uout&             uout,
    const bool&         image_method,
    const int&          n_rows // default val = -1
)
{
    const uint Mstart = 0;
    const uint Nstart = 0;
    const uint Mend = gamma.rows();
    const uint Nend = gamma.cols();

    UVLM::Types::Vector3 temp_uout;
    const uint ii = 0;
    // surface contribution
    for (unsigned int i=Mstart; i<Mend; ++i)
    {
        for (unsigned int j=Nstart; j<Nend; ++j)
        {
            temp_uout.setZero();
            UVLM::BiotSavart::vortex_ring(target_triad,
                                          zeta[0].template block<2,2>(i,j),
                                          zeta[1].template block<2,2>(i,j),
                                          zeta[2].template block<2,2>(i,j),
                                          gamma(i,j),
                                          temp_uout);
            uout[0](i, j) = temp_uout(0);
            uout[1](i, j) = temp_uout(1);
            uout[2](i, j) = temp_uout(2);
        }
    }
    // wake contribution
    // n_rows controls the number of panels that are included
    // in the final result. Usually for unsteady wake, the value
    // will be 1 when computing AIC coeffs.
    const uint mstar = (n_rows == -1) ? gamma_star.rows():n_rows;
    const uint i = Mend - 1;
    for (uint j=Nstart; j<Nend; ++j)
    {
        for (uint i_star=0; i_star<mstar; ++i_star)
        {
            temp_uout.setZero();
            UVLM::BiotSavart::vortex_ring(target_triad,
                                          zeta_star[0].template block<2,2>(i_star, j),
                                          zeta_star[1].template block<2,2>(i_star, j),
                                          zeta_star[2].template block<2,2>(i_star, j),
                                          gamma_star(i_star, j),
                                          temp_uout);
            uout[0](i, j) += temp_uout(0);
            uout[1](i, j) += temp_uout(1);
            uout[2](i, j) += temp_uout(2);
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
    const bool&         image_method,
    const t_normals&    normal,
    const UVLM::Types::Real vortex_radius
)
{
    const unsigned int rows_collocation = target_surface[0].rows();
    const unsigned int cols_collocation = target_surface[0].cols();

    UVLM::Types::VecMatrixX temp_uout;
    UVLM::Types::allocate_VecMat(temp_uout, zeta, -1);

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
                        uout(collocation_counter, surface_counter) +=
                            temp_uout[0](i_surf, j_surf)*normal[0](i_col, j_col) +
                            temp_uout[1](i_surf, j_surf)*normal[1](i_col, j_col) +
                            temp_uout[2](i_surf, j_surf)*normal[2](i_col, j_col);
                    }
                }
        }
    }
}





template <typename t_zeta,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_tsurface,
          typename t_uout,
          typename t_normals>
void UVLM::BiotSavart::multisurface_steady_wake
(
    const t_zeta&       zeta,
    const t_zeta_star&  zeta_star,
    const t_gamma&      gamma,
    const t_gamma_star& gamma_star,
    const t_tsurface&   target_surface,
    const bool&         horseshoe,
    t_uout&             uout,
    const bool&         image_method,
    const t_normals&    normal,
    const UVLM::Types::Real vortex_radius
)
{
    const unsigned int rows_collocation = target_surface[0].rows();
    const unsigned int cols_collocation = target_surface[0].cols();

    UVLM::Types::VecMatrixX temp_uout;
    UVLM::Types::allocate_VecMat(temp_uout, zeta, -1);

    int collocation_counter = -1;
    int surface_counter;
    for (unsigned int i_col=0; i_col<rows_collocation; ++i_col)
    {
        for (unsigned int j_col=0; j_col<cols_collocation; ++j_col)
        {
            ++collocation_counter;
            UVLM::Types::Vector3 target_triad;
            target_triad << target_surface[0](i_col, j_col),
                            target_surface[1](i_col, j_col),
                            target_surface[2](i_col, j_col);
            UVLM::Types::initialise_VecMat(temp_uout, 0.0);

            UVLM::BiotSavart::surface_with_steady_wake(zeta,
                                                       zeta_star,
                                                       gamma,
                                                       gamma_star,
                                                       target_triad,
                                                       horseshoe,
                                                       temp_uout
                                                      );

            unsigned int surf_rows = gamma.rows();
            unsigned int surf_cols = gamma.cols();

            surface_counter = -1;
            for (unsigned int i_surf=0; i_surf<surf_rows; ++i_surf)
            {
                for (unsigned int j_surf=0; j_surf<surf_cols; ++j_surf)
                {
                    ++surface_counter;
                    uout(collocation_counter, surface_counter) +=
                        temp_uout[0](i_surf, j_surf)*normal[0](i_col, j_col) +
                        temp_uout[1](i_surf, j_surf)*normal[1](i_col, j_col) +
                        temp_uout[2](i_surf, j_surf)*normal[2](i_col, j_col);
                }
            }
        }
    }
}

template <typename t_zeta,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_tsurface,
          typename t_uout,
          typename t_normals>
void UVLM::BiotSavart::multisurface_unsteady_wake
(
    const t_zeta&       zeta,
    const t_zeta_star&  zeta_star,
    const t_gamma&      gamma,
    const t_gamma_star& gamma_star,
    const t_tsurface&   target_surface,
    t_uout&             uout,
    const bool&         image_method,
    const t_normals&    normal,
    const int&          n_rows // default val = -1
)
{
    const unsigned int rows_collocation = target_surface[0].rows();
    const unsigned int cols_collocation = target_surface[0].cols();

    UVLM::Types::VecMatrixX temp_uout;
    UVLM::Types::allocate_VecMat(temp_uout, zeta, -1);

    int collocation_counter = -1;
    int surface_counter;
    for (unsigned int i_col=0; i_col<rows_collocation; ++i_col)
    {
        for (unsigned int j_col=0; j_col<cols_collocation; ++j_col)
        {
            ++collocation_counter;
            UVLM::Types::Vector3 target_triad;
            target_triad << target_surface[0](i_col, j_col),
                            target_surface[1](i_col, j_col),
                            target_surface[2](i_col, j_col);
            UVLM::Types::initialise_VecMat(temp_uout, 0.0);

            UVLM::BiotSavart::surface_with_unsteady_wake(zeta,
                                                         zeta_star,
                                                         gamma,
                                                         gamma_star,
                                                         target_triad,
                                                         temp_uout,
                                                         image_method,
                                                         n_rows
                                                        );

            unsigned int surf_rows = gamma.rows();
            unsigned int surf_cols = gamma.cols();

            surface_counter = -1;
            for (unsigned int i_surf=0; i_surf<surf_rows; ++i_surf)
            {
                for (unsigned int j_surf=0; j_surf<surf_cols; ++j_surf)
                {
                    ++surface_counter;
                    uout(collocation_counter, surface_counter) +=
                        temp_uout[0](i_surf, j_surf)*normal[0](i_col, j_col) +
                        temp_uout[1](i_surf, j_surf)*normal[1](i_col, j_col) +
                        temp_uout[2](i_surf, j_surf)*normal[2](i_col, j_col);
                }
            }
        }
    }
}

// template <typename t_zeta,
//           typename t_gamma,
//           typename t_zeta_col,
//           typename t_u_ind>
// void UVLM::BiotSavart::multisurface_on_multisurface
// (
//     const t_zeta& zeta,
//     const t_gamma& gamma,
//     const t_zeta_col& zeta_col,
//     const bool image_method,
//     t_u_ind& u_ind
// )
// {
//
// }

template <typename t_zeta,
          typename t_gamma,
          typename t_zeta_col,
          typename t_u_ind>
void UVLM::BiotSavart::whole_surface_on_surface
(
    const t_zeta& zeta,
    const t_gamma& gamma,
    const t_zeta_col& zeta_col,
    t_u_ind& u_ind,
    const bool image_method
)
{
    const uint col_n_M = zeta_col[0].rows();
    const uint col_n_N = zeta_col[0].cols();
    const uint n_M = zeta[0].rows();
    const uint n_N = zeta[0].cols();
    UVLM::Types::Vector3 target_triad;
    UVLM::Types::Vector3 uout;
    for (uint col_i_M=0; col_i_M<col_n_M; ++col_i_M)
    {
        for (uint col_j_N=0; col_j_N<col_n_N; ++col_j_N)
        {
            target_triad << zeta_col[0](col_i_M, col_j_N),
                            zeta_col[1](col_i_M, col_j_N),
                            zeta_col[2](col_i_M, col_j_N);
            uout.setZero();
            UVLM::BiotSavart::whole_surface
            (
                zeta,
                gamma,
                target_triad,
                uout,
                image_method
            );
            u_ind[0](col_i_M, col_j_N) += uout(0);
            u_ind[1](col_i_M, col_j_N) += uout(1);
            u_ind[2](col_i_M, col_j_N) += uout(2);
        }
    }
}


template <typename t_zeta,
          typename t_gamma,
          typename t_ttriad,
          typename t_uout>
void UVLM::BiotSavart::whole_surface
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

    uout.setZero();
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
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_uout>
void UVLM::BiotSavart::total_induced_velocity_on_wake
(
    const t_zeta&       zeta,
    const t_zeta_star&  zeta_star,
    const t_gamma&      gamma,
    const t_gamma_star& gamma_star,
    t_uout&             uout,
    const bool&         image_method,
    const UVLM::Types::Real vortex_radius
)
{
    const uint n_surf = zeta.size();
    for (uint col_i_surf=0; col_i_surf<n_surf; ++col_i_surf)
    {
        for (uint i_surf=0; i_surf<n_surf; ++i_surf)
        {
            // wake on wake
            UVLM::BiotSavart::whole_surface_on_surface
            (
                zeta_star[i_surf],
                gamma_star[i_surf],
                zeta_star[col_i_surf],
                uout[col_i_surf],
                image_method
            );
            // surface on wake
            UVLM::BiotSavart::whole_surface_on_surface
            (
                zeta[i_surf],
                gamma[i_surf],
                zeta_star[col_i_surf],
                uout[col_i_surf],
                image_method
            );
        }
    }
}
