#pragma once

#include "EigenInclude.h"
#include "types.h"


namespace UVLM
{
    namespace Unsteady
    {
        namespace Utils
        {
            template <typename t_zeta,
                      typename t_zeta_dot,
                      typename t_uext,
                      typename t_rbm_velocity,
                      typename t_centre_rot,
                      typename t_uext_out>
            void compute_resultant_grid_velocity
            (
                t_zeta& zeta,
                t_zeta_dot& zeta_dot,
                t_uext& uext,
                t_rbm_velocity& rbm_velocity,
                t_centre_rot& centre_rot,
                t_uext_out& uext_out
            );
            template <typename t_zeta,
                      typename t_zeta_dot,
                      typename t_uext,
                      typename t_rbm_velocity,
                      typename t_uext_out,
                      typename t_solid_vel>
            void compute_resultant_grid_velocity_solid_vel
            (
                t_zeta& zeta,
                t_zeta_dot& zeta_dot,
                t_uext& uext,
                t_rbm_velocity& rbm_velocity,
                t_uext_out& uext_out,
                t_solid_vel& solid_vel
            );
            template <typename t_zeta,
                      typename t_zeta_star,
                      typename t_gamma,
                      typename t_gamma_star,
                      typename t_uext,
                      typename t_uext_star,
                      typename t_uext_star_total,
                      typename t_rbm_velocity,
                      typename t_extra_gamma_star,
                      typename t_extra_zeta_star>
            void convect_unsteady_wake
            (
                const UVLM::Types::UVMopts& options,
                const t_zeta& zeta,
                t_zeta_star& zeta_star,
                const t_gamma& gamma,
                t_gamma_star& gamma_star,
                const t_uext& uext,
                const t_uext_star& uext_star,
                t_uext_star_total& uext_star_total,
                const t_rbm_velocity& rbm_velocity,
                t_extra_gamma_star& extra_gamma_star,
                t_extra_zeta_star& extra_zeta_star
            );
        }
    }
}

template <typename t_zeta,
          typename t_zeta_dot,
          typename t_uext,
          typename t_rbm_velocity,
          typename t_centre_rot,
          typename t_uext_out>
void UVLM::Unsteady::Utils::compute_resultant_grid_velocity
(
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_rbm_velocity& rbm_velocity,
    t_centre_rot& centre_rot,
    t_uext_out& uext_out
)
{
    const uint n_surf = zeta.size();
    UVLM::Types::Vector3 w_cross_zeta;
    UVLM::Types::Vector3 zeta_temp;
    UVLM::Types::initialise_VecVecMat(uext_out);

    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        const uint n_col = zeta[i_surf][0].cols();
        const uint n_row = zeta[i_surf][0].rows();
        for (uint i_col=0; i_col<n_col; ++i_col)
        {
            for (uint i_row=0; i_row<n_row; ++i_row)
            {
                zeta_temp << zeta[i_surf][0](i_row, i_col) - centre_rot(0),
                             zeta[i_surf][1](i_row, i_col) - centre_rot(1),
                             zeta[i_surf][2](i_row, i_col) - centre_rot(2);
                w_cross_zeta =
                    rbm_velocity.template block<3,1> (3, 0).cross(zeta_temp);
                for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                {
                    uext_out[i_surf][i_dim](i_row, i_col) =
                                              uext[i_surf][i_dim](i_row, i_col)
                                            - w_cross_zeta(i_dim)
                                            - zeta_dot[i_surf][i_dim](i_row, i_col)
                                            - rbm_velocity(i_dim);
                }
            }
        }
    }
}

template <typename t_zeta,
          typename t_zeta_dot,
          typename t_uext,
          typename t_rbm_velocity,
          typename t_uext_out,
          typename t_solid_vel>
void UVLM::Unsteady::Utils::compute_resultant_grid_velocity_solid_vel
(
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_rbm_velocity& rbm_velocity,
    t_uext_out& uext_out,
    t_solid_vel& solid_vel
)
{
    const uint n_surf = zeta.size();
    UVLM::Types::Vector3 w_cross_zeta;
    UVLM::Types::Vector3 zeta_temp;
    UVLM::Types::initialise_VecVecMat(uext_out);
    UVLM::Types::initialise_VecVecMat(solid_vel);

    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        const uint n_col = zeta[i_surf][0].cols();
        const uint n_row = zeta[i_surf][0].rows();
        for (uint i_col=0; i_col<n_col; ++i_col)
        {
            for (uint i_row=0; i_row<n_row; ++i_row)
            {
                zeta_temp << zeta[i_surf][0](i_row, i_col),
                             zeta[i_surf][1](i_row, i_col),
                             zeta[i_surf][2](i_row, i_col);
                w_cross_zeta =
                    rbm_velocity.template block<3,1> (3, 0).cross(zeta_temp);
                for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                {
                    solid_vel[i_surf][i_dim](i_row, i_col) =
                                             w_cross_zeta(i_dim)
                                            + zeta_dot[i_surf][i_dim](i_row, i_col)
                                            + rbm_velocity(i_dim);
                    uext_out[i_surf][i_dim](i_row, i_col) =
                                              uext[i_surf][i_dim](i_row, i_col)
                                            - solid_vel[i_surf][i_dim](i_row, i_col);
                }
            }
        }
    }
}

// wake convection
// the UVMopts flag convection_scheme determines how the
// wake is convected.
// convection_scheme == 0 => prescribed and fixed
// convection_scheme == 1 => prescribed following deformations of the wing
// convection_scheme == 2 => free, convection based on u_ext
// convection_scheme == 3 => free, convection based on u_ext and induced velocities.
template <typename t_zeta,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_uext,
          typename t_uext_star,
          typename t_uext_star_total,
          typename t_rbm_velocity,
          typename t_extra_gamma_star,
          typename t_extra_zeta_star>
void UVLM::Unsteady::Utils::convect_unsteady_wake
(
    const UVLM::Types::UVMopts& options,
    const t_zeta& zeta,
    t_zeta_star& zeta_star,
    const t_gamma& gamma,
    t_gamma_star& gamma_star,
    const t_uext& uext,
    const t_uext_star& uext_star,
    t_uext_star_total& uext_star_total,
    const t_rbm_velocity& rbm_velocity,
    t_extra_gamma_star& extra_gamma_star,
    t_extra_zeta_star& extra_zeta_star
)
{
    const uint n_surf = options.NumSurfaces;

    if (options.convection_scheme == 0)
    {
        // TODO: move this to a function
        for (uint i_surf=0; i_surf<n_surf; ++i_surf)
        {
            for (uint i_dim=0; i_dim<3; ++i_dim)
            {
                extra_zeta_star[i_surf][i_dim].template topRows<1>() = zeta_star[i_surf][i_dim].template bottomRows<1>();
            }
            extra_gamma_star[i_surf].template topRows<1>() = gamma_star[i_surf].template bottomRows<1>();
        }

        UVLM::Wake::General::displace_VecMat(gamma_star);
        UVLM::Types::copy_VecVecMat(uext_star, uext_star_total);
    } else if (options.convection_scheme == 1)
    {
        std::cerr << "convection_scheme == "
                  << options.convection_scheme
                  << " is not yet implemented in the UVLM solver"
                  << std::endl;
    } else if (options.convection_scheme == 2)
    {
        UVLM::Types::VecVecMatrixX zeros;
        UVLM::Types::allocate_VecVecMat(zeros, uext_star);
        // total stream velocity
        UVLM::Types::Vector6 rbm_no_omega = UVLM::Types::Vector6::Zero();
        rbm_no_omega.template head<3>() = rbm_velocity.template head<3>();
        UVLM::Types::Vector3 centre_rot = UVLM::Types::Vector3::Zero();

        UVLM::Unsteady::Utils::compute_resultant_grid_velocity
        (
            zeta_star,
            zeros,
            uext_star,
            rbm_no_omega,
            centre_rot,
            uext_star_total
        );
        // convection with uext + delta u (perturbation)
        // (no u_induced)
        UVLM::Wake::Discretised::convect(zeta_star,
                                         uext_star_total,
                                         options.dt);

        // TODO: move this to a function
        for (uint i_surf=0; i_surf<n_surf; ++i_surf)
        {
            for (uint i_dim=0; i_dim<3; ++i_dim)
            {
                extra_zeta_star[i_surf][i_dim].template topRows<1>() = zeta_star[i_surf][i_dim].template bottomRows<1>();
            }
            extra_gamma_star[i_surf].template topRows<1>() = gamma_star[i_surf].template bottomRows<1>();
        }

        // displace both zeta and gamma
        UVLM::Wake::General::displace_VecMat(gamma_star);
        UVLM::Wake::General::displace_VecVecMat(zeta_star);

        // copy last row of zeta into zeta_star
        UVLM::Wake::Discretised::generate_new_row
        (
            zeta_star,
            zeta
        );
    } else if (options.convection_scheme == 3)
    {
        // convection with uext + delta u (perturbation) + u_ind
        UVLM::Types::VecVecMatrixX u_convection;
        UVLM::Types::allocate_VecVecMat
        (
            u_convection,
            uext_star
        );
        UVLM::Types::VecVecMatrixX uext_star_total;
        UVLM::Types::allocate_VecVecMat(uext_star_total, uext_star);
        UVLM::Types::VecVecMatrixX zeros;
        UVLM::Types::allocate_VecVecMat(zeros, uext_star);
        // total stream velocity
        UVLM::Types::Vector6 rbm_no_omega = UVLM::Types::Vector6::Zero();
        rbm_no_omega.template head<3>() = rbm_velocity.template head<3>();
        UVLM::Types::Vector3 centre_rot = UVLM::Types::Vector3::Zero();

        UVLM::Unsteady::Utils::compute_resultant_grid_velocity
        (
            zeta_star,
            zeros,
            uext_star,
            rbm_no_omega,
            centre_rot,
            uext_star_total
        );
        // induced velocity by vortex rings
        UVLM::BiotSavart::total_induced_velocity_on_wake
        (
            zeta,
            zeta_star,
            gamma,
            gamma_star,
            u_convection,
            options.ImageMethod,
            options.vortex_radius_wake_ind
        );
        // remove first row of convection velocities
        for (uint i_surf=0; i_surf<n_surf; ++i_surf)
        {
            for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
            {
                u_convection[i_surf][i_dim].template topRows<1>().setZero();
            }
        }

        // u_convection = u_convection + uext_star;
        for (uint i_surf=0; i_surf<n_surf; ++i_surf)
        {
            for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
            {
                u_convection[i_surf][i_dim] = u_convection[i_surf][i_dim] +
                                              uext_star_total[i_surf][i_dim];
            }
        }

        UVLM::Wake::Discretised::convect(zeta_star,
                                         u_convection,
                                         options.dt);

        // TODO: move this to a function
        for (uint i_surf=0; i_surf<n_surf; ++i_surf)
        {
            for (uint i_dim=0; i_dim<3; ++i_dim)
            {
                extra_zeta_star[i_surf][i_dim].template topRows<1>() = zeta_star[i_surf][i_dim].template bottomRows<1>();
            }
            extra_gamma_star[i_surf].template topRows<1>() = gamma_star[i_surf].template bottomRows<1>();
        }

        // displace both zeta and gamma
        UVLM::Wake::General::displace_VecMat(gamma_star);
        UVLM::Wake::General::displace_VecVecMat(zeta_star);

        // copy last row of zeta into zeta_star
        UVLM::Wake::Discretised::generate_new_row
        (
            zeta_star,
            zeta
        );
    } else
    {
        std::cerr << "convection_scheme == "
                  << options.convection_scheme
                  << " is not supported by the UVLM solver. \n"
                  << "Supported options are from [0->3]"
                  << std::endl;
    }
    return;
}
