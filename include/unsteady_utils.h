#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "sources.h"


/**
 * @brief Namespace for the UVLM library.
 */
namespace UVLM
{
    /**
     * @brief Namespace including unsteady simulation related functions.
     */
    namespace Unsteady
    {    
        /**
         * @brief Namespace including utility functions for the unsteady namespace.
         */
        namespace Utils
        {
            /**
             * @brief Computes the resultant grid velocity considering rigid body motion, structural displacements, and flow velocities.
             *
             * @tparam t_zeta Type of the zeta (surface grid points) data structure.
             * @tparam t_zeta_dot Type of the zeta_dot (collocation point) data structure.
             * @tparam t_uext Type of the uext (external velocities) data structure.
             * @tparam t_rbm_velocity Type of the rbm_velocity (rigid body motion velocity) data structure.
             * @tparam t_centre_rot Type of the centre_rot (center of rotation) data structure.
             * @tparam t_uext_out Type of the uext_out (output external velocities) data structure.
             *
             * @param zeta Corner point coordinates of the vortex ring panels located on the surfaces.
             * @param zeta_dot Collocation point coordinates of the panels located on the surfaces.
             * @param uext External velocities.
             * @param rbm_velocity Rigid body motion velocity.
             * @param centre_rot Center of rotation of the body.
             * @param uext_out Output resultant grid velocity considering rigid body motion.
             */
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
            /**
             * @brief Computes the resultant grid velocity considering rigid body motion, translation, and solid velocity.
             *
             * @tparam t_zeta Type of the zeta (collocation points) data structure.
             * @tparam t_zeta_dot Type of the zeta_dot (collocation point velocities) data structure.
             * @tparam t_uext Type of the uext (external velocities) data structure.
             * @tparam t_rbm_velocity Type of the rbm_velocity (rigid body motion velocity) data structure.
             * @tparam t_uext_out Type of the uext_out (output external velocities) data structure.
             * @tparam t_solid_vel Type of the solid_vel (solid velocity) data structure.
             *
             * @param zeta Collocation points of the lifting surface.
             * @param zeta_dot Velocities of collocation points.
             * @param uext External velocities.
             * @param rbm_velocity Rigid body motion velocity.
             * @param uext_out Output resultant grid velocity considering rigid body motion and solid velocity.
             * @param solid_vel Solid velocity.
            */
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

            template <typename t_zeta_star,
                      typename t_gamma_star,
                      typename t_extra_gamma_star,
                      typename t_extra_zeta_star>
            void store_last_wake_panel_information
            (
                t_zeta_star& zeta_star,
                t_gamma_star& gamma_star,
                t_extra_gamma_star& extra_gamma_star,
                t_extra_zeta_star& extra_zeta_star,
                const uint n_surf
            );

            template <typename t_struct_lifting_surfaces>
            void free_wake_convection_lifting
            (
                UVLM::Types::VecVecMatrixX& u_convection,
                t_struct_lifting_surfaces& lifting_surfaces,
                const UVLM::Types::UVMopts& options
            );
            template <typename t_struct_lifting_surfaces>
            void free_wake_final_convection
            (  
                const UVLM::Types::UVMopts& options,
                UVLM::Types::VecVecMatrixX& u_convection,
                t_struct_lifting_surfaces& lifting_surfaces
            );


            template <typename t_struct_lifting_surfaces>
            void convect_unsteady_wake
            (
                const UVLM::Types::UVMopts& options,
                t_struct_lifting_surfaces& lifting_surfaces
            );
            template <typename t_struct_lifting_surfaces,
                    typename t_struct_phantom_surf,
                    typename t_struct_nl_body>
            void convect_unsteady_wake
            (
                const UVLM::Types::UVMopts& options,
                t_struct_lifting_surfaces& lifting_surfaces,
                t_struct_phantom_surf& phantom_surfaces,
                t_struct_nl_body& nl_body
            );

            template <typename t_zeta_star,
                      typename t_struct_nl_body>
            void induced_velocity_from_sources_on_wake
            (
                t_zeta_star& zeta_star,
                t_struct_nl_body& nl_body,
                UVLM::Types::VecVecMatrixX& u_induced_out
            );
            
            template <typename t_sigma_flat,
                        typename t_aic,
                        typename t_u_ind>
            void calculate_induced_velocity_col
            (
                const t_sigma_flat& sigma_flat,
                const t_aic& u_induced_x,
                const t_aic& u_induced_y,
                const t_aic& u_induced_z,
                t_u_ind& u_induced_col
            );

            template <typename t_rbm_velocity,
                        typename t_uext_star,
                        typename t_zeta_star,
                        typename t_uext_star_total>
            void get_uext_star_total
            (
                const t_rbm_velocity& rbm_velocity,
                const t_uext_star& uext_star,
                const t_zeta_star& zeta_star,
                t_uext_star_total& uext_star_total
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
    // New vector needed here
    UVLM::Types::Vector6 vec_rbm_vel_g;
    vec_rbm_vel_g << rbm_velocity[0], rbm_velocity[1], rbm_velocity[2], rbm_velocity[3], rbm_velocity[4], rbm_velocity[5];
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
                zeta_temp << zeta[i_surf][0](i_row, i_col) - centre_rot[0],
                             zeta[i_surf][1](i_row, i_col) - centre_rot[1],
                             zeta[i_surf][2](i_row, i_col) - centre_rot[2];
                    w_cross_zeta =
                        vec_rbm_vel_g.template block<3,1> (3, 0).cross(zeta_temp);
                    for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                    {
                        uext_out[i_surf][i_dim](i_row, i_col) =
                                                uext[i_surf][i_dim](i_row, i_col)
                                                - w_cross_zeta(i_dim)
                                                - zeta_dot[i_surf][i_dim](i_row, i_col)
                                                - vec_rbm_vel_g(i_dim);
                    }
                }
        }
        
    }
    // }
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
    UVLM::Types::Vector6 vec_rbm_vel_g;
    vec_rbm_vel_g << rbm_velocity[0], rbm_velocity[1], rbm_velocity[2], rbm_velocity[3], rbm_velocity[4], rbm_velocity[5];
    
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
                        vec_rbm_vel_g.template block<3,1> (3, 0).cross(zeta_temp);
                    for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                    {
                        solid_vel[i_surf][i_dim](i_row, i_col) =
                                                w_cross_zeta(i_dim)
                                                + zeta_dot[i_surf][i_dim](i_row, i_col)
                                                + vec_rbm_vel_g(i_dim);
                        uext_out[i_surf][i_dim](i_row, i_col) =
                                                uext[i_surf][i_dim](i_row, i_col)
                                                - solid_vel[i_surf][i_dim](i_row, i_col);
                    }
            }
        }
    }
    // }
}


// wake convection
// the UVMopts flag convection_scheme determines how the
// wake is convected.
// convection_scheme == 0 => prescribed and fixed
// convection_scheme == 1 => prescribed following deformations of the wing
// convection_scheme == 2 => free, convection based on u_ext
// convection_scheme == 3 => free, convection based on u_ext and induced velocities.
template <typename t_struct_lifting_surfaces>
void UVLM::Unsteady::Utils::convect_unsteady_wake
(
    const UVLM::Types::UVMopts& options,
    t_struct_lifting_surfaces& lifting_surfaces
)
{
    const uint n_surf = options.NumSurfaces;

    if (options.convection_scheme == 0)
    {

        UVLM::Unsteady::Utils::store_last_wake_panel_information(lifting_surfaces.zeta_star, lifting_surfaces.gamma_star, lifting_surfaces.extra_gamma_star, lifting_surfaces.extra_zeta_star, n_surf);

        UVLM::Wake::General::displace_VecMat(lifting_surfaces.gamma_star);
        UVLM::Types::copy_VecVecMat(lifting_surfaces.uext_star, lifting_surfaces.uext_star_total);
    } else if (options.convection_scheme == 1)
    {
        std::cerr << "convection_scheme == "
                  << options.convection_scheme
                  << " is not yet implemented in the UVLM solver"
                  << std::endl;
    } else if (options.convection_scheme == 2)
    {
        UVLM::Types::VecVecMatrixX zeros;
        UVLM::Types::allocate_VecVecMat(zeros, lifting_surfaces.uext_star);
        // total stream velocity
        
        UVLM::Types::Vector6 rbm_no_omega = UVLM::Types::Vector6::Zero();

        rbm_no_omega.template head<3>() = lifting_surfaces.rbm_vel_g.template head<3>();
        UVLM::Types::Vector3 centre_rot = UVLM::Types::Vector3::Zero();

        UVLM::Unsteady::Utils::compute_resultant_grid_velocity
        (
            lifting_surfaces.zeta_star,
            zeros,
            lifting_surfaces.uext_star,
            rbm_no_omega,
            centre_rot,
            lifting_surfaces.uext_star_total
        );
        // convection with uext + delta u (perturbation)
        // (no u_induced)
        UVLM::Wake::Discretised::convect(lifting_surfaces.zeta_star,
                                         lifting_surfaces.uext_star_total,
                                         options.dt);

        UVLM::Unsteady::Utils::store_last_wake_panel_information(lifting_surfaces.zeta_star, lifting_surfaces.gamma_star, lifting_surfaces.extra_gamma_star, lifting_surfaces.extra_zeta_star, n_surf);

        // displace both zeta and gamma
        UVLM::Wake::General::displace_VecMat(lifting_surfaces.gamma_star);
        UVLM::Wake::General::displace_VecVecMat(lifting_surfaces.zeta_star);

        // copy last row of zeta into zeta_star
        UVLM::Wake::Discretised::generate_new_row
        (
            lifting_surfaces.zeta_star,
            lifting_surfaces.zeta
        );
    } else if (options.convection_scheme == 3)
    {
        // convection with uext + delta u (perturbation) + u_ind
        UVLM::Types::VecVecMatrixX u_convection;
        UVLM::Types::allocate_VecVecMat
        (
            u_convection,
            lifting_surfaces.uext_star
        );
        UVLM::Unsteady::Utils::free_wake_convection_lifting
        (
            u_convection,
            lifting_surfaces,
            options
        );
        UVLM::Unsteady::Utils::free_wake_final_convection
        (
            options,    
            u_convection,
            lifting_surfaces
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


template <typename t_struct_lifting_surfaces>
void  UVLM::Unsteady::Utils::free_wake_final_convection
(  
    const UVLM::Types::UVMopts& options,
    UVLM::Types::VecVecMatrixX& u_convection,
    t_struct_lifting_surfaces& lifting_surfaces
)
{
        // remove first row of convection velocities
        for (uint i_surf=0; i_surf<lifting_surfaces.n_surf; ++i_surf)
        {
            for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
            {
                u_convection[i_surf][i_dim].template topRows<1>().setZero();
            }
        }

        // u_convection = u_convection + uext_star;
        for (uint i_surf=0; i_surf<lifting_surfaces.n_surf; ++i_surf)
        {
            for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
            {
                u_convection[i_surf][i_dim] = u_convection[i_surf][i_dim] +
                                              lifting_surfaces.uext_star_total[i_surf][i_dim];
            }
        }

        UVLM::Wake::Discretised::convect(lifting_surfaces.zeta_star,
                                         u_convection,
                                         options.dt);
                                         
        UVLM::Unsteady::Utils::store_last_wake_panel_information(lifting_surfaces.zeta_star, lifting_surfaces.gamma_star, lifting_surfaces.extra_gamma_star, lifting_surfaces.extra_zeta_star, lifting_surfaces.n_surf);


        // displace both zeta and gamma
        UVLM::Wake::General::displace_VecMat(lifting_surfaces.gamma_star);
        UVLM::Wake::General::displace_VecVecMat(lifting_surfaces.zeta_star);

        // copy last row of zeta into zeta_star
        UVLM::Wake::Discretised::generate_new_row
        (
            lifting_surfaces.zeta_star,
            lifting_surfaces.zeta
        );
}


template <typename t_struct_lifting_surfaces>
void UVLM::Unsteady::Utils::free_wake_convection_lifting
(
    UVLM::Types::VecVecMatrixX& u_convection,
    t_struct_lifting_surfaces& lifting_surfaces,
    const UVLM::Types::UVMopts& options
)
{
    UVLM::Types::allocate_VecVecMat(lifting_surfaces.uext_star_total, lifting_surfaces.uext_star);
    UVLM::Types::VecVecMatrixX zeros;
    UVLM::Types::allocate_VecVecMat(zeros, lifting_surfaces.uext_star);
    // total stream velocity
    UVLM::Types::Vector6 vec_rbm_vel_g;
    vec_rbm_vel_g << lifting_surfaces.rbm_vel_g[0], lifting_surfaces.rbm_vel_g[1], lifting_surfaces.rbm_vel_g[2], lifting_surfaces.rbm_vel_g[3], lifting_surfaces.rbm_vel_g[4], lifting_surfaces.rbm_vel_g[5];

    UVLM::Types::Vector6 rbm_no_omega = UVLM::Types::Vector6::Zero();
    rbm_no_omega.template head<3>() = vec_rbm_vel_g.template head<3>();
    UVLM::Types::Vector3 centre_rot = UVLM::Types::Vector3::Zero();

    UVLM::Unsteady::Utils::compute_resultant_grid_velocity
    (
        lifting_surfaces.zeta_star,
        zeros,
        lifting_surfaces.uext_star,
        rbm_no_omega,
        centre_rot,
        lifting_surfaces.uext_star_total
    );

    // induced velocity by vortex rings
    UVLM::BiotSavart::total_induced_velocity_on_wake
    (
        lifting_surfaces.zeta,
        lifting_surfaces.zeta_star,
        lifting_surfaces.gamma,
        lifting_surfaces.gamma_star,
        u_convection,
        options.ImageMethod,
        options.vortex_radius_wake_ind
    );
}


template <typename t_zeta_star,
          typename t_gamma_star,
          typename t_extra_gamma_star,
          typename t_extra_zeta_star>
void UVLM::Unsteady::Utils::store_last_wake_panel_information
(
    t_zeta_star& zeta_star,
    t_gamma_star& gamma_star,
    t_extra_gamma_star& extra_gamma_star,
    t_extra_zeta_star& extra_zeta_star,
    const uint n_surf
)
{
        for (uint i_surf=0; i_surf<n_surf; ++i_surf)
        {
            for (uint i_dim=0; i_dim<3; ++i_dim)
            {
                extra_zeta_star[i_surf][i_dim].template topRows<1>() = zeta_star[i_surf][i_dim].template bottomRows<1>();
            }
            extra_gamma_star[i_surf].template topRows<1>() = gamma_star[i_surf].template bottomRows<1>();
        }
}


template <typename t_struct_lifting_surfaces,
          typename t_struct_phantom_surf,
          typename t_struct_nl_body>
void UVLM::Unsteady::Utils::convect_unsteady_wake
(
    const UVLM::Types::UVMopts& options,
    t_struct_lifting_surfaces& lifting_surfaces,
    t_struct_phantom_surf& phantom_surfaces,
    t_struct_nl_body& nl_body
)
{
    const uint n_surf = options.NumSurfaces;
   
    if (options.convection_scheme != 3 and options.convection_scheme != 4)
    {
        UVLM::Unsteady::Utils::convect_unsteady_wake
        (
            options,
            lifting_surfaces
        );
    }
    else
    {
        // Free Wake: convection with uext + delta u (perturbation) + u_ind
        UVLM::Types::VecVecMatrixX u_convection;
        UVLM::Types::allocate_VecVecMat
        (
            u_convection,
            lifting_surfaces.uext_star
        );
        UVLM::Types::VecVecMatrixX u_convection_lifting;
        UVLM::Types::allocate_VecVecMat
        (
            u_convection_lifting,
            lifting_surfaces.uext_star
        );
        UVLM::Unsteady::Utils::free_wake_convection_lifting
        (
            u_convection,
            lifting_surfaces,
            options
        );

        //induced velociy by phantom vortex rings
        UVLM::Types::VecVecMatrixX u_convection_phantom;
        UVLM::Types::allocate_VecVecMat
        (
            u_convection_phantom,
            lifting_surfaces.uext_star
        );
        // update gamma phantom star
        // TO-DO: Check if update is needed here
        // phantom_surfaces.update_gamma_wake(zeta_star, gamma_star);
        UVLM::BiotSavart::total_induced_velocity_of_phantom_panels_on_wake
        (
            phantom_surfaces.zeta, //phantom
            phantom_surfaces.zeta_star,
            lifting_surfaces.zeta_star, //phantom + lifting
            phantom_surfaces.gamma, //phantom
            phantom_surfaces.gamma_star, //phantom
            u_convection_phantom,
            options.ImageMethod,
            options.vortex_radius_wake_ind
        );
        UVLM::Triads::VecVecMatrix_addition(u_convection_lifting, u_convection_phantom, u_convection);
        if (!options.phantom_wing_test)
        {
                //Calculate induced velocity by sources on wake
                UVLM::Types::VecVecMatrixX u_convection_nonlifting;
                UVLM::Unsteady::Utils::induced_velocity_from_sources_on_wake
        (
            lifting_surfaces.zeta_star,
            nl_body,
            u_convection_nonlifting
        );
        u_convection += u_convection_nonlifting;

        if (options.convection_scheme == 4)
        {
            //  For ignoring induced velocities on wake near junctions to prevent local instabilities
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                {
                    for (uint i_row = 0; i_row <  u_convection[i_surf][0].rows(); i_row++)
                        {
                            // TODO: Check if also induced phantom velocity has to be subtracted?
                            u_convection[i_surf][i_dim](i_row, 0) -= u_convection_nonlifting[i_surf][i_dim](i_row, 0);
                        }
                    }
                }
            }
        }
        
        UVLM::Unsteady::Utils::free_wake_final_convection
        (
            options,   
            u_convection,
            lifting_surfaces
        );  

        phantom_surfaces.update_wake(lifting_surfaces.zeta_star);
    }
    return;
}

template <typename t_zeta_star,
		  typename t_struct_nl_body>
void UVLM::Unsteady::Utils::induced_velocity_from_sources_on_wake
(
	t_zeta_star& zeta_star,
	t_struct_nl_body& nl_body,
    UVLM::Types::VecVecMatrixX& u_induced_out
)
{
	// Get surface vectors of zeta_star (account for points instead of panels)
	// determine convection velocity u_ind from non lifting surfaces
	UVLM::Types::VecVecMatrixX normals_star, longitudinals_star, perpendiculars_star;
	
	UVLM::Types::allocate_VecVecMat(normals_star, zeta_star);
	UVLM::Types::allocate_VecVecMat(longitudinals_star, zeta_star);
	UVLM::Types::allocate_VecVecMat(perpendiculars_star, zeta_star);
	UVLM::Geometry::generate_surface_vectors_wake(zeta_star, normals_star, longitudinals_star, perpendiculars_star);
		
	// Allocate matrices for source influence
	uint Ktotal_star = UVLM::Matrix::get_total_VecVecMat_size(normals_star);
	UVLM::Types::MatrixX u_induced_x = UVLM::Types::MatrixX::Zero(Ktotal_star, nl_body.Ktotal);
	UVLM::Types::MatrixX u_induced_y = UVLM::Types::MatrixX::Zero(Ktotal_star, nl_body.Ktotal);
	UVLM::Types::MatrixX u_induced_z = UVLM::Types::MatrixX::Zero(Ktotal_star, nl_body.Ktotal);

	// Get induced velocities by sources
	UVLM::Matrix::AIC_sources(nl_body.zeta,
								zeta_star,
								nl_body.longitudinals,
								nl_body.perpendiculars,
								nl_body.normals,
								longitudinals_star,
								perpendiculars_star,
								normals_star,
								u_induced_x,
								u_induced_y,
								u_induced_z,
                                false);
                                
	// add induced velocity by sources
	UVLM::Types::VectorX sigma_flat;
	UVLM::Matrix::reconstruct_vector_values_from_VecMatrixX(nl_body.sigma,
									sigma_flat,
									nl_body.zeta_col);

	// UVLM::Types::VecVecMatrixX u_induced_star_sources;
	UVLM::Types::allocate_VecVecMat(u_induced_out, zeta_star);
	UVLM::Unsteady::Utils::calculate_induced_velocity_col(sigma_flat,
													u_induced_x,
													u_induced_y,
													u_induced_z,
													u_induced_out);
    // convert induced velocities to global coordinate frame
    uint M, N;
    const uint n_surf = u_induced_out.size();
    UVLM::Types::Vector3 u_induced_col_global, longitudinal_panel, normal_panel, perpendicular_panel;
    for (uint i_surf=0; i_surf<n_surf; ++i_surf)
    {
        M = u_induced_out[i_surf][0].rows();
        N = u_induced_out[i_surf][0].cols();
        for (uint i_col=0; i_col<M; ++i_col)
        {
            for (uint j_col=0; j_col<N; ++j_col)
            {
                u_induced_col_global = UVLM::Types::Vector3(u_induced_out[i_surf][0](i_col, j_col),
                                                            u_induced_out[i_surf][1](i_col, j_col),
                                                            u_induced_out[i_surf][2](i_col, j_col));
				longitudinal_panel = UVLM::Types::Vector3(longitudinals_star[i_surf][0](i_col, j_col), longitudinals_star[i_surf][1](i_col, j_col), longitudinals_star[i_surf][2](i_col, j_col));
				perpendicular_panel = UVLM::Types::Vector3(perpendiculars_star[i_surf][0](i_col, j_col), perpendiculars_star[i_surf][1](i_col, j_col), perpendiculars_star[i_surf][2](i_col, j_col));
			    normal_panel = UVLM::Types::Vector3(normals_star[i_surf][0](i_col, j_col), normals_star[i_surf][1](i_col, j_col), normals_star[i_surf][2](i_col, j_col));
                UVLM::Geometry::convert_to_global_coordinate_system(u_induced_col_global,
                                longitudinal_panel,
                                perpendicular_panel,
                                normal_panel);	
                for (uint dim=0; dim<UVLM::Constants::NDIM; dim++)
                {
                    u_induced_out[i_surf][dim](i_col, j_col)= u_induced_col_global[dim];
                }

            }
        }
    }
}

template <typename t_sigma_flat,
			typename t_aic,
			typename t_u_ind>
void UVLM::Unsteady::Utils::calculate_induced_velocity_col
(
	const t_sigma_flat& sigma_flat,
	const t_aic& u_induced_x,
	const t_aic& u_induced_y,
	const t_aic& u_induced_z,
	t_u_ind& u_induced_col
)
{
	uint n_collocation_points = u_induced_x.rows();
	uint n_panels = u_induced_x.cols();
	UVLM::Types::MatrixX u_induced_col_flat = UVLM::Types::MatrixX::Zero(3,n_collocation_points);
	for (uint i_col=0; i_col<n_collocation_points; i_col++)
	{
		for (uint j_source=0; j_source<n_panels; j_source++)
		{
			u_induced_col_flat(0,i_col) += u_induced_x(i_col, j_source)* sigma_flat(j_source);
			u_induced_col_flat(1,i_col) += u_induced_y(i_col, j_source)* sigma_flat(j_source);
			u_induced_col_flat(2,i_col) += u_induced_z(i_col, j_source)* sigma_flat(j_source);
			
		}

	}

    UVLM::Matrix::reconstruct_VecVecMatrixX_values_from_MatrixX(u_induced_col_flat,
                                        u_induced_col,
                                        u_induced_col);

}


template <typename t_rbm_velocity,
          typename t_uext_star,
          typename t_zeta_star,
          typename t_uext_star_total>
void UVLM::Unsteady::Utils::get_uext_star_total
(
    const t_rbm_velocity& rbm_velocity,
    const t_uext_star& uext_star,
    const t_zeta_star& zeta_star,
    t_uext_star_total& uext_star_total
)
{
        UVLM::Types::VecVecMatrixX zeros;
        UVLM::Types::allocate_VecVecMat(zeros, uext_star);
        // total stream velocity
        
        UVLM::Types::Vector6 rbm_no_omega = UVLM::Types::Vector6::Zero();
        UVLM::Types::Vector6 vec_rbm_vel_g;
        vec_rbm_vel_g << rbm_velocity[0], rbm_velocity[1], rbm_velocity[2],
                         rbm_velocity[3], rbm_velocity[4], rbm_velocity[5];
    
        rbm_no_omega.template head<3>() = vec_rbm_vel_g.template head<3>();
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
}