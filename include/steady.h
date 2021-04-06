#pragma once

#include "types.h"
#include "constants.h"
#include "triads.h"
#include "mapping.h"
#include "geometry.h"
#include "biotsavart.h"
#include "matrix.h"
#include "phantom.h"
#include "wake.h"
#include "postproc.h"
#include "linear_solver.h"

#include "debugutils.h"

#include <iostream>

// DECLARATIONS
namespace UVLM
{
    namespace Steady
    {
        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_uext,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_forces,
                  typename t_rbm_vel_g>
        void solver
        (
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_forces& forces,
            t_rbm_vel_g& rbm_vel_g,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );
        template <typename t_zeta,
                  typename t_uext,
                  typename t_sigma,
                  typename t_forces>
        void solver_nonlifting_body
        (
            t_zeta& zeta,
            t_uext& uext,
            t_sigma& sigma,
            t_forces& forces,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );

        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_uext,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_forces,
                  typename t_flag_zeta_phantom,
                  typename t_rbm_vel_g,
                  typename t_zeta_nonlifting,
                  typename t_uext_nonlifting,
                  typename t_sigma_nonlifting,
                  typename t_forces_nonlifting>
        void solver_lifting_and_nonlifting_bodies
        (
            t_zeta& zeta,
            t_zeta_dot& zeta_dot,
            t_uext& uext,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_forces& forces,
            t_flag_zeta_phantom& flag_zeta_phantom,
            t_rbm_vel_g& rbm_vel_g,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions,
            t_zeta_nonlifting& zeta_nonlifting,
            t_uext_nonlifting& uext_nonlifting,
            t_sigma_nonlifting& sigma_nonlifting,
            t_forces_nonlifting& forces_nonlifting,
            const UVLM::Types::VMopts& options_nonlifting
        );
        template <typename t_zeta,
                  typename t_zeta_col,
                  typename t_uext_col,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_normals>
        void solve_horseshoe
        (
            t_zeta& zeta,
            t_zeta_col& zeta_col,
            t_uext_col& uext_col,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_normals& normals,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );
        template <typename t_zeta,
                  typename t_zeta_col,
                  typename t_uext_col,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_normals>
        void solve_discretised
        (
            t_zeta& zeta,
            t_zeta_col& zeta_col,
            t_uext_col& uext_col,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_normals& normals,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );
        template <typename t_zeta,
                  typename t_zeta_col,
                  typename t_uext_col,
                  typename t_sigma,
                  typename t_u_induced_col,
                  typename t_longitudinals,
                  typename t_perpendiculars,
                  typename t_normals>
        void solve_discretised_nonlifting_body
        (
            t_zeta& zeta,
            t_zeta_col zeta_col,
            t_uext_col& uext_col,
            t_sigma& sigma,
			t_u_induced_col&u_induced_col,
            t_longitudinals& longitudinals,
            t_perpendiculars& perpendiculars,
            t_normals& normals,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );


        template <typename t_zeta_lifting,
                  typename t_zeta_col_lifting,
                  typename t_zeta_star_lifting,
                  typename t_gamma_lifting,
                  typename t_gamma_star_lifting,
                  typename t_uext_col_lifting,
                  typename t_surface_vec_lifting,
                  typename t_zeta_nonlifting,
                  typename t_zeta_col_nonlifting,
                  typename t_uext_col_nonlifting,
                  typename t_u_induced_col_nonlifting,
                  typename t_sigma_nonlifting,
                  typename t_surface_vec_nonlifting,
                  typename t_zeta_phantom,
                  typename t_zeta_star_phantom,
                  typename t_surface_vec_phantom,
                  typename t_flag_zeta_phantom>
        void solve_discretised_lifting_and_nonlifting
        (
            const UVLM::Types::VMopts& options,
            const UVLM::Types::VMopts& options_nonlifting,
            const UVLM::Types::FlightConditions& flightconditions,
            const uint n_surf,
            const uint n_surf_nonlifting,
            const uint Ktotal,
            const uint Ktotal_lifting,
            const uint Ktotal_nonlifting,
            const uint Ktotal_phantom,
            t_zeta_lifting& zeta,
            t_zeta_col_lifting& zeta_col,
            t_zeta_star_lifting& zeta_star,
            t_gamma_lifting& gamma,
            t_gamma_star_lifting& gamma_star,
            t_uext_col_lifting& uext_col,
            t_surface_vec_lifting& normals,
            t_surface_vec_lifting& longitudinals,
            t_surface_vec_lifting& perpendiculars,
            t_zeta_nonlifting& zeta_nonlifting,
            t_zeta_col_nonlifting& zeta_col_nonlifting,
            t_uext_col_nonlifting& uext_col_nonlifting,
            t_u_induced_col_nonlifting& u_induced_col_nonlifting,
            t_sigma_nonlifting& sigma_nonlifting,
            t_surface_vec_nonlifting& normals_nonlifting,
            t_surface_vec_nonlifting& longitudinals_nonlifting,
            t_surface_vec_nonlifting& perpendiculars_nonlifting,
            t_zeta_phantom& zeta_phantom,
            t_zeta_star_phantom& zeta_star_phantom,
            t_surface_vec_phantom& normals_phantom,
            t_surface_vec_phantom& longitudinals_phantom,
            t_surface_vec_phantom& perpendiculars_phantom,  
            const t_flag_zeta_phantom& flag_zeta_phantom 
        );

        template <typename t_zeta,
                  typename t_zeta_col,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_uext_col,
          typename t_surface_vec>
        void wake_roll_up_lifting
        (
            t_zeta& zeta,
            t_zeta_col& zeta_col,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_uext_col& uext_total_col,
            t_surface_vec& normals,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );

        template <typename t_zeta,
                  typename t_zeta_col,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_uext_col,
                  typename t_zeta_nonlifting,
                  typename t_zeta_col_nonlifting,
                  typename t_uext_col_nonlifting,
                  typename t_u_induced_col_nonlifting,
                  typename t_sigma,
                  typename t_surface_vec_lifting,
                  typename t_surface_vec_nonlifting,
                  typename t_zeta_phantom,
                  typename t_surface_vec_phantom,
                  typename t_flag_zeta_phantom>
        void wake_roll_up_lifting_and_nonlifting
        (    
            t_zeta& zeta,
            t_zeta_col& zeta_col,
            t_zeta_star& zeta_star,
            t_gamma& gamma,
            t_gamma_star& gamma_star,
            t_uext_col& uext_total_col,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions,
            t_zeta_nonlifting& zeta_nonlifting,
            t_zeta_col_nonlifting& zeta_col_nonlifting,
            t_uext_col_nonlifting& uext_col_nonlifting,
            t_u_induced_col_nonlifting& u_induced_col_nonlifting,
            t_sigma& sigma,
            t_surface_vec_lifting& longitudinals_lifting,
            t_surface_vec_lifting& perpendiculars_lifting,
            t_surface_vec_lifting& normals_lifting,
            t_surface_vec_nonlifting& longitudinals_nonlifting,
            t_surface_vec_nonlifting& perpendiculars_nonlifting,
            t_surface_vec_nonlifting& normals_nonlifting,
            t_zeta_phantom& zeta_phantom,
            t_surface_vec_phantom& longitudinals_phantom,
            t_surface_vec_phantom& perpendiculars_phantom,
            t_surface_vec_phantom& normals_phantom,
            const UVLM::Types::VMopts& options_nonlifting,
            const uint Ktotal,
            const uint Ktotal_lifting,
            const uint Ktotal_nonlifting,
            const uint Ktotal_phantom,
            const t_flag_zeta_phantom& flag_zeta_phantom
        );

        template <typename t_zeta_star,
                  typename t_zeta_star_previous>       
        bool convergence_check_wake
        (
            const uint i_rollup,
            t_zeta_star& zeta_star,
            t_zeta_star_previous& zeta_star_previous,
            const double zeta_star_norm_first,
            const double rollup_tolerance
        );
    }
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta,
          typename t_zeta_dot,
          typename t_uext,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_forces,
          typename t_rbm_vel_g>
void UVLM::Steady::solver
(
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_forces& forces,
    t_rbm_vel_g& rbm_vel_g,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // Generate collocation points info
    //  Declaration
    UVLM::Types::VecVecMatrixX zeta_col;
    UVLM::Types::VecVecMatrixX uext_col;
    UVLM::Types::VecVecMatrixX uext_total;
    UVLM::Types::VecVecMatrixX uext_total_col;

    //  Allocation and mapping
    UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
    UVLM::Geometry::generate_colocationMesh(uext, uext_col);

    UVLM::Types::allocate_VecVecMat(uext_total, uext);
    UVLM::Types::copy_VecVecMat(uext, uext_total);
    UVLM::Types::allocate_VecVecMat(uext_total_col, uext, -1);

    // panel normals
    UVLM::Types::VecVecMatrixX normals;
    UVLM::Types::allocate_VecVecMat(normals, zeta_col);
    UVLM::Geometry::generate_surfaceNormal(zeta, normals);

    UVLM::Types::Vector3 u_steady;
    u_steady << uext[0][0](0,0),
                uext[0][1](0,0),
                uext[0][2](0,0);
    double delta_x = u_steady.norm()*options.dt;

    // total stream velocity
    UVLM::Unsteady::Utils::compute_resultant_grid_velocity
    (
        zeta,
        zeta_dot,
        uext,
        rbm_vel_g,
        uext_total
    );

    UVLM::Geometry::generate_colocationMesh(uext_total, uext_total_col);

    // if options.horseshoe, it is finished.
    if (options.horseshoe)
    {
        // solve the steady horseshoe problem
        UVLM::Steady::solve_horseshoe
        (
            zeta,
            zeta_col,
            uext_total_col,
            zeta_star,
            gamma,
            gamma_star,
            normals,
            options,
            flightconditions
        );

        UVLM::PostProc::calculate_static_forces
        (
            zeta,
            zeta_star,
            gamma,
            gamma_star,
            uext_total,
            forces,
            options,
            flightconditions
        );
        UVLM::Wake::Horseshoe::to_discretised(zeta_star,
                                              gamma_star,
                                              delta_x);
        return;
    }


    // create Wake
    // UVLM::Wake::Horseshoe::init(zeta, zeta_star, flightconditions);
    // UVLM::Wake::Horseshoe::to_discretised(zeta_star,
    //                                       gamma_star,
    //                                       delta_x);

    UVLM::Steady::solve_discretised
    (
        zeta,
        zeta_col,
        uext_total_col,
        zeta_star,
        gamma,
        gamma_star,
        normals,
        options,
        flightconditions
    );

    UVLM::Steady::wake_roll_up_lifting
    (
        zeta,
        zeta_col,
        zeta_star,
        gamma,
        gamma_star,
        uext_total_col,
        normals,
        options,
        flightconditions
    );
    UVLM::PostProc::calculate_static_forces_unsteady
    (
        zeta,
        zeta_dot,
        zeta_star,
        gamma,
        gamma_star,
        uext,
        rbm_vel_g,
        forces,
        options,
        flightconditions
    );
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta,
          typename t_uext,
          typename t_sigma,
          typename t_forces>
void UVLM::Steady::solver_nonlifting_body
(
    t_zeta& zeta,
    t_uext& uext,
    t_sigma& sigma,
    t_forces& forces,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // Generate collocation points info
	
    //  Declaration
    UVLM::Types::VecVecMatrixX zeta_col;
    UVLM::Types::VecVecMatrixX uext_col;
    UVLM::Types::VecVecMatrixX uext_total;
    UVLM::Types::VecVecMatrixX uext_total_col;
    UVLM::Types::VecVecMatrixX u_induced_col;

    //  Allocation and mapping
    UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
    UVLM::Geometry::generate_colocationMesh(uext, uext_col);
	
    UVLM::Types::allocate_VecVecMat(uext_total, uext);
    UVLM::Types::copy_VecVecMat(uext, uext_total);
    UVLM::Types::allocate_VecVecMat(uext_total_col, uext, -1);
	UVLM::Types::allocate_VecVecMat(u_induced_col, uext_col);

    // panel normals, longitudinal and perpendicular vectors
    UVLM::Types::VecVecMatrixX normals;
    UVLM::Types::allocate_VecVecMat(normals, zeta_col);
    UVLM::Types::VecVecMatrixX longitudinals;
    UVLM::Types::allocate_VecVecMat(longitudinals, zeta_col);
    UVLM::Types::VecVecMatrixX perpendiculars;
    UVLM::Types::allocate_VecVecMat(perpendiculars, zeta_col);
    UVLM::Geometry::generate_surface_vectors(zeta, normals, longitudinals, perpendiculars);

    UVLM::Steady::solve_discretised_nonlifting_body
    (
        zeta,
        zeta_col,
        uext_col,
        sigma,
        u_induced_col,
        longitudinals,
        perpendiculars,
        normals,
        options,
        flightconditions
    );

    UVLM::PostProc::calculate_static_forces_nonlifting_body
    (
        zeta,
        sigma,
        normals,
        longitudinals,
        perpendiculars,
        uext_col,
        u_induced_col,
        forces,
        options,
        flightconditions
    );
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/

template <typename t_zeta,
            typename t_zeta_dot,
            typename t_uext,
            typename t_zeta_star,
            typename t_gamma,
            typename t_gamma_star,
            typename t_forces,
            typename t_flag_zeta_phantom,
            typename t_rbm_vel_g,
            typename t_zeta_nonlifting,
            typename t_uext_nonlifting,
            typename t_sigma_nonlifting,
            typename t_forces_nonlifting>
void UVLM::Steady::solver_lifting_and_nonlifting_bodies
(
    t_zeta& zeta,
    t_zeta_dot& zeta_dot,
    t_uext& uext,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_forces& forces,
    t_flag_zeta_phantom& flag_zeta_phantom,
    t_rbm_vel_g& rbm_vel_g,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    t_zeta_nonlifting& zeta_nonlifting,
    t_uext_nonlifting& uext_nonlifting,
    t_sigma_nonlifting& sigma_nonlifting,
    t_forces_nonlifting& forces_nonlifting,
    const UVLM::Types::VMopts& options_nonlifting
)
{
    // Generate collocation points info for lifting surfaces
    //  Declaration
    UVLM::Types::VecVecMatrixX zeta_col;
    UVLM::Types::VecVecMatrixX uext_col;
    UVLM::Types::VecVecMatrixX uext_total;
    UVLM::Types::VecVecMatrixX uext_total_col;

    //  Allocation and mapping
    UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
    UVLM::Geometry::generate_colocationMesh(uext, uext_col);

    UVLM::Types::allocate_VecVecMat(uext_total, uext);
    UVLM::Types::copy_VecVecMat(uext, uext_total);
    UVLM::Types::allocate_VecVecMat(uext_total_col, uext, -1);

    // panel normals lifting surface
    UVLM::Types::VecVecMatrixX normals;
    UVLM::Types::allocate_VecVecMat(normals, zeta_col);   
    UVLM::Types::VecVecMatrixX longitudinals ;
    UVLM::Types::allocate_VecVecMat(longitudinals , zeta_col );
    UVLM::Types::VecVecMatrixX perpendiculars ;
    UVLM::Types::allocate_VecVecMat(perpendiculars , zeta_col ); 
    UVLM::Geometry::generate_surface_vectors(zeta, normals, longitudinals, perpendiculars);

    UVLM::Types::Vector3 u_steady;
    u_steady << uext[0][0](0,0),
                uext[0][1](0,0),
                uext[0][2](0,0);
    double delta_x = u_steady.norm()*options.dt;

    // total stream velocity
    UVLM::Unsteady::Utils::compute_resultant_grid_velocity
    (
        zeta,
        zeta_dot,
        uext,
        rbm_vel_g,
        uext_total
    );

    UVLM::Geometry::generate_colocationMesh(uext_total, uext_total_col);

    // -------- Phantom Panels -------
      
    // Phantom Cells required if true in flag zeta phantom array
    const bool phantom_cell_required = UVLM::Phantom::check_for_true_in_bool_vec_mat(flag_zeta_phantom);
    
    // Generate Panels, Collocation points and Surface Vec Phantom
    UVLM::Types::VecVecMatrixX zeta_phantom, zeta_star_phantom, zeta_col_phantom,normals_phantom, longitudinals_phantom, perpendiculars_phantom;
    std::cout<<"\n---------CHECK FOR PHANTOM PANELS!---------\n";
    if (phantom_cell_required)
    {
        UVLM::Phantom::create_phantom_zeta(zeta, zeta_phantom, flag_zeta_phantom);
        UVLM::Phantom::create_phantom_zeta_star(zeta_phantom, zeta_star, zeta_star_phantom); 
        UVLM::Geometry::generate_colocationMesh(zeta_phantom, zeta_col_phantom);
        UVLM::Types::allocate_VecVecMat(normals_phantom, zeta_phantom, -1);   
        UVLM::Types::allocate_VecVecMat(longitudinals_phantom , zeta_phantom,-1 );
        UVLM::Types::allocate_VecVecMat(perpendiculars_phantom, zeta_phantom, -1); 
        UVLM::Geometry::generate_surface_vectors(zeta_phantom, normals_phantom, longitudinals_phantom, perpendiculars_phantom);
    }
    // Generate collocation points info and panel vectors for nonlifting bodies
    //  Declaration
    UVLM::Types::VecVecMatrixX zeta_col_nonlifting;
    UVLM::Types::VecVecMatrixX uext_col_nonlifting;
    UVLM::Types::VecVecMatrixX uext_total_nonlifting;
    UVLM::Types::VecVecMatrixX uext_total_col_nonlifting;
    UVLM::Types::VecVecMatrixX u_induced_col_nonlifting;
    
    //  Allocation and mapping
    UVLM::Geometry::generate_colocationMesh(zeta_nonlifting, zeta_col_nonlifting);
    UVLM::Geometry::generate_colocationMesh(uext_nonlifting, uext_col_nonlifting);

    UVLM::Types::allocate_VecVecMat(uext_total_nonlifting, uext_nonlifting);
    UVLM::Types::copy_VecVecMat(uext_nonlifting, uext_total_nonlifting);
    UVLM::Types::allocate_VecVecMat(uext_total_col_nonlifting, uext_nonlifting, -1);	
	UVLM::Types::allocate_VecVecMat(u_induced_col_nonlifting, uext_col_nonlifting);

    // panel normals, longitudinal and perpendicular vectors
    UVLM::Types::VecVecMatrixX normals_nonlifting;
    UVLM::Types::allocate_VecVecMat(normals_nonlifting, zeta_col_nonlifting);
    UVLM::Types::VecVecMatrixX longitudinals_nonlifting;
    UVLM::Types::allocate_VecVecMat(longitudinals_nonlifting, zeta_col_nonlifting);
    UVLM::Types::VecVecMatrixX perpendiculars_nonlifting;
    UVLM::Types::allocate_VecVecMat(perpendiculars_nonlifting, zeta_col_nonlifting);
    UVLM::Geometry::generate_surface_vectors(zeta_nonlifting, normals_nonlifting, longitudinals_nonlifting, perpendiculars_nonlifting);
   
    // ########################################
    const uint n_surf = options.NumSurfaces;
    const uint n_surf_nonlifting = options_nonlifting.NumSurfaces;
    // size of rhs lifting surface 
    uint Ktotal_lifting = UVLM::Matrix::get_total_VecVecMat_size(uext_col);
    uint Ktotal_nonlifting = UVLM::Matrix::get_total_VecVecMat_size(uext_col_nonlifting);
    uint Ktotal_phantom = UVLM::Matrix::get_total_VecVecMat_size(normals_phantom);
    uint Ktotal = Ktotal_lifting + Ktotal_nonlifting;
    if (options.horseshoe)
    {
    //To-Do: Check if special solver for horseshoe needed
    UVLM::Steady::solve_discretised_lifting_and_nonlifting
    (
        options,
        options_nonlifting,
        flightconditions,
        n_surf,
        n_surf_nonlifting,
        Ktotal,
        Ktotal_lifting,
        Ktotal_nonlifting,
        Ktotal_phantom,
        zeta,
        zeta_col,
        zeta_star,
        gamma,
        gamma_star,
        uext_total_col,
        normals,
        longitudinals,
        perpendiculars,
        zeta_nonlifting,
        zeta_col_nonlifting,
        uext_col_nonlifting,
        u_induced_col_nonlifting,
        sigma_nonlifting,
        normals_nonlifting,
        longitudinals_nonlifting,
        perpendiculars_nonlifting,
        zeta_phantom,
        zeta_star_phantom,
        normals_phantom,
        longitudinals_phantom,
        perpendiculars_phantom,
        flag_zeta_phantom
        );
        if (options.Steady) 
        {
            int in_n_rows = -1;
            UVLM::Wake::Horseshoe::circulation_transfer(gamma,
                                                        gamma_star,
                                                        in_n_rows);
        }

        UVLM::Wake::Horseshoe::to_discretised(zeta_star,
                                              gamma_star,
                                              delta_x);
    }
    else
    {
        std::cout << "\n -------------- Vortex ring problems! -------------- \n";
        UVLM::Steady::solve_discretised_lifting_and_nonlifting
        (
        options,
        options_nonlifting,
        flightconditions,
        n_surf,
        n_surf_nonlifting,
        Ktotal,
        Ktotal_lifting,
        Ktotal_nonlifting,
        Ktotal_phantom,
        zeta,
        zeta_col,
        zeta_star,
        gamma,
        gamma_star,
        uext_total_col,
        normals,
        longitudinals,
        perpendiculars,
        zeta_nonlifting,
        zeta_col_nonlifting,
        uext_col_nonlifting,
        u_induced_col_nonlifting,
        sigma_nonlifting,
        normals_nonlifting,
        longitudinals_nonlifting,
        perpendiculars_nonlifting,
        zeta_phantom,
        zeta_star_phantom,
        normals_phantom,
        longitudinals_phantom,
        perpendiculars_phantom,
        flag_zeta_phantom
        );
        
        //To Do Use struct
        /*UVLM::Steady::wake_roll_up_lifting_and_nonlifting
        (
            zeta,
            zeta_col,
            zeta_star,
            gamma,
            gamma_star,
            uext_total_col,
            options,
            flightconditions,
            zeta_nonlifting,
            zeta_col_nonlifting,
            uext_col_nonlifting,
            u_induced_col_nonlifting,
            sigma_nonlifting,
            longitudinals,
            perpendiculars,
            normals,
            longitudinals_nonlifting,
            perpendiculars_nonlifting,
            normals_nonlifting,
            zeta_phantom,
            longitudinals_phantom,
            perpendiculars_phantom,
            normals_phantom,
            options_nonlifting,
            Ktotal,
            Ktotal_lifting,
            Ktotal_nonlifting,
            Ktotal_phantom,
            flag_zeta_phantom
        );*/
    }
    std::cout << "\n Calculate static forces nonlifting!\n";
    UVLM::PostProc::calculate_static_forces_nonlifting_body
    (
    zeta_nonlifting,
    sigma_nonlifting,
    normals_nonlifting,
    longitudinals_nonlifting,
    perpendiculars_nonlifting,
    uext_col_nonlifting,
    u_induced_col_nonlifting,
    forces_nonlifting,
        options_nonlifting,
        flightconditions
    );

    UVLM::PostProc::calculate_static_forces_unsteady
    (
        zeta,
    zeta_dot,
    zeta_star,
    gamma,
    gamma_star,
    uext,
    rbm_vel_g,
    forces,
    options,
            flightconditions
    );
}
/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta,
          typename t_zeta_col,
          typename t_uext_col,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_normals>
void UVLM::Steady::solve_horseshoe
(
    t_zeta& zeta,
    t_zeta_col& zeta_col,
    t_uext_col& uext_col,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_normals& normals,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // wake generation for horseshoe initialisation
    UVLM::Wake::Horseshoe::init(zeta,
                                zeta_star,
                                flightconditions);

    const uint n_surf = options.NumSurfaces;
    // size of rhs
    uint Ktotal = UVLM::Matrix::get_total_VecVecMat_size(uext_col);
    // RHS generation
    UVLM::Types::VectorX rhs;
    UVLM::Matrix::RHS(zeta_col,
                      zeta_star,
                      uext_col,
                      gamma_star,
                      normals,
                      options,
                      rhs,
                      Ktotal);

    // AIC generation
    UVLM::Types::MatrixX aic = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
    UVLM::Matrix::AIC(zeta,
                      zeta_col,
                      zeta_star,
                      uext_col,
                      normals,
                      options,
                      true,
                      aic);

    UVLM::Types::VectorX gamma_flat;
    UVLM::Matrix::deconstruct_gamma(gamma,
                                    gamma_flat,
                                    zeta_col);

    UVLM::LinearSolver::solve_system
    (
        aic,
        rhs,
        options,
        gamma_flat
    );

    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(gamma_flat,
                                    gamma,
                                    zeta_col);

    // copy gamma from trailing edge to wake if steady solution
    UVLM::Wake::Horseshoe::circulation_transfer(gamma,
                                                gamma_star);

}



/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta,
          typename t_zeta_col,
          typename t_uext_col,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_normals>
void UVLM::Steady::solve_discretised
(
    t_zeta& zeta,
    t_zeta_col& zeta_col,
    t_uext_col& uext_col,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_normals& normals,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    const uint n_surf = options.NumSurfaces;
    // size of rhs
    uint Ktotal = UVLM::Matrix::get_total_VecVecMat_size(uext_col);

    UVLM::Types::VectorX rhs;
    UVLM::Types::MatrixX aic = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
    // RHS generation
    UVLM::Matrix::RHS(zeta_col,
                      zeta_star,
                      uext_col,
                      gamma_star,
                      normals,
                      options,
                      rhs,
                      Ktotal);

    // AIC generation
    UVLM::Matrix::AIC(zeta,
                      zeta_col,
                      zeta_star,
                      uext_col,
                      normals,
                      options,
                      false,
                      aic);

    // linear system solution
    UVLM::Types::VectorX gamma_flat;
    UVLM::Matrix::deconstruct_gamma(gamma,
                                    gamma_flat,
                                    zeta_col);


    UVLM::LinearSolver::solve_system
    (
        aic,
        rhs,
        options,
        gamma_flat
    );

    // gamma flat to gamma
    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(gamma_flat,
                                    gamma,
                                    zeta_col);

    // copy gamma from trailing edge to wake
    int in_n_rows = -1;
    // if (!options.Steady) {in_n_rows = 1;}
    // UVLM::Wake::Horseshoe::circulation_transfer(gamma,
    //                                             gamma_star,
    //                                             in_n_rows);
    if (options.Steady) {
        UVLM::Wake::Horseshoe::circulation_transfer(gamma,
                                                gamma_star,
                                                in_n_rows);
    }
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_zeta,
          typename t_zeta_col,
          typename t_uext_col,
          typename t_sigma,
          typename t_u_induced_col,
          typename t_longitudinals,
          typename t_perpendiculars,
          typename t_normals>
void UVLM::Steady::solve_discretised_nonlifting_body
(
    t_zeta& zeta,
    t_zeta_col zeta_col,
    t_uext_col& uext_col,
    t_sigma& sigma,
    t_u_induced_col&u_induced_col,
    t_longitudinals& longitudinals,
    t_perpendiculars& perpendiculars,
    t_normals& normals,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    const uint n_surf = options.NumSurfaces;
    // size of rhs
    uint Ktotal = UVLM::Matrix::get_total_VecVecMat_size(uext_col);
    UVLM::Types::VectorX rhs;
    UVLM::Types::MatrixX aic_sources_x = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
    UVLM::Types::MatrixX aic_sources_y = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
    UVLM::Types::MatrixX aic_sources_z = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
    // RHS generation
    UVLM::Matrix::RHS_nonlifting_body(uext_col,
                              normals,
                              rhs,
                              Ktotal,
                              n_surf);
    // AIC generation
    UVLM::Matrix::AIC_sources(zeta,
                              zeta_col,
                              longitudinals,
                              perpendiculars,
                              normals,
                              longitudinals,
                              perpendiculars,
                              normals,
                              options,
							  aic_sources_x,
							  aic_sources_y,
							  aic_sources_z);

    // linear system solution
    UVLM::Types::VectorX sigma_flat;
    UVLM::Matrix::deconstruct_gamma(sigma,
                                    sigma_flat,
                                    zeta_col);

    UVLM::LinearSolver::solve_system
    (
        aic_sources_z,
        rhs,
        options,
        sigma_flat
    );

	UVLM::PostProc::calculate_induced_velocity_col(sigma_flat,
												   aic_sources_x,
												   aic_sources_y,
												   aic_sources_z,
												   u_induced_col);

    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(sigma_flat,
                                    sigma,
                                    uext_col);

}

template <typename t_zeta_lifting,
          typename t_zeta_col_lifting,
          typename t_zeta_star_lifting,
          typename t_gamma_lifting,
          typename t_gamma_star_lifting,
          typename t_uext_col_lifting,
          typename t_surface_vec_lifting,
          typename t_zeta_nonlifting,
          typename t_zeta_col_nonlifting,
          typename t_uext_col_nonlifting,
          typename t_u_induced_col_nonlifting,
          typename t_sigma_nonlifting,
          typename t_surface_vec_nonlifting,
          typename t_zeta_phantom,
          typename t_zeta_star_phantom,
          typename t_surface_vec_phantom,
          typename t_flag_zeta_phantom>
void UVLM::Steady::solve_discretised_lifting_and_nonlifting
(
    const UVLM::Types::VMopts& options,
    const UVLM::Types::VMopts& options_nonlifting,
    const UVLM::Types::FlightConditions& flightconditions,
    const uint n_surf,
    const uint n_surf_nonlifting,
    const uint Ktotal,
    const uint Ktotal_lifting,
    const uint Ktotal_nonlifting,
    const uint Ktotal_phantom,
    t_zeta_lifting& zeta,
    t_zeta_col_lifting& zeta_col,
    t_zeta_star_lifting& zeta_star,
    t_gamma_lifting& gamma,
    t_gamma_star_lifting& gamma_star,
    t_uext_col_lifting& uext_col,
    t_surface_vec_lifting& normals,
    t_surface_vec_lifting& longitudinals,
    t_surface_vec_lifting& perpendiculars,
    t_zeta_nonlifting& zeta_nonlifting,
    t_zeta_col_nonlifting& zeta_col_nonlifting,
    t_uext_col_nonlifting& uext_col_nonlifting,
    t_u_induced_col_nonlifting& u_induced_col_nonlifting,
    t_sigma_nonlifting& sigma_nonlifting,
    t_surface_vec_nonlifting& normals_nonlifting,
    t_surface_vec_nonlifting& longitudinals_nonlifting,
    t_surface_vec_nonlifting& perpendiculars_nonlifting,
    t_zeta_phantom& zeta_phantom,
    t_zeta_star_phantom& zeta_star_phantom,
    t_surface_vec_phantom& normals_phantom,
    t_surface_vec_phantom& longitudinals_phantom,
    t_surface_vec_phantom& perpendiculars_phantom,  
    const t_flag_zeta_phantom& flag_zeta_phantom  
)
{

    // RHS generation
    std::cout << "\n ----------- RHS generation ----------- \n";
    std::cout << " \n --generate RHS lifting -- \n";
    UVLM::Types::VectorX rhs_lifting;
    UVLM::Types::VectorX rhs_nonlifting;
    UVLM::Matrix::RHS(zeta_col,
                      zeta_star,
                      uext_col,
                      gamma_star,
                      normals,
                      options,
                      rhs_lifting,
                      Ktotal_lifting);

    UVLM::Matrix::RHS_nonlifting_body(uext_col_nonlifting,
                                      normals_nonlifting,
                                      rhs_nonlifting,
                                      Ktotal_nonlifting,
                                      n_surf_nonlifting);

    UVLM::Types::VectorX rhs_phantom, rhs;
    rhs_phantom.setZero(Ktotal_phantom);

        UVLM::Types::VectorX rhs_lifting_and_nonlifting = UVLM::Types::join_vectors(rhs_lifting, rhs_nonlifting);
        rhs = UVLM::Types::join_vectors(rhs_lifting_and_nonlifting, rhs_phantom);

    // AIC generation
    std::cout << "\n ----------- AIC generation ----------- \n";    
    UVLM::Types::MatrixX aic_lifting = UVLM::Types::MatrixX::Zero(Ktotal_lifting, Ktotal_lifting);
    UVLM::Types::MatrixX aic_nonlifting_x = UVLM::Types::MatrixX::Zero(Ktotal_nonlifting, Ktotal_nonlifting);
    UVLM::Types::MatrixX aic_nonlifting_y = UVLM::Types::MatrixX::Zero(Ktotal_nonlifting, Ktotal_nonlifting);
    UVLM::Types::MatrixX aic_nonlifting_z = UVLM::Types::MatrixX::Zero(Ktotal_nonlifting, Ktotal_nonlifting);
    
    UVLM::Types::MatrixX aic = UVLM::Types::MatrixX::Zero(Ktotal+Ktotal_phantom, Ktotal+Ktotal_phantom);
    UVLM::Matrix::aic_combined(Ktotal,
                                Ktotal_lifting,
                                Ktotal_nonlifting, 
                                Ktotal_phantom, 
                                options,
                                options_nonlifting,
                                zeta,
                                zeta_col,
                                zeta_star,
                                uext_col,
                                normals,
                                longitudinals,
                                perpendiculars,
                                aic_lifting,
                                zeta_nonlifting,
                                zeta_col_nonlifting,
                                uext_col_nonlifting,
                                normals_nonlifting,
                                longitudinals_nonlifting,
                                perpendiculars_nonlifting, 
                                zeta_phantom,
                                zeta_star_phantom,
                                normals_phantom,
                                longitudinals_phantom,
                                perpendiculars_phantom,      
                                flag_zeta_phantom,      
                                aic_nonlifting_x,
                                aic_nonlifting_y,
                                aic_nonlifting_z,
                                aic);


    // linear system solution
    UVLM::Types::VectorX gamma_and_sigma_flat = UVLM::Types::VectorX::Zero(Ktotal+Ktotal_phantom);
    UVLM::LinearSolver::solve_system
    (
        aic,
        rhs,
        options,
        gamma_and_sigma_flat
    );
    // split gamma and sigma Vector
    UVLM::Types::VectorX gamma_flat = gamma_and_sigma_flat.head(Ktotal_lifting);
    UVLM::Types::VectorX gamma_phantom_flat = gamma_and_sigma_flat.tail(Ktotal_phantom);    
    UVLM::Types::VectorX sigma_flat = gamma_and_sigma_flat.block(Ktotal_lifting, 0, Ktotal_nonlifting, 1);  

    // gamma flat to gamma
    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(gamma_flat,
                                    gamma,
                                    zeta_col);
    // copy gamma from trailing edge to wake
    int in_n_rows = -1;
    if (options.Steady) {
        UVLM::Wake::Horseshoe::circulation_transfer(gamma,
                                                gamma_star,
                                                in_n_rows);
    }

    // ---- Post-processing nonlifting body ----
	UVLM::PostProc::calculate_induced_velocity_col(sigma_flat,
												   aic_nonlifting_x,
												   aic_nonlifting_y,
												   aic_nonlifting_z,
												   u_induced_col_nonlifting);
   
    UVLM::Matrix::reconstruct_gamma(sigma_flat,
                                    sigma_nonlifting,
                                    zeta_col_nonlifting);
}

template <typename t_zeta,
          typename t_zeta_col,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_uext_col,
          typename t_surface_vec>
void UVLM::Steady::wake_roll_up_lifting
(
    t_zeta& zeta,
    t_zeta_col& zeta_col,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_uext_col& uext_total_col,
    t_surface_vec& normals,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    double zeta_star_norm_first = 0.0;
    double zeta_star_norm_previous = 0.0;
    double zeta_star_norm = 0.0;
    unsigned int N;
    bool convergence;

    UVLM::Types::VecVecMatrixX zeta_star_previous;
    if (options.n_rollup != 0)
    {
        zeta_star_norm_first = UVLM::Types::norm_VecVec_mat(zeta_star);
        zeta_star_norm_previous = zeta_star_norm_first;
        zeta_star_norm = 0.0;
        UVLM::Types::allocate_VecVecMat(zeta_star_previous, zeta_star);
        UVLM::Types::copy_VecVecMat(zeta_star, zeta_star_previous);
    }

    // ROLLUP LOOP--------------------------------------------------------
    for (uint i_rollup=0; i_rollup<options.n_rollup; ++i_rollup)
    {
        // determine convection velocity u_ind
        UVLM::Types::VecVecMatrixX u_ind;
        UVLM::Types::allocate_VecVecMat(u_ind,
                                        zeta_star);
        // induced velocity by vortex rings
        UVLM::BiotSavart::total_induced_velocity_on_wake(
            zeta,
            zeta_star,
            gamma,
            gamma_star,
            u_ind,
            options.ImageMethod,
            options.vortex_radius_wake_ind);

        // Do not move the vertices in the TE
        UVLM::Wake::Discretised::no_movement_of_vortices_at_TE(u_ind);

        // convect based on u_ind for all the grid.
        UVLM::Wake::Discretised::convect(zeta_star,
                                         u_ind,
                                         options.dt);

        // generate AIC again
        if (i_rollup%options.rollup_aic_refresh == 0)
        {
            UVLM::Steady::solve_discretised
            (
                zeta,
                zeta_col,
                uext_total_col,
                zeta_star,
                gamma,
                gamma_star,
                normals,
                options,
                flightconditions
            );
            
         }
        // convergence check -------------------
        zeta_star_norm = UVLM::Types::norm_VecVec_mat(zeta_star);
        convergence = convergence_check_wake(i_rollup,
                                             zeta_star,
                                             zeta_star_previous,
                                             zeta_star_norm_first,
                                             options.rollup_tolerance);
        if (convergence){break;}
        zeta_star_norm_previous = zeta_star_norm;
        UVLM::Types::copy_VecVecMat(zeta_star, zeta_star_previous);
    }
}
template <typename t_zeta,
          typename t_zeta_col,
          typename t_zeta_star,
          typename t_gamma,
          typename t_gamma_star,
          typename t_uext_col,
          typename t_zeta_nonlifting,
          typename t_zeta_col_nonlifting,
          typename t_uext_col_nonlifting,
          typename t_u_induced_col_nonlifting,
          typename t_sigma,
          typename t_surface_vec_lifting,
          typename t_surface_vec_nonlifting,
          typename t_zeta_phantom,
          typename t_surface_vec_phantom,
          typename t_flag_zeta_phantom>
void UVLM::Steady::wake_roll_up_lifting_and_nonlifting
(
    t_zeta& zeta,
    t_zeta_col& zeta_col,
    t_zeta_star& zeta_star,
    t_gamma& gamma,
    t_gamma_star& gamma_star,
    t_uext_col& uext_total_col,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    t_zeta_nonlifting& zeta_nonlifting,
    t_zeta_col_nonlifting& zeta_col_nonlifting,
    t_uext_col_nonlifting& uext_col_nonlifting,
    t_u_induced_col_nonlifting& u_induced_col_nonlifting,
    t_sigma& sigma,
    t_surface_vec_lifting& longitudinals_lifting,
    t_surface_vec_lifting& perpendiculars_lifting,
    t_surface_vec_lifting& normals_lifting,
    t_surface_vec_nonlifting& longitudinals_nonlifting,
    t_surface_vec_nonlifting& perpendiculars_nonlifting,
    t_surface_vec_nonlifting& normals_nonlifting,
    t_zeta_phantom& zeta_phantom,
    t_surface_vec_phantom& longitudinals_phantom,
    t_surface_vec_phantom& perpendiculars_phantom,
    t_surface_vec_phantom& normals_phantom,
    const UVLM::Types::VMopts& options_nonlifting,
    const uint Ktotal,
    const uint Ktotal_lifting,
    const uint Ktotal_nonlifting,
    const uint Ktotal_phantom,
    const t_flag_zeta_phantom& flag_zeta_phantom
)
{
    double zeta_star_norm_first = 0.0;
    double zeta_star_norm_previous = 0.0;
    double zeta_star_norm = 0.0;
    unsigned int N;
    bool convergence;

    UVLM::Types::VecVecMatrixX zeta_star_previous;
    if (options.n_rollup != 0)
    {
        zeta_star_norm_first = UVLM::Types::norm_VecVec_mat(zeta_star);
        zeta_star_norm_previous = zeta_star_norm_first;
        zeta_star_norm = 0.0;

        UVLM::Types::allocate_VecVecMat(zeta_star_previous, zeta_star);
        UVLM::Types::copy_VecVecMat(zeta_star, zeta_star_previous);
    }
    

    // ROLLUP LOOP--------------------------------------------------------
    for (uint i_rollup=0; i_rollup<options.n_rollup; ++i_rollup)
    {
        // determine convection velocity u_ind
        UVLM::Types::VecVecMatrixX u_ind;
        UVLM::Types::allocate_VecVecMat(u_ind,
                                        zeta_star);
        // induced velocity by vortex rings
        UVLM::BiotSavart::total_induced_velocity_on_wake(
            zeta,
            zeta_star,
            gamma,
            gamma_star,
            u_ind,
            options.ImageMethod,
            options.vortex_radius_wake_ind);

        // Get surface vectors of zeta_star (account for points instead of panels)
        UVLM::Types::VecVecMatrixX normals_star;
        UVLM::Types::allocate_VecVecMat(normals_star, zeta_star);   
        UVLM::Types::VecVecMatrixX longitudinals_star;
        UVLM::Types::allocate_VecVecMat(longitudinals_star, zeta_star);
        UVLM::Types::VecVecMatrixX perpendiculars_star;
        UVLM::Types::allocate_VecVecMat(perpendiculars_star , zeta_star); 
        UVLM::Geometry::generate_surface_vectors_wake(zeta_star, normals_star, longitudinals_star, perpendiculars_star);
            
            // Allocate matrices for source influence
        uint Ktotal_star = UVLM::Matrix::get_total_VecVecMat_size(normals_star);
        UVLM::Types::MatrixX u_induced_x = UVLM::Types::MatrixX::Zero(Ktotal_star, Ktotal_nonlifting);
        UVLM::Types::MatrixX u_induced_y = UVLM::Types::MatrixX::Zero(Ktotal_star, Ktotal_nonlifting);
        UVLM::Types::MatrixX u_induced_z = UVLM::Types::MatrixX::Zero(Ktotal_star, Ktotal_nonlifting);

        // Get induced velocities by sources
        UVLM::Matrix::AIC_sources(zeta_nonlifting,
                                    zeta_star,
                                    longitudinals_nonlifting,
                                    perpendiculars_nonlifting,
                                    normals_nonlifting,
                                    longitudinals_star,
                                    perpendiculars_star,
                                    normals_star,
                                    options_nonlifting,
                                    u_induced_x,
                                    u_induced_y,
                                    u_induced_z);
        // add induced velocity by sources
        UVLM::Types::VectorX sigma_flat;
        UVLM::Matrix::deconstruct_gamma(sigma,
                                        sigma_flat,
                                        zeta_col_nonlifting);

        UVLM::Types::VecVecMatrixX u_induced_star_sources;
        UVLM::Types::allocate_VecVecMat(u_induced_star_sources, u_ind);
	    UVLM::PostProc::calculate_induced_velocity_col(sigma_flat,
												       u_induced_x,
												       u_induced_y,
												       u_induced_z,
												       u_induced_star_sources);
        u_ind += u_induced_star_sources;

        // Do not move the vertices in the TE
        UVLM::Wake::Discretised::no_movement_of_vortices_at_TE(u_ind);


        // convect based on u_ind for all the grid.
        UVLM::Wake::Discretised::convect(zeta_star,
                                         u_ind,
                                         options.dt);

        // convergence check -------------------
        zeta_star_norm = UVLM::Types::norm_VecVec_mat(zeta_star);
        convergence = convergence_check_wake(i_rollup,
                                             zeta_star,
                                             zeta_star_previous,
                                             zeta_star_norm_first,
                                             options.rollup_tolerance);
        if (convergence){break;}
        zeta_star_norm_previous = zeta_star_norm;
        UVLM::Types::copy_VecVecMat(zeta_star, zeta_star_previous);
    }
}

template <typename t_zeta_star,
          typename t_zeta_star_previous>       
bool UVLM::Steady::convergence_check_wake
(
    const uint i_rollup,
    t_zeta_star& zeta_star,
    t_zeta_star_previous& zeta_star_previous,
    const double zeta_star_norm_first,
    const double rollup_tolerance
)
{
    if (i_rollup != 0)
    {
        double eps = std::abs(UVLM::Types::norm_VecVec_mat(zeta_star - zeta_star_previous))/zeta_star_norm_first;
        std::cout << "    UVLM: Rollup iteration: " << i_rollup << ". Error: " << eps << std::endl;
        if (eps < rollup_tolerance)
        {
            return true; //convergence check passed --> break
        }
    }
    return false;
    

}