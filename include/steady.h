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
        template <typename t_struct_lifting_surfaces>
        void solver
        (
            t_struct_lifting_surfaces& lifting_surfaces,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );
        template <typename t_struct_nl_body>
        void solver_nonlifting_body
        (
            t_struct_nl_body& nl_body,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );

        template <typename t_struct_lifting_surface,
                  typename t_struct_phantom_surf,
                  typename t_struct_nl_body>
        void solver_lifting_and_nonlifting_bodies
        (
            t_struct_lifting_surface& lifting_surfaces,
            t_struct_phantom_surf& phantom_surfaces,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions,
            t_struct_nl_body& nl_body
        );

        template <typename t_struct_lifting_surface>
        void solve_horseshoe
        (
            t_struct_lifting_surface& lifting_surfaces,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );
        template <typename t_struct_lifting_surface>
        void solve_discretised
        (
            t_struct_lifting_surface& lifting_surfaces,
            const UVLM::Types::VMopts& options
        );
        template <typename t_struct_nl_body>
        void solve_discretised_nonlifting_body
        (
            t_struct_nl_body& nl_body,
            const UVLM::Types::VMopts& options
        );

        template <typename t_struct_lifting_surfaces,
                  typename t_struct_nl_body,
                  typename t_struct_phantom_surf>
        void solve_discretised_lifting_and_nonlifting
        (
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions,  
            t_struct_lifting_surfaces& lifting_surfaces,
            t_struct_nl_body& nl_body,
            t_struct_phantom_surf& phantom_surfaces
        );

        template <typename t_struct_lifting_surface>
        void wake_roll_up_lifting
        (
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions,           
            t_struct_lifting_surface& lifting_surfaces
        );
        template <typename t_struct_lifting_surface,
                  typename t_struct_nonlifting_surface>
        void wake_roll_up
        (
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions,            
            t_struct_lifting_surface& lifting_surfaces,
            t_struct_nonlifting_surface& nl_body
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
template <typename t_struct_lifting_surfaces>
void UVLM::Steady::solver
(
    t_struct_lifting_surfaces& lifting_surfaces,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    lifting_surfaces.get_surface_parameters();    

     // total stream velocity
     //TO-DO: Adjust for nonlifting bodies as well! and create function in combination with uext_total_col
    UVLM::Unsteady::Utils::compute_resultant_grid_velocity
    (
        lifting_surfaces.zeta,
        lifting_surfaces.zeta_dot,
        lifting_surfaces.u_ext,
        options.rbm_vel_g,
        options.centre_rot_g,
        lifting_surfaces.uext_total
    );

    // To-DO: Interpolate to collocation points
    // To-Do: Check if Interpolation is needed afterall
    UVLM::Geometry::generate_colocationMesh(lifting_surfaces.uext_total, lifting_surfaces.uext_total_col);

    // if options.horseshoe, it is finished.
    if (options.horseshoe)
    {
        // solve the steady horseshoe problem
        UVLM::Steady::solve_horseshoe
        (
            lifting_surfaces,
            options,
            flightconditions
        );

        UVLM::PostProc::calculate_static_forces
        (
            lifting_surfaces.zeta,
            lifting_surfaces.zeta_star,
            lifting_surfaces.gamma,
            lifting_surfaces.gamma_star,
            lifting_surfaces.uext_total,
            lifting_surfaces.forces,
            options,
            flightconditions
        );
        double delta_x = options.dt * UVLM::Types::norm_Vec(lifting_surfaces.u_ext[0][0](0,0),
                                                            lifting_surfaces.u_ext[0][1](0,0),
                                                            lifting_surfaces.u_ext[0][2](0,0));
        UVLM::Wake::Horseshoe::to_discretised(lifting_surfaces.zeta_star,
                                              lifting_surfaces.gamma_star,
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
        lifting_surfaces,
        options
    );
    
    // UVLM::Steady::wake_roll_up_lifting
    // (
    //     options,
    //     flightconditions,
    //     lifting_surfaces
    // );
    UVLM::PostProc::calculate_static_forces_unsteady
    (
        lifting_surfaces.zeta,
        lifting_surfaces.zeta_dot,
        lifting_surfaces.zeta_star,
        lifting_surfaces.gamma,
        lifting_surfaces.gamma_star,
        lifting_surfaces.u_ext,
        options.rbm_vel_g,
        lifting_surfaces.forces,
        options,
        flightconditions
    );
    
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_struct_nl_body>
void UVLM::Steady::solver_nonlifting_body
(
    t_struct_nl_body& nl_body,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{    
    nl_body.get_surface_parameters();
    UVLM::Steady::solve_discretised_nonlifting_body
    (
        nl_body,
        options
    );
    // export_data_to_csv_file("zeta_x.csv", nl_body.zeta[0][0]);
    // export_data_to_csv_file("zeta_x_col.csv", nl_body.zeta_col[0][0]);

    UVLM::PostProc::calculate_static_forces_nonlifting_body
    (
        nl_body,
        flightconditions
    );
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_struct_lifting_surface,
          typename t_struct_phantom_surf,
          typename t_struct_nl_body>
void UVLM::Steady::solver_lifting_and_nonlifting_bodies
(
    t_struct_lifting_surface& lifting_surfaces,
    t_struct_phantom_surf& phantom_surfaces,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    t_struct_nl_body& nl_body
)
{
    // Lifting surfaces
    lifting_surfaces.get_surface_parameters();    
     // total stream velocity
    UVLM::Unsteady::Utils::compute_resultant_grid_velocity
    (
        lifting_surfaces.zeta,
        lifting_surfaces.zeta_dot,
        lifting_surfaces.u_ext,
        options.rbm_vel_g,
        options.centre_rot_g,
        lifting_surfaces.uext_total
    );
    
    UVLM::Geometry::generate_colocationMesh(lifting_surfaces.uext_total, lifting_surfaces.uext_total_col);

    nl_body.get_surface_parameters(options.phantom_wing_test);

    // -------- Phantom Panels -------   
    phantom_surfaces.get_surface_parameters();
    phantom_surfaces.update_wake(lifting_surfaces.zeta_star);

    // ########################################
    UVLM::Steady::solve_discretised_lifting_and_nonlifting
    (
        options,
        flightconditions,
        lifting_surfaces,
        nl_body,
        phantom_surfaces
    );

    // UVLM::PostProc::calculate_static_forces_nonlifting_body
    // (
    //     nl_body,
    //     flightconditions
    // );
        UVLM::PostProc::calculate_static_forces_unsteady
        (
            lifting_surfaces,
            phantom_surfaces,
            options,
            flightconditions
        );
    }


/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_struct_lifting_surface>
void UVLM::Steady::solve_horseshoe
(
    t_struct_lifting_surface& lifting_surfaces,
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // wake generation for horseshoe initialisation
    UVLM::Wake::Horseshoe::init(lifting_surfaces.zeta,
                                lifting_surfaces.zeta_star,
                                flightconditions);
    lifting_surfaces.get_aerodynamic_solver_inputs(options);

    UVLM::Types::VectorX gamma_flat;
    UVLM::Matrix::deconstruct_gamma(lifting_surfaces.gamma,
                                    gamma_flat,
                                    lifting_surfaces.zeta_col);

    UVLM::LinearSolver::solve_system
    (
        lifting_surfaces.aic,
        lifting_surfaces.rhs,
        options,
        gamma_flat
    );

    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(gamma_flat,
                                    lifting_surfaces.gamma,
                                    lifting_surfaces.zeta_col);

    // copy gamma from trailing edge to wake if steady solution
    UVLM::Wake::Horseshoe::circulation_transfer(lifting_surfaces.gamma,
                                                lifting_surfaces.gamma_star);

}



/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_struct_lifting_surface>
void UVLM::Steady::solve_discretised
(
    t_struct_lifting_surface& lifting_surfaces,
    const UVLM::Types::VMopts& options
)
{
    lifting_surfaces.get_aerodynamic_solver_inputs(options);
    // linear system solution
    UVLM::Types::VectorX gamma_flat;
    UVLM::Matrix::deconstruct_gamma(lifting_surfaces.gamma,
                                    gamma_flat,
                                    lifting_surfaces.zeta_col);
    UVLM::LinearSolver::solve_system
    (
        lifting_surfaces.aic,
        lifting_surfaces.rhs,
        options,
        gamma_flat
    );

    // gamma flat to gamma
    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(gamma_flat,
                                    lifting_surfaces.gamma,
                                    lifting_surfaces.zeta_col);

    if (options.Steady) 
    {
        // copy gamma from trailing edge to wake
        int in_n_rows = -1;
        UVLM::Wake::Horseshoe::circulation_transfer(lifting_surfaces.gamma,
                                                lifting_surfaces.gamma_star,
                                                in_n_rows);
    }
}

/*-----------------------------------------------------------------------------

-----------------------------------------------------------------------------*/
template <typename t_struct_nl_body>
void UVLM::Steady::solve_discretised_nonlifting_body
(
    t_struct_nl_body& nl_body,
    const UVLM::Types::VMopts& options
)
{
    // Get AIC and RHS 
    nl_body.get_aerodynamic_solver_inputs();
    // linear system solution
    UVLM::Types::VectorX sigma_flat;
    UVLM::Matrix::deconstruct_gamma(nl_body.sigma,
                                    sigma_flat,
                                    nl_body.zeta_col);

    UVLM::LinearSolver::solve_system
    (
        nl_body.aic_sources_z,
        nl_body.rhs,
        options,
        sigma_flat
    );

	UVLM::PostProc::calculate_induced_velocity_col(sigma_flat,
												   nl_body.aic_sources_x,
												   nl_body.aic_sources_y,
												   nl_body.aic_sources_z,
												   nl_body.u_induced_col);
    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(sigma_flat,
                                    nl_body.sigma,
                                    nl_body.uext_col);

}
template <typename t_struct_lifting_surfaces,
            typename t_struct_nl_body,
            typename t_struct_phantom_surf>
void UVLM::Steady::solve_discretised_lifting_and_nonlifting
(
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    t_struct_lifting_surfaces& lifting_surfaces,
    t_struct_nl_body& nl_body,
    t_struct_phantom_surf& phantom_surfaces
)
{
    if (options.horseshoe)
    {
        // wake generation for horseshoe initialisation
        UVLM::Wake::Horseshoe::init(lifting_surfaces.zeta,
                                    lifting_surfaces.zeta_star,
                                    flightconditions);        
        phantom_surfaces.update_wake(lifting_surfaces.zeta_star);
    }
    // Setup Lifting and Lifting Surfaces
    lifting_surfaces.get_aerodynamic_solver_inputs(options);
    if (!options.Steady)
    {
        UVLM::Matrix::RHS_lifting_unsteady(lifting_surfaces.zeta_col,
            lifting_surfaces.zeta_star,
            lifting_surfaces.uext_col,
            lifting_surfaces.gamma_star,
            lifting_surfaces.normals,
            options,
            lifting_surfaces.rhs,
            lifting_surfaces.Ktotal,
            phantom_surfaces.gamma_star,
            phantom_surfaces.zeta_star
        );
    }
    nl_body.get_aerodynamic_solver_inputs(options.phantom_wing_test);

    // RHS generation
    // Create RHS phantom and merge all RHS vectors
    UVLM::Types::VectorX rhs_phantom, rhs;
    rhs_phantom.setZero(phantom_surfaces.Ktotal);
    if (!options.phantom_wing_test)
    {
        UVLM::Types::VectorX rhs_lifting_and_nonlifting = UVLM::Types::join_vectors(lifting_surfaces.rhs, nl_body.rhs);
        rhs = UVLM::Types::join_vectors(rhs_lifting_and_nonlifting, rhs_phantom);
    }
    else
    {
        rhs = UVLM::Types::join_vectors(lifting_surfaces.rhs, rhs_phantom);     
    }
    uint Ktotal = nl_body.Ktotal + lifting_surfaces.Ktotal + phantom_surfaces.Ktotal;
    // AIC generation
    UVLM::Types::MatrixX aic = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
    UVLM::Matrix::aic_combined(lifting_surfaces,
                                nl_body, 
                                phantom_surfaces,   
                                aic,
                                options);

    // linear system solution
    UVLM::Types::VectorX gamma_and_sigma_flat = UVLM::Types::VectorX::Zero(Ktotal);
    UVLM::LinearSolver::solve_system
    (
        aic,
        rhs,
        options,
        gamma_and_sigma_flat
    );
    // split gamma and sigma Vector
    UVLM::Types::VectorX gamma_flat = gamma_and_sigma_flat.head(lifting_surfaces.Ktotal);
    UVLM::Types::VectorX gamma_phantom_flat = gamma_and_sigma_flat.tail(phantom_surfaces.Ktotal);    
    UVLM::Types::VectorX sigma_flat = gamma_and_sigma_flat.block(lifting_surfaces.Ktotal, 0, nl_body.Ktotal, 1);  
    
    // gamma flat to gamma
    // probably could be done better with a Map
    UVLM::Matrix::reconstruct_gamma(gamma_flat,
                                    lifting_surfaces.gamma,
                                    lifting_surfaces.zeta_col);    
    UVLM::Matrix::reconstruct_gamma(gamma_phantom_flat,
                                    phantom_surfaces.gamma,
                                    phantom_surfaces.zeta_col);
    
    if (options.Steady) 
    {  
        // copy gamma from trailing edge to wake
        int in_n_rows = -1;
        UVLM::Wake::Horseshoe::circulation_transfer(lifting_surfaces.gamma,
                                                    lifting_surfaces.gamma_star,
                                                    in_n_rows);
    }

    // ---- Post-processing nonlifting body ----
    if (nl_body.Ktotal > 0)
    {
	UVLM::PostProc::calculate_induced_velocity_col(sigma_flat,
												   nl_body.aic_sources_x,
												   nl_body.aic_sources_y,
												   nl_body.aic_sources_z,
												   nl_body.u_induced_col);

    // ToDo: Integrate in function below!
	UVLM::Types::MatrixX aic_nonlifting_on_lifting_z = UVLM::Types::MatrixX::Zero(lifting_surfaces.Ktotal, nl_body.Ktotal);
    UVLM::Types::MatrixX aic_nonlifting_on_lifting_x = UVLM::Types::MatrixX::Zero(lifting_surfaces.Ktotal, nl_body.Ktotal);
    UVLM::Types::MatrixX aic_nonlifting_on_lifting_y = UVLM::Types::MatrixX::Zero(lifting_surfaces.Ktotal, nl_body.Ktotal);
    UVLM::Matrix::AIC_sources(nl_body.zeta,
                              lifting_surfaces.zeta_col,
                              nl_body.longitudinals,
                              nl_body.perpendiculars,
                              nl_body.normals,
                              lifting_surfaces.longitudinals, //collocation
                              lifting_surfaces.perpendiculars, //collocation
                              lifting_surfaces.normals, //collocation
							  aic_nonlifting_on_lifting_x,
							  aic_nonlifting_on_lifting_y,
                              aic_nonlifting_on_lifting_z,
                              false);  
    lifting_surfaces.get_induced_col_from_sources(sigma_flat, 
                                                    aic_nonlifting_on_lifting_x,
                                                    aic_nonlifting_on_lifting_y,
                                                    aic_nonlifting_on_lifting_z, 
                                                    nl_body.Ktotal);
                                    
    UVLM::Matrix::reconstruct_gamma(sigma_flat,
                                    nl_body.sigma,
                                    nl_body.zeta_col);
    }
}

template <typename t_struct_lifting_surface>
void UVLM::Steady::wake_roll_up_lifting
(
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    t_struct_lifting_surface& lifting_surfaces
)
{
    // TODO: Combine with wake_roll_up
    double zeta_star_norm_first = 0.0;
    double zeta_star_norm_previous = 0.0;
    double zeta_star_norm = 0.0;
    bool convergence;

    UVLM::Types::VecVecMatrixX zeta_star_previous;
    if (options.n_rollup != 0)
    {
        zeta_star_norm_first = UVLM::Types::norm_VecVec_mat(lifting_surfaces.zeta_star);
        zeta_star_norm_previous = zeta_star_norm_first;
        zeta_star_norm = 0.0;
        UVLM::Types::allocate_VecVecMat(zeta_star_previous, lifting_surfaces.zeta_star);
        UVLM::Types::copy_VecVecMat(lifting_surfaces.zeta_star, zeta_star_previous);
    }

    // ROLLUP LOOP--------------------------------------------------------
    for (uint i_rollup=0; i_rollup<options.n_rollup; ++i_rollup)
    {
        // determine convection velocity u_ind
        UVLM::Types::VecVecMatrixX u_ind;
        UVLM::Types::allocate_VecVecMat(u_ind,
                                        lifting_surfaces.zeta_star);
        // induced velocity by vortex rings
        UVLM::BiotSavart::total_induced_velocity_on_wake(
            lifting_surfaces.zeta,
            lifting_surfaces.zeta_star,
            lifting_surfaces.gamma,
            lifting_surfaces.gamma_star,
            u_ind,
            options.ImageMethod,
            options.vortex_radius_wake_ind);

        // Do not move the vertices in the TE
        UVLM::Wake::Discretised::no_movement_of_vortices_at_TE(u_ind);

        // convect based on u_ind for all the grid.
        UVLM::Wake::Discretised::convect(lifting_surfaces.zeta_star,
                                         u_ind,
                                         options.dt);

        // generate AIC again
        if (i_rollup%options.rollup_aic_refresh == 0)
        {
            UVLM::Steady::solve_discretised
            (
                lifting_surfaces,
                options
            );
            
         }
        // convergence check -------------------
        zeta_star_norm = UVLM::Types::norm_VecVec_mat(lifting_surfaces.zeta_star);
        convergence = convergence_check_wake(i_rollup,
                                             lifting_surfaces.zeta_star,
                                             zeta_star_previous,
                                             zeta_star_norm_first,
                                             options.rollup_tolerance);
        if (convergence){break;}
        zeta_star_norm_previous = zeta_star_norm;
        UVLM::Types::copy_VecVecMat(lifting_surfaces.zeta_star, zeta_star_previous);
    }
}
template <typename t_struct_lifting_surface,
          typename t_struct_nonlifting_surface>
void UVLM::Steady::wake_roll_up
(
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    t_struct_lifting_surface& lifting_surfaces,
    t_struct_nonlifting_surface& nl_body
)
{
    double zeta_star_norm_first = 0.0;
    double zeta_star_norm_previous = 0.0;
    double zeta_star_norm = 0.0;
    unsigned int N;
    bool convergence;
    int n_rollup = options.n_rollup;
    UVLM::Types::VecVecMatrixX zeta_star_previous;
    if (n_rollup != 0)
    {
        zeta_star_norm_first = UVLM::Types::norm_VecVec_mat(lifting_surfaces.zeta_star);
        zeta_star_norm_previous = zeta_star_norm_first;
        zeta_star_norm = 0.0;
        zeta_star_previous = UVLM::Types::allocate_VecVecMat(lifting_surfaces.zeta_star); 
        UVLM::Types::copy_VecVecMat(lifting_surfaces.zeta_star, zeta_star_previous);
    }
    
    // ROLLUP LOOP--------------------------------------------------------
    for (uint i_rollup=0; i_rollup<n_rollup; ++i_rollup)
    {

        // determine convection velocity u_ind from lifting surfaces
        UVLM::Types::VecVecMatrixX u_ind;
        UVLM::Types::allocate_VecVecMat(u_ind, lifting_surfaces.zeta_star);
        
        // induced velocity by vortex rings
        UVLM::BiotSavart::total_induced_velocity_on_wake(
            lifting_surfaces.zeta,
            lifting_surfaces.zeta_star,
            lifting_surfaces.gamma,
            lifting_surfaces.gamma_star,
            u_ind,
            options.ImageMethod,
            options.vortex_radius_wake_ind);

        if (!options.only_lifting)
        {
            // Get surface vectors of zeta_star (account for points instead of panels)
            // determine convection velocity u_ind from non lifting surfaces
            UVLM::Types::VecVecMatrixX normals_star, longitudinals_star, perpendiculars_star;
            
            UVLM::Types::allocate_VecVecMat(normals_star, lifting_surfaces.zeta_star);
            UVLM::Types::allocate_VecVecMat(longitudinals_star, lifting_surfaces.zeta_star);
            UVLM::Types::allocate_VecVecMat(perpendiculars_star, lifting_surfaces.zeta_star);
            UVLM::Geometry::generate_surface_vectors_wake(lifting_surfaces.zeta_star, normals_star, longitudinals_star, perpendiculars_star);
                
            // Allocate matrices for source influence
            uint Ktotal_star = UVLM::Matrix::get_total_VecVecMat_size(normals_star);
            UVLM::Types::MatrixX u_induced_x = UVLM::Types::MatrixX::Zero(Ktotal_star, nl_body.Ktotal);
            UVLM::Types::MatrixX u_induced_y = UVLM::Types::MatrixX::Zero(Ktotal_star, nl_body.Ktotal);
            UVLM::Types::MatrixX u_induced_z = UVLM::Types::MatrixX::Zero(Ktotal_star, nl_body.Ktotal);

            // Get induced velocities by sources
            UVLM::Matrix::AIC_sources(nl_body.zeta,
                                        lifting_surfaces.zeta_star,
                                        nl_body.longitudinals,
                                        nl_body.perpendiculars,
                                        nl_body.normals,
                                        longitudinals_star,
                                        perpendiculars_star,
                                        normals_star,
                                        u_induced_x,
                                        u_induced_y,
                                        u_induced_z);
            // add induced velocity by sources
            UVLM::Types::VectorX sigma_flat;
            UVLM::Matrix::deconstruct_gamma(nl_body.sigma,
                                            sigma_flat,
                                            nl_body.zeta_col);

            UVLM::Types::VecVecMatrixX u_induced_star_sources;
            UVLM::Types::allocate_VecVecMat(u_induced_star_sources, u_ind);
	        UVLM::PostProc::calculate_induced_velocity_col(sigma_flat,
                                                            u_induced_x,
                                                            u_induced_y,
                                                            u_induced_z,
                                                            u_induced_star_sources);
            u_ind += u_induced_star_sources;
        }

        // Do not move the vertices in the TE
        UVLM::Wake::Discretised::no_movement_of_vortices_at_TE(u_ind);


        // convect based on u_ind for all the grid.
        UVLM::Wake::Discretised::convect(lifting_surfaces.zeta_star,
                                         u_ind,
                                         options.dt);

        // convergence check -------------------
        zeta_star_norm = UVLM::Types::norm_VecVec_mat(lifting_surfaces.zeta_star);
        convergence = convergence_check_wake(i_rollup,
                                             lifting_surfaces.zeta_star,
                                             zeta_star_previous,
                                             zeta_star_norm_first,
                                             options.rollup_tolerance);
        if (convergence){break;}
        zeta_star_norm_previous = zeta_star_norm;
        UVLM::Types::copy_VecVecMat(lifting_surfaces.zeta_star, zeta_star_previous);
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
