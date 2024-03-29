/**
 * @file steady.h
 * @brief This file contains the implementation of the steady-state solver for the UVLM (Unsteady Vortex Lattice Method) library.
 */

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

/**
 * @brief Namespace for the UVLM library.
 */
namespace UVLM
{
    /**
     * @brief Namespace for steady-state solver functionality within the UVLM library.
     */
    namespace Steady
    {
        /**
         * @brief Solver for steady-state aerodynamics of lifting surfaces.
         *
         * This function solves the steady-state aerodynamics problem for lifting surfaces.
         *
         * @tparam t_struct_lifting_surfaces Type of lifting surfaces data structure.
         * @param lifting_surfaces The struct containing information of the lifting surfaces to be solved.
         * @param options Options for the solver.
         * @param flightconditions Flight conditions for the solver.
         */
        template <typename t_struct_lifting_surfaces>
        void solver(
            t_struct_lifting_surfaces& lifting_surfaces,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );

        /**
         * @brief Solver for steady-state aerodynamics of non-lifting bodies.
         *
         * This function solves the steady-state aerodynamics problem for non-lifting bodies.
         *
         * @tparam t_struct_nl_body Type of non-lifting body data structure.
         * @param nl_body The struct containing information of the non-lifting body to be solved.
         * @param options Options for the solver.
         * @param flightconditions Flight conditions for the solver.
         */
        template <typename t_struct_nl_body>
        void solver_nonlifting_body(
            t_struct_nl_body& nl_body,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );

        /**
         * @brief Solver for steady-state aerodynamics of the VLM solver for lifting surfaces coupled with the 
         *        linear source panel method for non-lifting bodies.
         *
         * This function solves the steady-state aerodynamics problem for a combination of lifting and non-lifting bodies.
         *
         * @tparam t_struct_lifting_surface Type of lifting surface data structure.
         * @tparam t_struct_phantom_surf Type of phantom surface data structure.
         * @tparam t_struct_nl_body Type of non-lifting body data structure.
         * @param lifting_surfaces The lifting surface struct.
         * @param phantom_surfaces Struct of phantom surfaces used for modeling fuselage-wing junctions.
         * @param options Options for the solver.
         * @param flightconditions Flight conditions for the solver.
         * @param nl_body The struct containing information of the non-lifting body to be solved.
         */
        template <typename t_struct_lifting_surface,
                  typename t_struct_phantom_surf,
                  typename t_struct_nl_body>
        void solver_lifting_and_nonlifting_bodies(
            t_struct_lifting_surface& lifting_surfaces,
            t_struct_phantom_surf& phantom_surfaces,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions,
            t_struct_nl_body& nl_body
        );

        /**
         * @brief Solver for the steady-state VLM using horseshoe vertices for the wake model.
         *
         * This function solves the steady-state VLM with horseshoe vertices used for the wake modeling.
         *
         * @tparam t_struct_lifting_surface Type of lifting surface data structure.
         * @param lifting_surfaces The lifting surface struct.
         * @param options Options for the solver.
         * @param flightconditions Flight conditions for the solver.
         */
        template <typename t_struct_lifting_surface>
        void solve_horseshoe(
            t_struct_lifting_surface& lifting_surfaces,
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );

        /**
         * @brief  Discretised solver for the linear source panel method for lifting surfaces.
         *
         * This function sets up the linear system of equations describing the aerodynamic problem lifting surfaces.
         *
         * @tparam t_struct_lifting_surface Type of lifting surface data structure.
         * @param lifting_surfaces The struct containing information of the lifting surface.
         * @param options Options for the solver.
         */
        template <typename t_struct_lifting_surface>
        void solve_discretised(
            t_struct_lifting_surface& lifting_surfaces,
            const UVLM::Types::VMopts& options
        );

        /**
         * @brief Discretised solver for the linear source panel method for non-lifting bodies.
         * 
         * This function sets up the linear system of equations describing the aerodynamic problem for non-lifting bodies.
         * 
         * @tparam t_struct_nl_body Type of non-lifting body data structure.
         * @param nl_body The non-lifting body to be solved with the discretized model.
         * @param options Options for the solver.
         */
        template <typename t_struct_nl_body>
        void solve_discretised_nonlifting_body(
            t_struct_nl_body& nl_body,
            const UVLM::Types::VMopts& options
        );

        /**
         * @brief Discretised solver for the steady-state VLM solver for lifting surfaces coupled with the 
         *        linear source panel method for non-lifting bodies.
         *
         * This function sets up the linear system of equations describing the aerodynamic problem for a 
         * combination of lifting and non-lifting bodies.
         *
         * @tparam t_struct_lifting_surfaces Type of lifting surface data structure.
         * @tparam t_struct_nl_body Type of non-lifting body data structure.
         * @tparam t_struct_phantom_surf Type of phantom surface data structure.
         * @param options Options for the solver.
         * @param flightconditions Flight conditions for the solver.
         * @param lifting_surfaces The lifting surfaces to be solved.
         * @param nl_body The non-lifting body to be solved.
         * @param phantom_surfaces Phantom surfaces used for modeling nearby wake effects.
         */
        template <typename t_struct_lifting_surfaces,
                  typename t_struct_nl_body,
                  typename t_struct_phantom_surf>
        void solve_discretised_lifting_and_nonlifting(
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions,  
            t_struct_lifting_surfaces& lifting_surfaces,
            t_struct_nl_body& nl_body,
            t_struct_phantom_surf& phantom_surfaces
        );

        /**
         * @brief Wake roll-up solver for lifting surfaces.
         *
         * This function performs wake roll-up for lifting surfaces to capture unsteady effects.
         *
         * @tparam t_struct_lifting_surface Type of lifting surface data structure.
         * @param options Options for the solver.
         * @param flightconditions Flight conditions for the solver.
         * @param lifting_surfaces The lifting surfaces to perform wake roll-up on.
         */
        template <typename t_struct_lifting_surface>
        void wake_roll_up_lifting(
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions,           
            t_struct_lifting_surface& lifting_surfaces
        );

        /**
         * @brief Wake roll-up solver for configurations including lifting surfaces coupled with non-lifting bodies.
         *
         * This function performs wake roll-up for lifting surfaces coupled with non-lifting bodies to capture unsteady effects.
         *
         * @tparam t_struct_lifting_surface Type of lifting surface data structure.
         * @tparam t_struct_nonlifting_surface Type of non-lifting surface data structure.
         * @param options Options for the solver.
         * @param flightconditions Flight conditions for the solver.
         * @param lifting_surfaces The lifting surfaces to perform wake roll-up on.
         * @param nl_body The struct of the non-lifting body influencing the wake roll-up on.
         */
        template <typename t_struct_lifting_surface,
                  typename t_struct_nonlifting_surface>
        void wake_roll_up(
            const UVLM::Types::VMopts& options,
            const UVLM::Types::FlightConditions& flightconditions,            
            t_struct_lifting_surface& lifting_surfaces,
            t_struct_nonlifting_surface& nl_body
        );

        /**
         * @brief Convergence check for wake roll-up.
         *
         * This function checks the convergence of wake roll-up iteration.
         *
         * @tparam t_zeta_star Type of the zeta_star data structure.
         * @tparam t_zeta_star_previous Type of the previous zeta_star data structure.
         * @param i_rollup Current roll-up iteration.
         * @param zeta_star Current zeta_star.
         * @param zeta_star_previous Previous zeta_star.
         * @param zeta_star_norm_first Norm of the initial zeta_star.
         * @param rollup_tolerance Convergence tolerance for roll-up iterations.
         * @return true if convergence is achieved, false otherwise.
         */
        template <typename t_zeta_star,
                  typename t_zeta_star_previous>       
        bool convergence_check_wake(
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
        lifting_surfaces.rbm_vel_g,
        lifting_surfaces.centre_rot,
        lifting_surfaces.uext_total
    );

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
        lifting_surfaces.rbm_vel_g,
        lifting_surfaces.centre_rot,
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
        lifting_surfaces.rbm_vel_g,
        lifting_surfaces.centre_rot,
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

    UVLM::PostProc::calculate_static_forces_nonlifting_body
    (
        nl_body,
        flightconditions
    );
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
    lifting_surfaces.get_aerodynamic_solver_inputs(options);

    UVLM::Types::VectorX gamma_flat;
    UVLM::Matrix::reconstruct_vector_values_from_VecMatrixX(lifting_surfaces.gamma,
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
    UVLM::Matrix::reconstruct_VecMatrixX_values_from_vector(gamma_flat,
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
    UVLM::Matrix::reconstruct_vector_values_from_VecMatrixX(lifting_surfaces.gamma,
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
    UVLM::Matrix::reconstruct_VecMatrixX_values_from_vector(gamma_flat,
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
    UVLM::Matrix::reconstruct_vector_values_from_VecMatrixX(nl_body.sigma,
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
    UVLM::Matrix::reconstruct_VecMatrixX_values_from_vector(sigma_flat,
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
        phantom_surfaces.update_wake(lifting_surfaces.zeta_star);
    }
    // Setup Lifting and Lifting Surfaces
    lifting_surfaces.get_aerodynamic_solver_inputs(options);
    if (!options.Steady)
    {
        UVLM::Matrix::RHS_unsteady_phantom_unsteady(lifting_surfaces.zeta_col,
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

    UVLM::Matrix::reconstruct_VecMatrixX_values_from_vector(gamma_flat,
                                    lifting_surfaces.gamma,
                                    lifting_surfaces.zeta_col);    
    UVLM::Matrix::reconstruct_VecMatrixX_values_from_vector(gamma_phantom_flat,
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

    if ((!options.phantom_wing_test) && (options.consider_u_ind_by_sources_for_lifting_forces))
    {
        // nonlifting influence on lifting surface forces (or almost no influence but messing up the junction force)
        lifting_surfaces.get_coordinates_center_vertices();
        UVLM::PostProc::calculate_source_induced_velocities_on_points(
            lifting_surfaces.coordinates_center_chordwise_vertices,
            nl_body,
            sigma_flat,
            lifting_surfaces.u_induced_by_sources_on_center_chordwise_vertices
        );

        UVLM::PostProc::calculate_source_induced_velocities_on_points(
            lifting_surfaces.coordinates_center_spanwise_vertices,
            nl_body,
            sigma_flat,
            lifting_surfaces.u_induced_by_sources_on_center_spanwise_vertices
        );
    }                
    UVLM::Matrix::reconstruct_VecMatrixX_values_from_vector(sigma_flat,
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
            UVLM::Matrix::reconstruct_vector_values_from_VecMatrixX(nl_body.sigma,
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
