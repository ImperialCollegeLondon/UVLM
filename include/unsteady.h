#pragma once

#include "types.h"
#include "triads.h"
#include "unsteady_utils.h"
#include "constants.h"
#include "mapping.h"
#include "geometry.h"
#include "biotsavart.h"
#include "matrix.h"
#include "postproc.h"
#include "steady.h"
#include "wake.h"

#include <iostream>

// DECLARATIONS
namespace UVLM
{
    namespace Unsteady
    {
        template <typename t_lifting_surfaces>
        void solver
        (
            const uint& i_iter,
            t_lifting_surfaces& lifting_surfaces_unsteady,
            const UVLM::Types::UVMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );        
        
        template <typename t_lifting_surfaces,
                  typename t_struct_nl_body,
                  typename t_phantom_surfaces>
        void solver_lifting_and_nonlifting
        (
            const uint& i_iter,
            t_lifting_surfaces& lifting_surfaces_unsteady,
            t_struct_nl_body& nl_body,
            t_phantom_surfaces& phantom_surfaces,
            const UVLM::Types::UVMopts& options,
            const UVLM::Types::FlightConditions& flightconditions
        );
    }
}

template <typename t_lifting_surfaces>
void UVLM::Unsteady::solver
(
    const uint& i_iter,
    t_lifting_surfaces& lifting_surfaces_unsteady,
    const UVLM::Types::UVMopts& options,
    const UVLM::Types::FlightConditions& flightconditions
)
{
    // SOLVE------------------------------------------
    const uint n_surf = options.NumSurfaces;
    const double dt = options.dt;
    // Generate collocation points info
    //  Declaration
    lifting_surfaces_unsteady.get_surface_parameters();

    UVLM::Unsteady::Utils::compute_resultant_grid_velocity_solid_vel
    (
        lifting_surfaces_unsteady.zeta,
        lifting_surfaces_unsteady.zeta_dot,
        lifting_surfaces_unsteady.u_ext,
        options.rbm_vel_g,
        lifting_surfaces_unsteady.uext_total,
        lifting_surfaces_unsteady.solid_vel
    );

    //  Allocation and mapping
    // Same in steady    
    UVLM::Geometry::generate_colocationMesh(lifting_surfaces_unsteady.uext_total,
                                            lifting_surfaces_unsteady.uext_total_col);

    // Unsteady specific
    UVLM::Types::VecVecMatrixX uext_star_total;
    UVLM::Types::allocate_VecVecMat(uext_star_total, lifting_surfaces_unsteady.uext_star);
    
    // for what is extra gamma star?
    UVLM::Types::VecMatrixX extra_gamma_star;
    UVLM::Types::VecVecMatrixX extra_zeta_star;
    extra_zeta_star.resize(n_surf);
    for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
    {
        extra_gamma_star.push_back(UVLM::Types::MatrixX());
        extra_gamma_star[i_surf].setZero(1, lifting_surfaces_unsteady.gamma_star[i_surf].cols());
        for (unsigned int i_dim=0; i_dim<3; ++i_dim)
        {
            extra_zeta_star[i_surf].push_back(UVLM::Types::MatrixX(1,
                                                                   lifting_surfaces_unsteady.gamma_star[i_surf].cols() + 1));
        }
    }
    
    UVLM::Types::VMopts steady_options = UVLM::Types::UVMopts2VMopts(options);

    // different than in steady
    if (options.convect_wake)
    {
        UVLM::Unsteady::Utils::convect_unsteady_wake
        (
            options,
            lifting_surfaces_unsteady.zeta,
            lifting_surfaces_unsteady.zeta_star,
            lifting_surfaces_unsteady.gamma,
            lifting_surfaces_unsteady.gamma_star,
            lifting_surfaces_unsteady.u_ext,
            lifting_surfaces_unsteady.uext_star,
            uext_star_total,
            extra_gamma_star,
            extra_zeta_star
        );
    }

    UVLM::Wake::Discretised::circulation_transfer(lifting_surfaces_unsteady.zeta,
                                                  lifting_surfaces_unsteady.zeta_star,
                                                  lifting_surfaces_unsteady.gamma,
                                                  lifting_surfaces_unsteady.gamma_star,
                                                  lifting_surfaces_unsteady.uext_total_col,
                                                  dt);
  
    if (!options.cfl1)
    {
        UVLM::Wake::Discretised::cfl_n1(options,
                                        lifting_surfaces_unsteady.zeta_star,
                                        lifting_surfaces_unsteady.gamma_star,
                                        extra_gamma_star,
                                        extra_zeta_star,
                                        lifting_surfaces_unsteady.dist_to_orig,
                                        uext_star_total,
                                        lifting_surfaces_unsteady.solid_vel,
                                        dt);
    }

    // we can use UVLM::Steady::solve_discretised if uext_col
    // is the total velocity including non-steady contributions.
    // TODO: Check if struct for unsteady surfaces can be transformed to steady to save mem and time?
    UVLM::Steady::solve_discretised
    (
        lifting_surfaces_unsteady,
        steady_options
    );

    if (options.quasi_steady)
    {
        UVLM::Wake::Horseshoe::circulation_transfer(lifting_surfaces_unsteady.gamma,
                                                    lifting_surfaces_unsteady.gamma_star,
                                                    -1);
    }
    
    // forces calculation
    // set forces to 0 just in case
    UVLM::Types::initialise_VecVecMat(lifting_surfaces_unsteady.forces);
    // static:
    UVLM::PostProc::calculate_static_forces_unsteady
    (
        lifting_surfaces_unsteady.zeta,
        lifting_surfaces_unsteady.zeta_dot,
        lifting_surfaces_unsteady.zeta_star,
        lifting_surfaces_unsteady.gamma,
        lifting_surfaces_unsteady.gamma_star,
        lifting_surfaces_unsteady.u_ext,
        options.rbm_vel_g,
        lifting_surfaces_unsteady.forces,
        steady_options,
        flightconditions
    );
    // TODO: Check if delete?
    UVLM::Types::initialise_VecVecMat(lifting_surfaces_unsteady.dynamic_forces);

}
    template <typename t_lifting_surfaces,
                typename t_struct_nl_body,
                typename t_phantom_surfaces>
    void UVLM::Unsteady::solver_lifting_and_nonlifting
    (
        const uint& i_iter,
        t_lifting_surfaces& lifting_surfaces_unsteady,
        t_struct_nl_body& nl_body,
        t_phantom_surfaces& phantom_surfaces,
        const UVLM::Types::UVMopts& options,
        const UVLM::Types::FlightConditions& flightconditions
    )
    {
      // SOLVE------------------------------------------
    const uint n_surf = options.NumSurfaces;
    const double dt = options.dt;
    // Generate collocation points info
    //  Declaration
    lifting_surfaces_unsteady.get_surface_parameters();

    UVLM::Unsteady::Utils::compute_resultant_grid_velocity_solid_vel
    (
        lifting_surfaces_unsteady.zeta,
        lifting_surfaces_unsteady.zeta_dot,
        lifting_surfaces_unsteady.u_ext,
        options.rbm_vel_g,
        lifting_surfaces_unsteady.uext_total,
        lifting_surfaces_unsteady.solid_vel
    );
    //  Allocation and mapping
    UVLM::Geometry::generate_colocationMesh(lifting_surfaces_unsteady.uext_total,
                                            lifting_surfaces_unsteady.uext_total_col);

    UVLM::Types::VecVecMatrixX uext_star_total;
    UVLM::Types::allocate_VecVecMat(uext_star_total, lifting_surfaces_unsteady.uext_star);

    UVLM::Types::VecMatrixX extra_gamma_star;
    UVLM::Types::VecVecMatrixX extra_zeta_star;
    extra_zeta_star.resize(n_surf);
    for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
    {
        extra_gamma_star.push_back(UVLM::Types::MatrixX());
        extra_gamma_star[i_surf].setZero(1, lifting_surfaces_unsteady.gamma_star[i_surf].cols());
        for (unsigned int i_dim=0; i_dim<3; ++i_dim)
        {
            extra_zeta_star[i_surf].push_back(UVLM::Types::MatrixX(1,
                                                                   lifting_surfaces_unsteady.gamma_star[i_surf].cols() + 1));
        }
    }
    
    if (!options.only_lifting)
    {
        // nonlifting parameters/geometry attributes        
        nl_body.get_surface_parameters();
        
        // phantom TO-DO: case if no phantom panels required? what is the input to the function below?
        phantom_surfaces.get_surface_parameters();
        phantom_surfaces.update_wake(lifting_surfaces_unsteady.zeta_star);
        phantom_surfaces.update_gamma_wake(lifting_surfaces_unsteady.zeta_star, lifting_surfaces_unsteady.gamma_star);
        phantom_surfaces.update_gamma(lifting_surfaces_unsteady.Ktotal,
                                     lifting_surfaces_unsteady.zeta_col,
                                     lifting_surfaces_unsteady.gamma);
    }
    UVLM::Types::VMopts steady_options = UVLM::Types::UVMopts2VMopts(options);

    // different than in steady
    if (options.convect_wake)
    {
        if (!options.only_lifting)
        {
            // To-Do: update gamma phantom            
            UVLM::Unsteady::Utils::convect_unsteady_wake
            (
                options,
                lifting_surfaces_unsteady.zeta,
                lifting_surfaces_unsteady.zeta_star,
                lifting_surfaces_unsteady.gamma,
                lifting_surfaces_unsteady.gamma_star,
                lifting_surfaces_unsteady.u_ext,
                lifting_surfaces_unsteady.uext_star,
                uext_star_total,
                extra_gamma_star,
                extra_zeta_star,
                phantom_surfaces,
                nl_body
            );
        }
        else
        {
            UVLM::Unsteady::Utils::convect_unsteady_wake
            (
                options,
                lifting_surfaces_unsteady.zeta,
                lifting_surfaces_unsteady.zeta_star,
                lifting_surfaces_unsteady.gamma,
                lifting_surfaces_unsteady.gamma_star,
                lifting_surfaces_unsteady.u_ext,
                lifting_surfaces_unsteady.uext_star,
                uext_star_total,
                extra_gamma_star,
                extra_zeta_star
            );
        }
        
    }
    UVLM::Wake::Discretised::circulation_transfer(lifting_surfaces_unsteady.zeta,
                                                  lifting_surfaces_unsteady.zeta_star,
                                                  lifting_surfaces_unsteady.gamma,
                                                  lifting_surfaces_unsteady.gamma_star,
                                                  lifting_surfaces_unsteady.uext_total_col,
                                                  dt);
    

    
    if (!options.cfl1)
    {
        UVLM::Wake::Discretised::cfl_n1(options,
                                        lifting_surfaces_unsteady.zeta_star,
                                        lifting_surfaces_unsteady.gamma_star,
                                        extra_gamma_star,
                                        extra_zeta_star,
                                        lifting_surfaces_unsteady.dist_to_orig,
                                        uext_star_total,
                                        lifting_surfaces_unsteady.solid_vel,
                                        dt);
    }

    // we can use UVLM::Steady::solve_discretised if uext_col
    // is the total velocity including non-steady contributions.
    if (!options.only_lifting)
    {
        // // nonlifting parameters/geometry attributes        
        // nl_body.get_surface_parameters();
        // // phantom TO-DO: case if no phantom panels required? what is the input to the function below?
        // phantom_surfaces.get_surface_parameters();
        phantom_surfaces.update_wake(lifting_surfaces_unsteady.zeta_star);
        phantom_surfaces.update_gamma_wake(lifting_surfaces_unsteady.zeta_star, lifting_surfaces_unsteady.gamma_star);
        
        UVLM::Steady::solve_discretised_lifting_and_nonlifting
        (
            steady_options,
            flightconditions,
            lifting_surfaces_unsteady,
            nl_body,
            phantom_surfaces
        );
        // UVLM::PostProc::calculate_static_forces_nonlifting_body
        // (
        //     nl_body,
        //     flightconditions
        // );
    }
    else
    {
        UVLM::Steady::solve_discretised
        (
            lifting_surfaces_unsteady,
            steady_options
        );
    }
    if (options.quasi_steady)
    {
        // TO-DO: Check if function has to be adjusted for NL Bodies
        UVLM::Wake::Horseshoe::circulation_transfer(lifting_surfaces_unsteady.gamma,
                                                    lifting_surfaces_unsteady.gamma_star,
                                                    -1);
    }

    // forces calculation
    // set forces to 0 just in case
    // TODO: Check if necessary
    UVLM::Types::initialise_VecVecMat(lifting_surfaces_unsteady.forces);
    // static:
    UVLM::Types::Vector3 centre_rot_g = UVLM::Types::Vector3::Zero();
    UVLM::PostProc::calculate_static_forces_unsteady
    (
        lifting_surfaces_unsteady,
        phantom_surfaces,
        steady_options,
        flightconditions
    );
    UVLM::Types::initialise_VecVecMat(lifting_surfaces_unsteady.dynamic_forces);
    
}