#pragma once

#include <cmath>
#include "EigenInclude.h"
#include "types.h"
#include "unsteady_utils.h"

namespace UVLM
{
    namespace PostProc
    {
        template <typename t_normal,
				  typename t_longitudinal,
				  typename t_perpendicular,
				  typename t_uext_col,
				  typename t_u_induced,
				  typename t_u_total_panel>
		void get_local_panel_velocity
		(
		const t_normal& normal,
		const t_longitudinal& longitudinal,
		const t_perpendicular& perpendicular,
        const t_uext_col& uext_col,
        const t_u_induced& u_induced_col,
        t_u_total_panel& u_total_panel
		);
        template <typename t_forces_collocation,
				  typename t_forces_node>
		void transform_forces_from_col_to_nodes
		(
            const t_forces_collocation& forces_collocation,
            t_forces_node& forces_node,            
            const bool rotational_body
		);
        template <typename t_forces_node,
				  typename t_forces_norm>
		void get_norm_forces
		(
            const t_forces_node& forces_node,
            t_forces_norm& forces_norm
		);
        template <typename t_forces,
                  typename t_span_seg_forces,
                  typename t_chord_seg_forces>
        void transfer_forces_to_nodes
        (
            t_forces& forces,
            t_span_seg_forces& span_seg_forces,
            t_chord_seg_forces& chord_seg_forces,
            const uint N,
            const uint M,
            const uint i_surf,
            const uint ignore_first_x_nodes_in_force_calculation
        );
        template <typename t_zeta,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_uext,
                  typename t_forces>
        void calculate_static_forces
        (
            const t_zeta& zeta,
            const t_zeta_star& zeta_star,
            const t_gamma& gamma,
            const t_gamma_star& gamma_star,
            const t_uext& uext,
            t_forces&  forces,
            const UVLM::Types::VMopts options,
            const UVLM::Types::FlightConditions& flightconditions
        )
        {
            // Set forces to 0
            UVLM::Types::initialise_VecVecMat(forces);

            // first calculate all the velocities at the corner points
            UVLM::Types::VecVecMatrixX velocities;
            UVLM::Types::allocate_VecVecMat(velocities, zeta);
            // free stream contribution
            UVLM::Types::copy_VecVecMat(uext, velocities);

            // not bothered with effciency.
            // if it is so critical, it could be improved
            const uint n_surf = zeta.size();
            UVLM::Types::Vector3 dl;
            UVLM::Types::Vector3 v;
            UVLM::Types::Vector3 f;
            UVLM::Types::Vector3 v_ind;
            UVLM::Types::Vector3 rp;
            uint start;
            uint end;
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const uint M = gamma[i_surf].rows();
                const uint N = gamma[i_surf].cols();

                for (uint i_M=0; i_M<M; ++i_M)
                {
                    for (uint i_N=0; i_N<N; ++i_N)
                    {
                        UVLM::Types::Vector3 r1;
                        UVLM::Types::Vector3 r2;
                        const unsigned int n_segment = 4;
                        for (unsigned int i_segment=0; i_segment<n_segment; ++i_segment)
                        {
                            if ((i_segment == 1) && (i_M == M - 1))
                            {
                                // trailing edge
                                continue;
                            }
                            unsigned int start = i_segment;
                            unsigned int end = (start + 1)%n_segment;
                            uint i_start = i_M + UVLM::Mapping::vortex_indices(start, 0);
                            uint j_start = i_N + UVLM::Mapping::vortex_indices(start, 1);
                            uint i_end = i_M + UVLM::Mapping::vortex_indices(end, 0);
                            uint j_end = i_N + UVLM::Mapping::vortex_indices(end, 1);


                            r1 << zeta[i_surf][0](i_start, j_start),
                                  zeta[i_surf][1](i_start, j_start),
                                  zeta[i_surf][2](i_start, j_start);
                            r2 << zeta[i_surf][0](i_end, j_end),
                                  zeta[i_surf][1](i_end, j_end),
                                  zeta[i_surf][2](i_end, j_end);
                            // position of the center point of the vortex filament
                            rp = 0.5*(r1 + r2);

                            // induced vel by vortices at vp
                            v_ind.setZero();
                            for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
                            {
                                UVLM::Types::VecMatrixX temp_uout;
                                UVLM::Types::allocate_VecMat(temp_uout,
                                                             zeta[ii_surf],
                                                             -1);
                                UVLM::BiotSavart::surface_with_steady_wake
                                (
                                    zeta[ii_surf],
                                    zeta_star[ii_surf],
                                    gamma[ii_surf],
                                    gamma_star[ii_surf],
                                    rp,
                                    options.horseshoe,
                                    temp_uout,
                                    options.ImageMethod,
                                    options.vortex_radius
                                );
                                v_ind(0) += temp_uout[0].sum();
                                v_ind(1) += temp_uout[1].sum();
                                v_ind(2) += temp_uout[2].sum();
                            }

                            dl = r2 - r1;

                            v << 0.5*(uext[i_surf][0](i_start, j_start) +
                                      uext[i_surf][0](i_end, j_end)),
                                 0.5*(uext[i_surf][1](i_start, j_start) +
                                      uext[i_surf][1](i_end, j_end)),
                                 0.5*(uext[i_surf][2](i_start, j_start) +
                                      uext[i_surf][2](i_end, j_end));

                            v = (v + v_ind).eval();

                            f = flightconditions.rho*gamma[i_surf](i_M, i_N)*v.cross(dl);

                            // transfer forces to matrix
                            // there are no moments
                            for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                            {
                                forces[i_surf][i_dim](i_start, j_start) +=
                                    0.5*f(i_dim);
                                forces[i_surf][i_dim](i_end, j_end) +=
                                    0.5*f(i_dim);
                            }
                        }
                    }
                }
            }
        }
		
        template <typename t_struct_nl_body>
        void calculate_static_forces_nonlifting_body
        (
            t_struct_nl_body& nl_body,
            const UVLM::Types::FlightConditions& flightconditions
        )
		{
            // Set forces to 0
            UVLM::Types::initialise_VecVecMat(nl_body.forces);

            // allocate 
            UVLM::Types::VecVecMatrixX velocities, forces_collocation;
            UVLM::Types::VecMatrixX pressure_coefficients, forces_norm, forces_norm_collocation;
            UVLM::Types::allocate_VecVecMat(velocities, nl_body.uext_col);            
            UVLM::Types::allocate_VecVecMat(forces_collocation, nl_body.uext_col);
			double freestream_velocity_squared;
			double induced_velocity_squared;
            UVLM::Types::allocate_VecMat(pressure_coefficients, nl_body.sigma);
            UVLM::Types::allocate_VecMat(forces_norm, nl_body.sigma, 1);
            UVLM::Types::allocate_VecMat(forces_norm_collocation, nl_body.sigma);

			UVLM::PostProc::get_local_panel_velocity(nl_body.normals,
													 nl_body.longitudinals,
													 nl_body.perpendiculars,
													 nl_body.uext_col,
													 nl_body.u_induced_col,
													 velocities
													 );
			const uint n_surf = nl_body.zeta.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const uint M = nl_body.sigma[i_surf].rows();
                const uint N = nl_body.sigma[i_surf].cols();
				
                for (uint i_M=0; i_M<M; ++i_M)
                {
                    for (uint i_N=0; i_N<N; ++i_N)
                    {
						freestream_velocity_squared = UVLM::Types::norm_Vec_squared(nl_body.uext_col[i_surf][0](i_M, i_N),
                                                                                    nl_body.uext_col[i_surf][1](i_M, i_N),
                                                                                    nl_body.uext_col[i_surf][2](i_M, i_N));
						induced_velocity_squared = UVLM::Types::norm_Vec_squared(velocities[i_surf][0](i_M, i_N),
                                                                                 velocities[i_surf][1](i_M, i_N), 
                                                                                 velocities[i_surf][2](i_M, i_N));                        
						
						nl_body.pressure_coefficients[i_surf](i_M, i_N) = 1.0 - induced_velocity_squared/freestream_velocity_squared;
						
						UVLM::Types::Real area = 0;
                        area = UVLM::Geometry::panel_area
                        (
                            nl_body.zeta[i_surf][0].template block<2,2>(i_M, i_N),
                            nl_body.zeta[i_surf][1].template block<2,2>(i_M, i_N),
                            nl_body.zeta[i_surf][2].template block<2,2>(i_M, i_N)
                        );
						// no moments considered 
						for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
						{
							forces_collocation[i_surf][i_dim](i_M, i_N) = 0.5*flightconditions.rho*freestream_velocity_squared
							                                  *nl_body.pressure_coefficients[i_surf](i_M, i_N)
															  *area*nl_body.normals[i_surf][i_dim](i_M, i_N);
						}
					}
				}
			}

            // export_data_to_csv_file("pressure_coefficient.csv", pressure_coefficients[0]);
            // Transform forces from collocation points to nodes
            UVLM::PostProc::transform_forces_from_col_to_nodes(forces_collocation, nl_body.forces, true);
            
        }

        template <typename t_zeta,
                  typename t_zeta_dot,
                  typename t_zeta_star,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_uext,
                  typename t_rbm_velocity,
                  typename t_centre_rot,
                  typename t_forces>
        void calculate_static_forces_unsteady
        (
            const t_zeta& zeta,
            const t_zeta_dot& zeta_dot,
            const t_zeta_star& zeta_star,
            const t_gamma& gamma,
            const t_gamma_star& gamma_star,
            const t_uext& uext,
            const t_rbm_velocity& rbm_velocity,
            const t_centre_rot& centre_rot,
            t_forces&  forces,
            const UVLM::Types::VMopts options,
            const UVLM::Types::FlightConditions& flightconditions
        )
        {
            // Set forces to 0
            UVLM::Types::initialise_VecVecMat(forces);

            // first calculate all the velocities at the corner points
            UVLM::Types::VecVecMatrixX velocities;
            UVLM::Types::allocate_VecVecMat(velocities, zeta);
            // free stream contribution
            UVLM::Types::copy_VecVecMat(uext, velocities);

            // u_ext taking into account unsteady contributions
            // TODO: get rid of rbm velocity input in function calculate static forces
            UVLM::Unsteady::Utils::compute_resultant_grid_velocity
            (
                zeta,
                zeta_dot,
                uext,
                rbm_velocity,
                centre_rot,
                velocities
            );
            // not bothered with effciency.
            // if it is so critical, it could be improved
            const uint n_surf = zeta.size();

            UVLM::Types::VecVecMatrixX span_seg_forces;
            UVLM::Types::VecVecMatrixX chord_seg_forces;

            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const uint M = gamma[i_surf].rows();
                const uint N = gamma[i_surf].cols();

                UVLM::Types::allocate_VecVecMat(span_seg_forces, 1, 3, M+1, N);
                UVLM::Types::allocate_VecVecMat(chord_seg_forces, 1, 3, M, N+1);

                // Computation of induced velocity in each segment
                #pragma omp parallel for collapse(2)
                for (uint i_M=0; i_M<M; ++i_M)
                {
                    for (uint i_N=0; i_N<N; ++i_N)
                    {
                        UVLM::Types::Vector3 dl;
                        UVLM::Types::Vector3 v;
                        UVLM::Types::Vector3 f;
                        UVLM::Types::Vector3 v_ind;
                        UVLM::Types::Vector3 rp;
                        UVLM::Types::Vector3 r1;
                        UVLM::Types::Vector3 r2;
                        UVLM::Types::Real delta_gamma;

                        // Spanwise vortices
                        r1 << zeta[i_surf][0](i_M, i_N),
                              zeta[i_surf][1](i_M, i_N),
                              zeta[i_surf][2](i_M, i_N);
                        r2 << zeta[i_surf][0](i_M, i_N+1),
                              zeta[i_surf][1](i_M, i_N+1),
                              zeta[i_surf][2](i_M, i_N+1);

                        // position of the center point of the vortex filament
                        rp = 0.5*(r1 + r2);

                        // induced vel by vortices at vp
                        v_ind.setZero();
                        for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
                        {
                            v_ind += UVLM::BiotSavart::whole_surface(zeta[ii_surf],
                                                                              gamma[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);

                            v_ind += UVLM::BiotSavart::whole_surface(zeta_star[ii_surf],
                                                                              gamma_star[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);
                        }

                        dl = r2-r1;

                        v << 0.5*(velocities[i_surf][0](i_M, i_N) +
                                  velocities[i_surf][0](i_M, i_N+1)),
                             0.5*(velocities[i_surf][1](i_M, i_N) +
                                  velocities[i_surf][1](i_M, i_N+1)),
                             0.5*(velocities[i_surf][2](i_M, i_N) +
                                  velocities[i_surf][2](i_M, i_N+1));

                        v = (v + v_ind).eval();
                         if (i_M == 0){
                            delta_gamma = -gamma[i_surf](i_M, i_N);
                        } else if (i_M == M){
                            // Might be needed if TE forces are computed
                            delta_gamma = gamma[i_surf](i_M-1, i_N);
                        } else {
                            delta_gamma = gamma[i_surf](i_M-1, i_N) - gamma[i_surf](i_M, i_N);
                        }

                        f = flightconditions.rho*delta_gamma*v.cross(dl);
                        span_seg_forces[0][0](i_M, i_N) = f(0);
                        span_seg_forces[0][1](i_M, i_N) = f(1);
                        span_seg_forces[0][2](i_M, i_N) = f(2);

                        // Chordwise vortice
                        r2 << zeta[i_surf][0](i_M+1, i_N),
                              zeta[i_surf][1](i_M+1, i_N),
                              zeta[i_surf][2](i_M+1, i_N);

                        // position of the center point of the vortex filament
                        rp = 0.5*(r1 + r2);

                        // induced vel by vortices at vp
                        v_ind.setZero();
                        for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
                        {
                            v_ind += UVLM::BiotSavart::whole_surface(zeta[ii_surf],
                                                                              gamma[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);

                            v_ind += UVLM::BiotSavart::whole_surface(zeta_star[ii_surf],
                                                                              gamma_star[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);
                        }

                        dl = r2-r1;

                        v << 0.5*(velocities[i_surf][0](i_M, i_N) +
                                  velocities[i_surf][0](i_M+1, i_N)),
                             0.5*(velocities[i_surf][1](i_M, i_N) +
                                  velocities[i_surf][1](i_M+1, i_N)),
                             0.5*(velocities[i_surf][2](i_M, i_N) +
                                  velocities[i_surf][2](i_M+1, i_N));

                        v = (v + v_ind).eval();

                        if (i_N == 0){
                            delta_gamma = gamma[i_surf](i_M, i_N);
                        } else if (i_N == N){
                            delta_gamma = -gamma[i_surf](i_M, i_N-1);
                        } else {
                            delta_gamma = gamma[i_surf](i_M, i_N) - gamma[i_surf](i_M, i_N-1);
                        }

                        f = flightconditions.rho*delta_gamma*v.cross(dl);
                        chord_seg_forces[0][0](i_M, i_N) = f(0);
                        chord_seg_forces[0][1](i_M, i_N) = f(1);
                        chord_seg_forces[0][2](i_M, i_N) = f(2);
                    }
                }

                // Influence of the last chordwise column of vortices
                UVLM::Types::Vector3 dl;
                UVLM::Types::Vector3 v;
                UVLM::Types::Vector3 f;
                UVLM::Types::Vector3 v_ind;
                UVLM::Types::Vector3 rp;
                UVLM::Types::Vector3 r1;
                UVLM::Types::Vector3 r2;
                UVLM::Types::Real delta_gamma;
                for (uint i_M=0; i_M<M; ++i_M){

                    r1 << zeta[i_surf][0](i_M, N),
                          zeta[i_surf][1](i_M, N),
                          zeta[i_surf][2](i_M, N);
                    r2 << zeta[i_surf][0](i_M+1, N),
                          zeta[i_surf][1](i_M+1, N),
                          zeta[i_surf][2](i_M+1, N);

                    // position of the center point of the vortex filament
                    rp = 0.5*(r1 + r2);

                    // induced vel by vortices at vp
                    v_ind.setZero();

                    for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
                    {
                        v_ind += UVLM::BiotSavart::whole_surface(zeta[ii_surf],
                                                                          gamma[ii_surf],
                                                                          rp,
                                                                          options.ImageMethod,
                                                                          options.vortex_radius);

                        v_ind += UVLM::BiotSavart::whole_surface(zeta_star[ii_surf],
                                                                          gamma_star[ii_surf],
                                                                          rp,
                                                                          options.ImageMethod,
                                                                          options.vortex_radius);
                    }

                    dl = r2-r1;

                    v << 0.5*(velocities[i_surf][0](i_M, N) +
                              velocities[i_surf][0](i_M+1, N)),
                         0.5*(velocities[i_surf][1](i_M, N) +
                              velocities[i_surf][1](i_M+1, N)),
                         0.5*(velocities[i_surf][2](i_M, N) +
                              velocities[i_surf][2](i_M+1, N));

                    v = (v + v_ind).eval();

                    delta_gamma = -gamma[i_surf](i_M, N-1);
                    f = flightconditions.rho*delta_gamma*v.cross(dl);
                    chord_seg_forces[0][0](i_M, N) = f(0);
                    chord_seg_forces[0][1](i_M, N) = f(1);
                    chord_seg_forces[0][2](i_M, N) = f(2);

                }

                // #pragma omp parallel for collapse(2) reduction(sum_Vector3: uout)
                // Transfer forces to nodes
                transfer_forces_to_nodes(
                    forces,
                    span_seg_forces,
                    chord_seg_forces,
                    N,
                    M,
                    i_surf,
                    options.ignore_first_x_nodes_in_force_calculation
                );
            }
        }
        template <typename t_forces,
                  typename t_span_seg_forces,
                  typename t_chord_seg_forces>
        void transfer_forces_to_nodes
        (
            t_forces& forces,
            t_span_seg_forces& span_seg_forces,
            t_chord_seg_forces& chord_seg_forces,
            const uint N,
            const uint M,
            const uint i_surf,
            const uint ignore_first_x_nodes_in_force_calculation
        )
        {
            
            for (uint i_M=0; i_M<M+1; ++i_M)
            {
                for (uint i_N=ignore_first_x_nodes_in_force_calculation; i_N<N+1; ++i_N)
                {
                    for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                    {
                        // Spanwise segments
                        if (i_N != 0){
                            forces[i_surf][i_dim](i_M, i_N) += 0.5*span_seg_forces[0][i_dim](i_M, i_N-1);
                        }
                        if (i_N != N){
                            forces[i_surf][i_dim](i_M, i_N) += 0.5*span_seg_forces[0][i_dim](i_M, i_N);
                        }

                        // Chordwise segments
                        if (i_M != 0){
                            forces[i_surf][i_dim](i_M, i_N) += 0.5*chord_seg_forces[0][i_dim](i_M-1, i_N);
                        }
                        if (i_M != M){
                            forces[i_surf][i_dim](i_M, i_N) += 0.5*chord_seg_forces[0][i_dim](i_M, i_N);
                        }
                    }
                }
            }
            
        }
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
        template <typename t_corner_points,
                  typename t_nl_body>
        void calculate_source_induced_velocities_on_points
        (
            const t_corner_points& corner_points,
            const t_nl_body& nl_body,
            UVLM::Types::VectorX& sigma_flat,
            UVLM::Types::VecVecMatrixX& u_induced_by_sources
        )
        {
            // Get surface vectors of zeta_star (account for points instead of panels)
            // determine convection velocity u_ind from non lifting surfaces
            UVLM::Types::VecVecMatrixX normals, longitudinals, perpendiculars;
            
            UVLM::Types::allocate_VecVecMat(normals, corner_points);
            UVLM::Types::allocate_VecVecMat(longitudinals, corner_points);
            UVLM::Types::allocate_VecVecMat(perpendiculars, corner_points);
            UVLM::Geometry::generate_surface_vectors_wake(corner_points, normals, longitudinals, perpendiculars);
                
            // Allocate matrices for source influence
            uint Ktotal = UVLM::Matrix::get_total_VecVecMat_size(normals);
            UVLM::Types::MatrixX aic_x = UVLM::Types::MatrixX::Zero(Ktotal, nl_body.Ktotal);
            UVLM::Types::MatrixX aic_y = UVLM::Types::MatrixX::Zero(Ktotal, nl_body.Ktotal);
            UVLM::Types::MatrixX aic_z = UVLM::Types::MatrixX::Zero(Ktotal, nl_body.Ktotal);

            // Get induced velocities by sources
            UVLM::Matrix::AIC_sources(nl_body.zeta,
                                        corner_points,
                                        nl_body.longitudinals,
                                        nl_body.perpendiculars,
                                        nl_body.normals,
                                        longitudinals,
                                        perpendiculars,
                                        normals,
                                        aic_x,
                                        aic_y,
                                        aic_z);
	        calculate_induced_velocity_col(sigma_flat,
                                            aic_x,
                                            aic_y,
                                            aic_z,
                                            u_induced_by_sources);
        }

        // TODO: Merge this function with the one above
        template <typename t_struct_lifting_surfaces,
                  typename t_struct_phantom>
        void calculate_static_forces_unsteady
        (
            t_struct_lifting_surfaces& lifting_surfaces,
            const t_struct_phantom& phantom_surfaces,
            const UVLM::Types::VMopts options,
            const UVLM::Types::FlightConditions& flightconditions
        )
        {
            // Set forces to 0
            UVLM::Types::initialise_VecVecMat(lifting_surfaces.forces);

            // first calculate all the velocities at the corner points
            UVLM::Types::VecVecMatrixX velocities;
            UVLM::Types::allocate_VecVecMat(velocities, lifting_surfaces.zeta);
            // free stream contribution
            UVLM::Types::copy_VecVecMat(lifting_surfaces.u_ext, velocities);

            // u_ext taking into account unsteady contributions
            // TODO: get rid of rbm velocity input in function calculate static forces
            UVLM::Unsteady::Utils::compute_resultant_grid_velocity
            (
                lifting_surfaces.zeta,
                lifting_surfaces.zeta_dot,
                lifting_surfaces.u_ext,
                lifting_surfaces.rbm_vel_g,
                lifting_surfaces.centre_rot,
                velocities
            );

            // not bothered with effciency.
            // if it is so critical, it could be improved
            const uint n_surf = lifting_surfaces.zeta.size();

            UVLM::Types::VecVecMatrixX span_seg_forces;
            UVLM::Types::VecVecMatrixX chord_seg_forces;

            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const uint M = lifting_surfaces.gamma[i_surf].rows();
                const uint N = lifting_surfaces.gamma[i_surf].cols();

                UVLM::Types::allocate_VecVecMat(span_seg_forces, 1, 3, M+1, N);
                UVLM::Types::allocate_VecVecMat(chord_seg_forces, 1, 3, M, N+1);

                // Computation of induced velocity in each segment
                #pragma omp parallel for collapse(2)
                for (uint i_M=0; i_M<M; ++i_M)
                {
                    for (uint i_N=0; i_N<N; ++i_N)
                    {
                        UVLM::Types::Vector3 dl;
                        UVLM::Types::Vector3 v;
                        UVLM::Types::Vector3 f;
                        UVLM::Types::Vector3 v_ind;
                        UVLM::Types::Vector3 rp;
                        UVLM::Types::Vector3 r1;
                        UVLM::Types::Vector3 r2;
                        UVLM::Types::Real delta_gamma;

                        // Spanwise vortices
                        r1 << lifting_surfaces.zeta[i_surf][0](i_M, i_N),
                              lifting_surfaces.zeta[i_surf][1](i_M, i_N),
                              lifting_surfaces.zeta[i_surf][2](i_M, i_N);
                        r2 << lifting_surfaces.zeta[i_surf][0](i_M, i_N+1),
                              lifting_surfaces.zeta[i_surf][1](i_M, i_N+1),
                              lifting_surfaces.zeta[i_surf][2](i_M, i_N+1);

                        // position of the center point of the vortex filament
                        rp = 0.5*(r1 + r2);

                        // induced vel by vortices at vp
                        v_ind.setZero();
                        for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
                        {
                            v_ind += UVLM::BiotSavart::whole_surface(lifting_surfaces.zeta[ii_surf],
                                                                              lifting_surfaces.gamma[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);

                            v_ind += UVLM::BiotSavart::whole_surface(lifting_surfaces.zeta_star[ii_surf],
                                                                              lifting_surfaces.gamma_star[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);
                            
                            v_ind += UVLM::BiotSavart::whole_surface(phantom_surfaces.zeta[ii_surf],
                                                                              phantom_surfaces.gamma[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);

                            v_ind += UVLM::BiotSavart::whole_surface(phantom_surfaces.zeta_star[ii_surf],
                                                                              phantom_surfaces.gamma_star[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);
                                                                              
                                                                    
                        }
                        if ((!options.phantom_wing_test) && (options.consider_u_ind_by_sources_for_lifting_forces))
                        {
                            for (uint idim=0; idim < UVLM::Constants::NDIM; idim++)
                            {                 
                                v_ind(idim) += lifting_surfaces.u_induced_by_sources_on_center_spanwise_vertices[i_surf][idim](i_M,i_N);
                            }

                        }

                        dl = r2-r1;

                        v << 0.5*(velocities[i_surf][0](i_M, i_N) +
                                  velocities[i_surf][0](i_M, i_N+1)),
                             0.5*(velocities[i_surf][1](i_M, i_N) +
                                  velocities[i_surf][1](i_M, i_N+1)),
                             0.5*(velocities[i_surf][2](i_M, i_N) +
                                  velocities[i_surf][2](i_M, i_N+1));

                        v = (v + v_ind).eval();
                        
                        if (i_M == 0){
                            delta_gamma = -lifting_surfaces.gamma[i_surf](i_M, i_N);
                        } else if (i_M == M){
                            // Might be needed if TE forces are computed
                            delta_gamma = lifting_surfaces.gamma[i_surf](i_M-1, i_N);
                        } else {
                            delta_gamma = lifting_surfaces.gamma[i_surf](i_M-1, i_N) - lifting_surfaces.gamma[i_surf](i_M, i_N);
                        }

                        f = flightconditions.rho*delta_gamma*v.cross(dl);
                        span_seg_forces[0][0](i_M, i_N) = f(0);
                        span_seg_forces[0][1](i_M, i_N) = f(1);
                        span_seg_forces[0][2](i_M, i_N) = f(2);

                        // Chordwise vortice
                        r2 << lifting_surfaces.zeta[i_surf][0](i_M+1, i_N),
                              lifting_surfaces.zeta[i_surf][1](i_M+1, i_N),
                              lifting_surfaces.zeta[i_surf][2](i_M+1, i_N);

                        // position of the center point of the vortex filament
                        rp = 0.5*(r1 + r2);

                        // induced vel by vortices at vp
                        v_ind.setZero();
                        for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
                        {
                            v_ind += UVLM::BiotSavart::whole_surface(lifting_surfaces.zeta[ii_surf],
                                                                              lifting_surfaces.gamma[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);

                            v_ind += UVLM::BiotSavart::whole_surface(lifting_surfaces.zeta_star[ii_surf],
                                                                              lifting_surfaces.gamma_star[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);
                                                                              
                            v_ind += UVLM::BiotSavart::whole_surface(phantom_surfaces.zeta[ii_surf],
                                                                              phantom_surfaces.gamma[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);

                            v_ind += UVLM::BiotSavart::whole_surface(phantom_surfaces.zeta_star[ii_surf],
                                                                              phantom_surfaces.gamma_star[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);
                                                                              
                        }
                                                
                        if ((!options.phantom_wing_test) && (options.consider_u_ind_by_sources_for_lifting_forces))
                        {
                            for (uint idim=0; idim < UVLM::Constants::NDIM; idim++)
                            {
                            
                                v_ind(idim) += lifting_surfaces.u_induced_by_sources_on_center_chordwise_vertices[i_surf][idim](i_M,i_N);
                            }
                        }

                        dl = r2-r1;

                        v << 0.5*(velocities[i_surf][0](i_M, i_N) +
                                  velocities[i_surf][0](i_M+1, i_N)),
                             0.5*(velocities[i_surf][1](i_M, i_N) +
                                  velocities[i_surf][1](i_M+1, i_N)),
                             0.5*(velocities[i_surf][2](i_M, i_N) +
                                  velocities[i_surf][2](i_M+1, i_N));

                        v = (v + v_ind).eval();

                        if (i_N == 0){
                            delta_gamma = lifting_surfaces.gamma[i_surf](i_M, i_N);
                        } else if (i_N == N){
                            delta_gamma = -lifting_surfaces.gamma[i_surf](i_M, i_N-1);
                        } else {
                            delta_gamma = lifting_surfaces.gamma[i_surf](i_M, i_N) - lifting_surfaces.gamma[i_surf](i_M, i_N-1);
                        }

                        f = flightconditions.rho*delta_gamma*v.cross(dl);
                        chord_seg_forces[0][0](i_M, i_N) = f(0);
                        chord_seg_forces[0][1](i_M, i_N) = f(1);
                        chord_seg_forces[0][2](i_M, i_N) = f(2);
                    }
                }

                // Influence of the last chordwise column of vortices
                UVLM::Types::Vector3 dl;
                UVLM::Types::Vector3 v;
                UVLM::Types::Vector3 f;
                UVLM::Types::Vector3 v_ind;
                UVLM::Types::Vector3 rp;
                UVLM::Types::Vector3 r1;
                UVLM::Types::Vector3 r2;
                UVLM::Types::Real delta_gamma;
                for (uint i_M=0; i_M<M; ++i_M){

                    r1 << lifting_surfaces.zeta[i_surf][0](i_M, N),
                          lifting_surfaces.zeta[i_surf][1](i_M, N),
                          lifting_surfaces.zeta[i_surf][2](i_M, N);
                    r2 << lifting_surfaces.zeta[i_surf][0](i_M+1, N),
                          lifting_surfaces.zeta[i_surf][1](i_M+1, N),
                          lifting_surfaces.zeta[i_surf][2](i_M+1, N);

                    // position of the center point of the vortex filament
                    rp = 0.5*(r1 + r2);

                    // induced vel by vortices at vp
                    v_ind.setZero();

                    for (uint ii_surf=0; ii_surf<n_surf; ++ii_surf)
                    {
                        v_ind += UVLM::BiotSavart::whole_surface(lifting_surfaces.zeta[ii_surf],
                                                                          lifting_surfaces.gamma[ii_surf],
                                                                          rp,
                                                                          options.ImageMethod,
                                                                          options.vortex_radius);

                        v_ind += UVLM::BiotSavart::whole_surface(lifting_surfaces.zeta_star[ii_surf],
                                                                          lifting_surfaces.gamma_star[ii_surf],
                                                                          rp,
                                                                          options.ImageMethod,
                                                                          options.vortex_radius);
                            v_ind += UVLM::BiotSavart::whole_surface(phantom_surfaces.zeta[ii_surf],
                                                                              phantom_surfaces.gamma[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);

                            v_ind += UVLM::BiotSavart::whole_surface(phantom_surfaces.zeta_star[ii_surf],
                                                                              phantom_surfaces.gamma_star[ii_surf],
                                                                              rp,
                                                                              options.ImageMethod,
                                                                              options.vortex_radius);
                                                                              
                                                                          
                    }

                    dl = r2-r1;

                    v << 0.5*(velocities[i_surf][0](i_M, N) +
                              velocities[i_surf][0](i_M+1, N)),
                         0.5*(velocities[i_surf][1](i_M, N) +
                              velocities[i_surf][1](i_M+1, N)),
                         0.5*(velocities[i_surf][2](i_M, N) +
                              velocities[i_surf][2](i_M+1, N));

                    v = (v + v_ind).eval();

                    delta_gamma = -lifting_surfaces.gamma[i_surf](i_M, N-1);
                    f = flightconditions.rho*delta_gamma*v.cross(dl);
                    chord_seg_forces[0][0](i_M, N) = f(0);
                    chord_seg_forces[0][1](i_M, N) = f(1);
                    chord_seg_forces[0][2](i_M, N) = f(2);

                }

                // #pragma omp parallel for collapse(2) reduction(sum_Vector3: uout)
                // Transfer forces to nodes
                transfer_forces_to_nodes(
                    lifting_surfaces.forces,
                    span_seg_forces,
                    chord_seg_forces,
                    N,
                    M,
                    i_surf,
                    options.ignore_first_x_nodes_in_force_calculation
                );
            }
        }

        // Forces is not set to 0, forces are added
        template <typename t_zeta,
                  typename t_zeta_star,
                  typename t_zeta_col,
                  typename t_gamma,
                  typename t_gamma_star,
                  typename t_gamma_dot,
                  typename t_normals,
                  typename t_forces>
        void calculate_dynamic_forces
        (
            const t_zeta& zeta,
            const t_zeta_star& zeta_star,
            const t_zeta_col& zeta_col,
            const t_gamma& gamma,
            const t_gamma_star& gamma_star,
            const t_gamma_dot& gamma_dot,
            const t_normals& normals,
            t_forces&  forces,
            const UVLM::Types::UVMopts options,
            const UVLM::Types::FlightConditions& flightconditions
        )
        {
            const UVLM::Types::Real dt = options.dt;
            const uint n_surf = zeta.size();

            UVLM::Types::VecVecMatrixX unsteady_force;
            UVLM::Types::allocate_VecVecMat(unsteady_force, forces, -1);

            // calculate unsteady forces
            // f_uns = rho*A*n*gamma_dot
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const uint n_rows = gamma[i_surf].rows();
                const uint n_cols = gamma[i_surf].cols();
                for (uint i=0; i<n_rows; ++i)
                {
                    for (uint j=0; j<n_cols; ++j)
                    {
                        // area calculation
                        UVLM::Types::Real area = 0;
                        area = UVLM::Geometry::panel_area
                        (
                            zeta[i_surf][0].template block<2,2>(i, j),
                            zeta[i_surf][1].template block<2,2>(i, j),
                            zeta[i_surf][2].template block<2,2>(i, j)
                        );

                        // rho*A*n*gamma_dot
                        for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                        {
                            unsteady_force[i_surf][i_dim](i, j) =
                            (
                                // 0.0*flightconditions.rho
                                // -flightconditions.rho
                                -flightconditions.rho
                               *area
                               *normals[i_surf][i_dim](i, j)
                               *gamma_dot[i_surf](i, j)
                            );
                        }

                        // transfer forces to vortex corners
                        UVLM::Types::Vector3 zeta_col_panel;
                        zeta_col_panel << zeta_col[i_surf][0](i, j),
                                          zeta_col[i_surf][1](i, j),
                                          zeta_col[i_surf][2](i, j);
                        UVLM::Types::Vector3 panel_force;
                        for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                        {
                            panel_force(i_dim) = unsteady_force[i_surf][i_dim](i, j);
                        }

                        for (uint ii=0; ii<2; ++ii)
                        {
                            if ((ii == 1) && (i == n_rows - 1))
                            {
                                // trailing edge
                                continue;
                            }
                            for (uint jj=0; jj<2; ++jj)
                            {
                                // forces
                                for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                                {
                                    forces[i_surf][i_dim](i + ii, j + jj) +=
                                        0.25*panel_force(i_dim);
                                }
                                // moments
                                // moment = r cross F
                                // UVLM::Types::Vector3 zeta_corner;
                                // zeta_corner << zeta[i_surf][0](i + ii, j + jj),
                                //                zeta[i_surf][1](i + ii, j + jj),
                                //                zeta[i_surf][2](i + ii, j + jj);

                                // UVLM::Types::Vector3 r;
                                // r = zeta_corner - zeta_col_panel;
                                // UVLM::Types::Vector3 moment;
                                // moment = 0.25*r.cross(panel_force);
                                //
                                // for (uint i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                                // {
                                //     uint i_moment = i_dim + 3;
                                //     forces[i_surf][i_moment](i + ii, j + jj) +=
                                //         moment(i_dim);
                                // }
                            }
                        }
                    }
                }
            }
        }

        // Forces is not set to 0, forces are added
        template <typename t_u_ext,
                  typename t_zeta,
                  typename t_zeta_dot,
                  typename t_normals,
                  typename t_rbm_vel,
                  typename t_incidence_angle>
        void calculate_incidence_angle
        (
            const t_u_ext& u_ext,
            const t_zeta& zeta,
            const t_zeta_dot& zeta_dot,
            const t_normals& normals,
            const t_rbm_vel& rbm_vel_g,
            const UVLM::Types::UVMopts& options,
            t_incidence_angle& incidence_angle
        )
        {
            // compute instantaneous velocity for every panel
            UVLM::Types::VecVecMatrixX velocities;
            UVLM::Types::allocate_VecVecMat(velocities, zeta);
            // free stream contribution
            UVLM::Types::copy_VecVecMat(u_ext, velocities);

            UVLM::Types::Vector3 centre_rot = UVLM::Types::Vector3::Zero();
            // u_ext taking into account unsteady contributions
            UVLM::Unsteady::Utils::compute_resultant_grid_velocity
            (
                zeta,
                zeta_dot,
                u_ext,
                rbm_vel_g,
                centre_rot,
                velocities
            );

            // Currently, the incidence computed is the one based on external velocities
            // thus this function provides something like the geometric angle of attack.
            // If we wanted the effective angle of attack (lifting line theory nomenclature)
            // we would need to include the induced velocity, I think.
            // It is not efficient to code it here so I will output this velocities directly
            // from the solver

            // stall angle will be computed as the angle between
            // the velocity at the leading edge and the chord line
            const uint n_surf = zeta.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const uint M = zeta[i_surf][0].rows() - 1;
                const uint N = zeta[i_surf][0].cols() - 1;

                // loop through spanwise panels
                for (uint i_N=0; i_N<N; ++i_N)
                {
                    UVLM::Types::Vector3 chord_line_g;
                    UVLM::Types::Vector3 local_vel;

                    // chord line
                    // trailing edge
                    chord_line_g(0) = 0.5*(zeta[i_surf][0](M, i_N) +
                                           zeta[i_surf][0](M, i_N + 1));
                    chord_line_g(1) = 0.5*(zeta[i_surf][1](M, i_N) +
                                           zeta[i_surf][1](M, i_N + 1));
                    chord_line_g(2) = 0.5*(zeta[i_surf][2](M, i_N) +
                                           zeta[i_surf][2](M, i_N + 1));

                    chord_line_g(0) -= 0.5*(zeta[i_surf][0](0, i_N) +
                                            zeta[i_surf][0](0, i_N + 1));
                    chord_line_g(1) -= 0.5*(zeta[i_surf][1](0, i_N) +
                                            zeta[i_surf][1](0, i_N + 1));
                    chord_line_g(2) -= 0.5*(zeta[i_surf][2](0, i_N) +
                                            zeta[i_surf][2](0, i_N + 1));

                    // velocity at the center of the vortex segment
                    // is avg of both vertices
                    local_vel(0) = 0.5*(velocities[i_surf][0](0, i_N) +
                                        velocities[i_surf][0](0, i_N + 1));
                    local_vel(1) = 0.5*(velocities[i_surf][1](0, i_N) +
                                        velocities[i_surf][1](0, i_N + 1));
                    local_vel(2) = 0.5*(velocities[i_surf][2](0, i_N) +
                                        velocities[i_surf][2](0, i_N + 1));

                    // we consider angle of attack the angle
                    // between these two vectors projected in the
                    // vertical plane containing the chord line

                    UVLM::Types::Vector3 plane_normal;
                    UVLM::Types::Vector3 vertical;
                    vertical << normals[i_surf][0](0, i_N),
                                normals[i_surf][1](0, i_N),
                                normals[i_surf][2](0, i_N);

                    // check that normal is up, if not, change the sign
                    // UVLM::Types::Real sign = 1;
                    // if (vertical(2) < 0)
                    // {
                    //     sign = -1;
                    // }

                    plane_normal = vertical.cross(chord_line_g).normalized();

                    UVLM::Types::Vector3 local_vel_projected;
                    local_vel_projected = local_vel -
                        (local_vel.dot(plane_normal))*plane_normal;

                    // UVLM::Types::Real angle;
                    // angle = std::atan2(
                    //     plane_normal.dot(chord_line_g.cross(local_vel_projected)),
                    //     chord_line_g.dot(local_vel_projected)
                    // );

                    // I think this is more intuitive, they give the same result
                    UVLM::Types::Real angle;
                    angle = std::acos(
                        chord_line_g.dot(local_vel_projected)/
                        chord_line_g.norm()/local_vel_projected.norm()
                    );

                    UVLM::Types::Real sign = 1.;
                    if (local_vel_projected.dot(vertical) < 0)
                    {
                        sign = -1.;
                    }

                    for (uint i_M=0; i_M<M; ++i_M)
                    {
                        incidence_angle[i_surf](i_M, i_N) = -angle*sign;
                    }
                }
            }
        }

        template <typename t_normal,
				  typename t_longitudinal,
				  typename t_perpendicular,
				  typename t_uext_col,
				  typename t_u_induced,
				  typename t_u_total_panel>
		void get_local_panel_velocity
		(
		const t_normal& normal,
		const t_longitudinal& longitudinal,
		const t_perpendicular& perpendicular,
        const t_uext_col& uext_col,
        const t_u_induced& u_induced_col,
        t_u_total_panel& u_total_panel
		)
		{
			const uint n_surf = uext_col.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
			{
				const uint M = uext_col[i_surf][0].rows();
                const uint N = uext_col[i_surf][0].cols();
                for (uint i_col=0; i_col<M; ++i_col)
                {
	                for (uint j_col=0; j_col<N; ++j_col)
					{
			            UVLM::Types::Vector3 vector_trans_tmp;
			            UVLM::Types::Vector3 u_induced_col_global = UVLM::Types::Vector3(u_induced_col[i_surf][0](i_col, j_col),
																						 u_induced_col[i_surf][1](i_col, j_col),
																						 u_induced_col[i_surf][2](i_col, j_col));
						UVLM::Types::Vector3 longitudinal_panel = UVLM::Types::Vector3(longitudinal[i_surf][0](i_col, j_col), longitudinal[i_surf][1](i_col, j_col), longitudinal[i_surf][2](i_col, j_col));
						UVLM::Types::Vector3 perpendicular_panel = UVLM::Types::Vector3(perpendicular[i_surf][0](i_col, j_col), perpendicular[i_surf][1](i_col, j_col), perpendicular[i_surf][2](i_col, j_col));
			            UVLM::Types::Vector3 normal_panel = UVLM::Types::Vector3(normal[i_surf][0](i_col, j_col), normal[i_surf][1](i_col, j_col), normal[i_surf][2](i_col, j_col));
						UVLM::Geometry::convert_to_global_coordinate_system(u_induced_col_global,
															longitudinal_panel,
															perpendicular_panel,
															normal_panel);
						for (uint dim=0; dim<UVLM::Constants::NDIM; dim++)
						{
							u_total_panel[i_surf][dim](i_col, j_col)= u_induced_col_global[dim] + uext_col[i_surf][dim](i_col, j_col);
						}
						
					}
				}
			}
        }

        template <typename t_normal,
				  typename t_longitudinal,
				  typename t_perpendicular,
				  typename t_u_induced>
		void transform_VecVecMatrix_to_global_coordinate_system
		(
		const t_normal& normal,
		const t_longitudinal& longitudinal,
		const t_perpendicular& perpendicular,
        t_u_induced& u_induced_col
		)
		{
            // ToDo: Outsource to geometry? integrate in function above
			const uint n_surf = u_induced_col.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
			{
				const uint M = u_induced_col[i_surf][0].rows();
                const uint N = u_induced_col[i_surf][0].cols();
                for (uint i_col=0; i_col<M; ++i_col)
                {
	                for (uint j_col=0; j_col<N; ++j_col)
					{
			            UVLM::Types::Vector3 vector_trans_tmp;
			            UVLM::Types::Vector3 u_induced_col_global = UVLM::Types::Vector3(u_induced_col[i_surf][0](i_col, j_col),
																						 u_induced_col[i_surf][1](i_col, j_col),
																						 u_induced_col[i_surf][2](i_col, j_col));
						UVLM::Types::Vector3 longitudinal_panel = UVLM::Types::Vector3(longitudinal[i_surf][0](i_col, j_col), longitudinal[i_surf][1](i_col, j_col), longitudinal[i_surf][2](i_col, j_col));
						UVLM::Types::Vector3 perpendicular_panel = UVLM::Types::Vector3(perpendicular[i_surf][0](i_col, j_col), perpendicular[i_surf][1](i_col, j_col), perpendicular[i_surf][2](i_col, j_col));
			            UVLM::Types::Vector3 normal_panel = UVLM::Types::Vector3(normal[i_surf][0](i_col, j_col), normal[i_surf][1](i_col, j_col), normal[i_surf][2](i_col, j_col));
						UVLM::Geometry::convert_to_global_coordinate_system(u_induced_col_global,
															longitudinal_panel,
															perpendicular_panel,
															normal_panel);
						for (uint dim=0; dim<UVLM::Constants::NDIM; dim++)
						{
							u_induced_col[i_surf][dim](i_col, j_col)= u_induced_col_global[dim];
						}
						
					}
				}
			}
		}

        template <typename t_forces_collocation,
				  typename t_forces_node>
		void transform_forces_from_col_to_nodes
		(
            const t_forces_collocation& forces_collocation,
            t_forces_node& forces_node,
            const bool rotational_body
		)
		{
			const uint n_surf = forces_node.size();
            uint considered_columns, column_shift;
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
			{
				uint M = forces_node[i_surf][0].rows();
                uint N = forces_node[i_surf][0].cols();


                for (uint i_row=0; i_row<M; ++i_row)
                {
                    if (rotational_body && i_row == M-1)
                    {
                        // skips last row for rotational bodies where zeta's last row contains the position 
                        // of the nodes already contained in the first row
                        break;
                    }
                    else
                    {
                        for (uint j_col=0; j_col<N; ++j_col)
                        {	
                            
                            if ((j_col ==0 || j_col ==N-1) && (rotational_body)) 
                            {
                                // skips first and last column as for rotanional body the force at the sharp edge should be zero
                                continue;
                            }
                            if (j_col == 0) 
                            {
                                considered_columns = 1;
                            column_shift = 0;
                        }
                        else if (j_col ==N-1)
                        {
                            considered_columns = 1;
                            column_shift = 1;
                        }
                        else
                        {
                            considered_columns = 2;
                            column_shift = 1;
                        }
                        for (uint dim=0; dim<UVLM::Constants::NDIM; dim++)
						{

                            if (i_row == M-2)
                            {
                                forces_node[i_surf][dim](i_row, j_col) = forces_collocation[i_surf][dim].block(i_row, j_col-column_shift, 1, considered_columns).sum();
                                forces_node[i_surf][dim](i_row, j_col) += forces_collocation[i_surf][dim].block(0, j_col-column_shift, 1, considered_columns).sum();
                                forces_node[i_surf][dim](i_row, j_col) *= 0.25;
                                }
                                else if (i_row == M-1)
                                {
                                    forces_node[i_surf][dim](i_row, j_col) = forces_node[i_surf][dim](0, j_col);
                                }
                                else
                                {
                                    forces_node[i_surf][dim](i_row, j_col) = forces_collocation[i_surf][dim].block(i_row, j_col-column_shift, 2, considered_columns).mean();
                                }	
                            }						
                        }
                    }
				}
			}
		}       
        template <typename t_forces_node,
				  typename t_forces_norm>
		void get_norm_forces
		(
            const t_forces_node& forces_node,
            t_forces_norm& forces_norm
		)
		{
            UVLM::Types::Vector3 vec_force;
            
			const uint n_surf = forces_node.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
			{
				const uint M = forces_node[i_surf][0].rows();
                const uint N = forces_node[i_surf][0].cols();
                
				const uint M_norm = forces_norm[i_surf].rows();
                const uint N_norm = forces_norm[i_surf].cols();
                for (uint i_row=0; i_row<M_norm; ++i_row)
                {
	                for (uint j_col=0; j_col<N_norm; ++j_col)
					{		

                        vec_force << forces_node[i_surf][0](i_row, j_col)
                                     ,forces_node[i_surf][1](i_row, j_col)
                                     ,forces_node[i_surf][2](i_row, j_col);
                        forces_norm[i_surf](i_row, j_col) = vec_force.norm();
					}
				}
			}
		}
    }
}
