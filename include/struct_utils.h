/**
 * @file struct_utils.h
 * @brief This file contains C++ structures for handling various aerodynamic surface types.
 */

#include "types.h"
#include "geometry.h"
#include "mapping.h"
#include "matrix.h"
#include "phantom.h"

#include <iostream>

namespace UVLM
{
    namespace StructUtils
    {
        /**
         * @struct surface
         * @brief Represents a basic aerodynamic surface.
        */
        struct surface
        {
            uint n_surf; // Number of surfaces.
            uint Ktotal; // Number of total collocation points.
            UVLM::Types::VecDimensions dimensions; // An array of dimensions (chordwise and spanwise) for each surface.
            UVLM::Types::VecVecMapX zeta; // A mapped matrix containing corner point coordinates for each surface grid.
            UVLM::Types::VecVecMapX u_ext; // A mapped matrix containing external velocities at each surface corner point.
            UVLM::Types::VecVecMapX forces; // A mapped matrix containing forces at each surface corner point.
                   
            UVLM::Types::VecVecMatrixX zeta_col; // A  matix containing collocation point coordinates located on each surface grid.
            UVLM::Types::VecVecMatrixX uext_col; // A matrix containing external flow velocities (free flow and gust) on the collocation points of each surface.
            UVLM::Types::VecVecMatrixX uext_total;// A  matrix containing all (flow, structural, rigid) external velocities at each surface corner point.
            UVLM::Types::VecVecMatrixX uext_total_col; // A matrix containing all (flow, structural, rigid) external velocities on the collocation points of each surface.
            UVLM::Types::VecVecMatrixX normals; // A matrix containing the normal vectors for the panels of each surface.
            UVLM::Types::VecVecMatrixX longitudinals; // A matrix containing the longitudinal vectors of the panels of each surface.
            UVLM::Types::VecVecMatrixX perpendiculars; //  A matrix containing the perpendicular vectors for the panels of each surface.
            UVLM::Types::VectorX rhs; // Boundary condition vector (UVLM B-Matrix) aka little right hand side.
            /**
             * @brief Constructor for a basic aerodynamic surface.
             *
             * @param n_surfaces Number of surfaces.
             * @param p_dimensions Pointer to array containing dimensions (chordwise and spanwise) for each surface.
             * @param p_zeta Pointer to matrix containing corner point coordinates for each surface grid.
             * @param p_u_ext Pointer to matrix containing external flow velocities (free flow and gust) at corner points of each surface grid.
             * @param p_forces Pointer to matrix containing forces at corner points of each surface grid.
             */
            surface 
            (
                uint n_surfaces,
                unsigned int** p_dimensions,
                double** p_zeta,
                double** p_u_ext,
                double** p_forces
            )
            {
                n_surf = n_surfaces;                
                UVLM::Mapping::transform_dimensions(n_surf,p_dimensions, dimensions);
                UVLM::Mapping::map_VecVecMat(dimensions, p_zeta, zeta, 1);  
                UVLM::Mapping::map_VecVecMat(dimensions, p_u_ext, u_ext, 1);   
                UVLM::Mapping::map_VecVecMat(dimensions, p_forces, forces, 1, 2*UVLM::Constants::NDIM);  
            }

            /**
             * @brief Get surface parameters such as collocation points, surface vectors, and allocate necessary matrices.
             */
            void get_surface_parameters()
            {
                // Get collocation points
                UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
                UVLM::Geometry::generate_colocationMesh(u_ext, uext_col);
                Ktotal = UVLM::Matrix::get_total_VecVecMat_size(uext_col);
                // Alocate
                UVLM::Types::allocate_VecVecMat(uext_total, u_ext);
                UVLM::Types::copy_VecVecMat(u_ext, uext_total);
                UVLM::Types::allocate_VecVecMat(uext_total_col, u_ext, -1);                
                // surface vectors                
                UVLM::Types::allocate_VecVecMat(normals, zeta_col);
                UVLM::Types::allocate_VecVecMat(longitudinals, zeta_col);
                UVLM::Types::allocate_VecVecMat(perpendiculars, zeta_col);
                UVLM::Geometry::generate_surface_vectors(zeta, normals, longitudinals, perpendiculars); 
            }
        };
        /**
         * @struct lifting_surface
         * @brief Represents a lifting aerodynamic surface (i.e. wings).
        */
        struct lifting_surface : surface
        {

            UVLM::Types::VecDimensions dimensions_star; // An array of dimensions (chordwise and spanwise) for wake surfaces.
            UVLM::Types::VecMapX gamma; // A mapped matrix of circulation strength for each vortex ring panel of each surface
            UVLM::Types::VecMapX gamma_star; // A mapped matrix of circulation strength for each wake vortex ring panel of each surface.
            UVLM::Types::VecVecMapX zeta_dot; // A mapped matrix containing corner point velocities for each surface grid.
            UVLM::Types::VecVecMapX zeta_star; // A mapped matrix containing corner point coordinates for each surface grid.
            UVLM::Types::MatrixX aic; // Matrix with aerodynamic influence coefficients considering all surfaces.
            UVLM::Types::VectorX rbm_vel_g; // An array of rigid body motion velocities.
            UVLM::Types::VectorX centre_rot; // An array of rotational velocities of rigid bodies around their centers.
            UVLM::Types::VecVecMatrixX coordinates_center_spanwise_vertices; // A matrix with coordinates of the center of each spanwise vortex.
            UVLM::Types::VecVecMatrixX coordinates_center_chordwise_vertices; // A matrix with coordinates of the center of each chordwise vortex.
            UVLM::Types::VecVecMatrixX u_induced_col_sources; // A matrix with the induced velocities by sources on the lifting surfaces collocation points.
            UVLM::Types::VecVecMatrixX u_induced_by_sources_on_center_chordwise_vertices; // A matrix with the induced velocities by sources on the center of each chordwise vortex.
            UVLM::Types::VecVecMatrixX u_induced_by_sources_on_center_spanwise_vertices; // A matrix with the induced velocities by sources on the center of each spanwise vortex.
            
            /**
             * @brief Constructor for a lifting aerodynamic surface.
             * 
             * This constructor takes input pointers from SHARPy and map them to Eigen matrices, 
             * and allocates space for other matrices used later. 
             *
             * @param n_surfaces Number of surfaces.
             * @param p_dimensions Pointer to matrix of surface dimensions.
             * @param p_zeta Pointer to matrix with corner point coordinates of each surface grid.
             * @param p_u_ext Pointer to matrix with external velocities at each surface corner point.
             * @param p_forces Pointer to matrix with forces at each surface corner point.
             * @param p_zeta_star Pointer to matrix with corner point coordinates of each wake surface grid.
             * @param p_zeta_dot Pointer to matrix with corner point velocities of each surface grid.
             * @param p_gamma Pointer to matrix with circulation strength for each vortex ring panel on each surface.
             * @param p_gamma_star Pointer to matrix with circulation strength for each wake vortex ring panel on each wake surface.
             * @param p_dimensions_star Pointer to matrix of wake surface dimensions.
             * @param p_rbm_vel Pointer to array with rigid body motion velocities.
             * @param p_centre_rot Pointer to array with rotational velocities of rigid bodies around their centers.
             */
            lifting_surface
            (
                uint n_surfaces,
                unsigned int** p_dimensions,
                double** p_zeta,
                double** p_u_ext,
                double** p_forces,
                double** p_zeta_star,
                double** p_zeta_dot,
                double** p_gamma,
                double** p_gamma_star,
                unsigned int** p_dimensions_star,
                double*  p_rbm_vel,
                double*  p_centre_rot
            ):surface{n_surfaces, p_dimensions, p_zeta, p_u_ext, p_forces}
            {    
                UVLM::Mapping::transform_dimensions(n_surf,p_dimensions_star,dimensions_star);                
                UVLM::Mapping::map_VecVecMat(dimensions_star,p_zeta_star,zeta_star,1);
                UVLM::Mapping::map_VecVecMat(dimensions, p_zeta_dot, zeta_dot, 1);
                UVLM::Mapping::map_VecMat(dimensions, p_gamma, gamma, 0);
                UVLM::Mapping::map_VecMat(dimensions_star, p_gamma_star, gamma_star, 0);
            
                UVLM::Types::allocate_VecVecMat(u_induced_col_sources, u_ext, -1);  
                UVLM::Types::allocate_VecVecMat(coordinates_center_spanwise_vertices, u_ext, -1);
                UVLM::Types::allocate_VecVecMat(coordinates_center_chordwise_vertices, coordinates_center_spanwise_vertices); 
                UVLM::Types::allocate_VecVecMat(u_induced_by_sources_on_center_chordwise_vertices, coordinates_center_spanwise_vertices); 
                UVLM::Types::allocate_VecVecMat(u_induced_by_sources_on_center_spanwise_vertices, coordinates_center_spanwise_vertices);  
                
                UVLM::Mapping::map_VecX(2*UVLM::Constants::NDIM, p_rbm_vel, rbm_vel_g);   
                UVLM::Mapping::map_VecX(2*UVLM::Constants::NDIM, p_centre_rot, centre_rot);   
            }

            /**
             * @brief Get surface parameters such as collocation points, surface vectors, and allocate necessary matrices.
             */
            void get_surface_parameters()
            {
                surface::get_surface_parameters();
            }
            /**
             * @brief Compute boundary conditions and AIC matrices..
             *
             * @param options UVLM options set in SHARPy.
             */
            void get_aerodynamic_solver_inputs(const UVLM::Types::VMopts& options)
            {
                rhs.resize(Ktotal); 
                UVLM::Matrix::RHS(zeta_col,
                                zeta_star,
                                uext_total_col,
                                gamma_star,
                                normals,
                                options,
                                rhs,
                                Ktotal);
                aic = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
                UVLM::Matrix::AIC(zeta,
                                zeta_col,
                                zeta_star,
                                uext_total_col,
                                normals,
                                options,
                                options.horseshoe,
                                aic);      
            }
            /**
             * @brief Get induced velocities at collocation points from sources.
             * 
             * Used later in the force calculation.
             *
             * @param sigma_flat Flattened source strengths.
             * @param aic_nonlifting_on_lifting_x AIC matrix for sources (nonlifting bodies) on vortex ring panels ( lifting surfaces) in the x-direction.
             * @param aic_nonlifting_on_lifting_y AIC matrix for sources (nonlifting bodies) on vortex ring panels ( lifting surfaces) in the y-direction.
             * @param aic_nonlifting_on_lifting_z AIC matrix for sources (nonlifting bodies) on vortex ring panels ( lifting surfaces) in the z-direction.
             * @param Ktotal_nonlifting Total number of collocation points located on non-lifting surfaces.
             */
            void get_induced_col_from_sources(UVLM::Types::VectorX& sigma_flat,
                                              UVLM::Types::MatrixX& aic_nonlifting_on_lifting_x,
                                              UVLM::Types::MatrixX& aic_nonlifting_on_lifting_y,
                                              UVLM::Types::MatrixX& aic_nonlifting_on_lifting_z,
                                              const uint Ktotal_nonlifting)
            {
                 
                    UVLM::Types::allocate_VecVecMat(u_induced_col_sources, uext_col);
                    UVLM::PostProc::calculate_induced_velocity_col(sigma_flat,
                                                            aic_nonlifting_on_lifting_x,
                                                            aic_nonlifting_on_lifting_y,
                                                            aic_nonlifting_on_lifting_z,
                                                            u_induced_col_sources);

                    UVLM::PostProc::transform_VecVecMatrix_to_global_coordinate_system
                    (
                    normals,
                    longitudinals,
                    perpendiculars,
                    u_induced_col_sources
                    );
            }
            /**
             * @brief Calculate the coordinates of the center vertices of the lifting surface.
             */
            void get_coordinates_center_vertices()
            {   
                // TODO: parallelization
                for (uint isurf=0; isurf < n_surf; isurf++)
                {
                    for (uint icol=0; icol < dimensions[isurf].second; icol++)
                    {
                        for (uint irow=0; irow < dimensions[isurf].first; irow++)
                        {                            
                            for (uint idim=0; idim < UVLM::Constants::NDIM; idim++)
                            {  

                                coordinates_center_chordwise_vertices[isurf][idim](irow, icol) = 0.5 * (zeta[isurf][idim](irow, icol) +zeta[isurf][idim](irow + 1, icol));
                                coordinates_center_spanwise_vertices[isurf][idim](irow, icol) = 0.5 * (zeta[isurf][idim](irow, icol) +zeta[isurf][idim](irow, icol + 1));
                            }
                        }

                    }

                }
            }


        };
        /**
         * @struct lifting_surface_unsteady
         * @brief Represents a lifting aerodynamic surface considering unsteady aerodynamic effects.
        */
        struct lifting_surface_unsteady : lifting_surface
        {
            UVLM::Types::VecVecMapX dynamic_forces; // A mapped matrix with unsteady forces at each surface corner point.
            UVLM::Types::VecVecMapX uext_star; // A mapped matrix with external velocities at each wake surface corner point.
            UVLM::Types::VecMapX dist_to_orig; // A mapped matrix containing of distances from the trailing edge of the wake vertices.
            UVLM::Types::VecVecMatrixX solid_vel; // A matrix with resulting grid velocities (no flow velocities) for each corner point of the lifting surfaces.
            UVLM::Types::VecVecMatrixX uext_star_total; // A matrix with total velocities acting on each wake surface corner point.
            UVLM::Types::VecVecMapX normals_unsteady; // A mapped matrix with the normal vectors of the panels of each surface.
                
            UVLM::Types::VecMatrixX extra_gamma_star; // A matrix to store last wake panel circulation information.
            UVLM::Types::VecVecMatrixX extra_zeta_star;  // A matrix to store last wake panel geometry information.
            
            /**
             * @brief Constructor for an unsteady lifting aerodynamic surface.
             *
             * @param n_surfaces Number of surfaces.
             * @param p_dimensions Pointer to matrix of surface dimensions.
             * @param p_zeta Pointer to matrix with corner point coordinates of each surface grid.
             * @param p_u_ext Pointer to matrix with external velocities at each surface corner point.
             * @param p_forces Pointer to matrix with forces at each surface corner point.
             * @param p_zeta_star Pointer to matrix with corner point coordinates of each wake surface grid.
             * @param p_zeta_dot Pointer to matrix with corner point velocities of each surface grid.
             * @param p_gamma Pointer to matrix with circulation strength for each vortex ring panel on each surface.
             * @param p_gamma_star Pointer to matrix with circulation strength for each wake vortex ring panel on each wake surface.
             * @param p_dimensions_star Pointer to matrix of wake surface dimensions.
             * @param p_rbm_vel Pointer to array with rigid body motion velocities.
             * @param p_centre_rot Pointer to array with rotational velocities of rigid bodies around their centers.
             * @param p_dist_to_orig Pointer to array of distances from the trailing edge of the wake vertices.
             * @param p_dynamic_forces Pointer to matrix with unsteady forces at each surface corner point.
             * @param p_uext_star Pointer to matrix with external velocities at each wake surface corner point.
             * @param p_normals_unsteady Pointer to matrix with the normal vectors of the panels of each surface.
             */
            lifting_surface_unsteady
            (
                uint n_surfaces,
                unsigned int** p_dimensions,
                double** p_zeta,
                double** p_u_ext,
                double** p_forces,
                double** p_zeta_star,
                double** p_zeta_dot,
                double** p_gamma,
                double** p_gamma_star,
                unsigned int** p_dimensions_star,
                double*  p_rbm_vel,
                double*  p_centre_rot,
                double**p_dist_to_orig,
                double**p_dynamic_forces,
                double** p_uext_star,
                double** p_normals_unsteady
            ):lifting_surface{n_surfaces, p_dimensions, p_zeta, p_u_ext, p_forces, p_zeta_star, p_zeta_dot, p_gamma, p_gamma_star, p_dimensions_star, p_rbm_vel, p_centre_rot}
            {    
                
                UVLM::Mapping::map_VecVecMat(dimensions,
                                            p_dynamic_forces,
                                            dynamic_forces,
                                            1,
                                            2*UVLM::Constants::NDIM);
                UVLM::Mapping::map_VecMat(dimensions_star,
                                            p_dist_to_orig,
                                            dist_to_orig,
                                            1);  
                UVLM::Mapping::map_VecVecMat(dimensions_star,
                                            p_uext_star,
                                            uext_star,
                                            1);   
                                            
                UVLM::Mapping::map_VecVecMat(dimensions,
                                                p_normals_unsteady,
                                                normals_unsteady,
                                                0);         
            }
            /**
             * @brief Get surface parameters such as collocation points, surface vectors, and allocate necessary matrices.
             */
            void get_surface_parameters()
            {
                lifting_surface::get_surface_parameters();
                UVLM::Types::allocate_VecVecMat(solid_vel, u_ext);
                UVLM::Geometry::generate_surfaceNormal(zeta, normals_unsteady);
            }
        };

        /**
         * @struct phantom_surface
         * @brief Represents a phantom surface used for fuselage-wing junction handling.
         */
        struct phantom_surface
        {

            bool phantom_cell_required; // Flag if a phantom cell is required for a specific surface.
            uint Ktotal; // Number of total collocation points.
            uint n_surf; // Number of surfaces.
            UVLM::Types::VecVecMatrixX zeta; // A matrix containing corner point coordinates for each surface grid.
            UVLM::Types::VecVecMatrixX zeta_col; // A  matix containing collocation point coordinates located on each surface grid.
            UVLM::Types::VecVecMatrixX zeta_star; // A matrix containing corner point coordinates for each wake surface grid.
            UVLM::Types::VecVecMatrixX normals; // A  matrix containing the normal vectors for the panels of each surface.
            UVLM::Types::VecVecMatrixX longitudinals; // A matrix containing the longitudinal vectors of the panels of each surface.
            UVLM::Types::VecVecMatrixX perpendiculars; // A matrix containing the perpendicular vectors for the panels of each surface. 
            UVLM::Types::VecMatrixX gamma; // A matrix of circulation strength for each vortex ring panel of each surface.
            UVLM::Types::VecMatrixX gamma_star; // A matrix of circulation strength for each wake vortex ring panel of each surface.
            UVLM::Types::VectorX gamma_flat; // A vector that contains the circulation strength value stored in gamma reshaped into a vector.
            UVLM::Types::MatrixXint flag_zeta_phantom; // A vector containing the number of the junction partner surface if existing.
                
            /**
             * @brief Constructor for a phantom surface.
             *
             * @param p_flag_zeta_phantom Pointer to  vector with the number of the junction partner surface if existing.
             * @param n_surfaces Number of surfaces.
             * @param zeta_lifting Mapped matrix with corner point coordinates of the lifting surface grid.
             * @param zeta_lifting_star Mapped matrix with corner point coordinates of the lifting wake surface grid.
             * @param dimensions_lifting Matrix of dimensions of the lifting surfaces surface dimensions.
             */
            phantom_surface
            (
                int* p_flag_zeta_phantom,
                uint n_surfaces,
                UVLM::Types::VecVecMapX zeta_lifting,
                UVLM::Types::VecVecMapX zeta_lifting_star,
                UVLM::Types::VecDimensions dimensions_lifting
            )
            {   
                n_surf = n_surfaces;
                // TODO: Change flag_zeta_phantom to VecMat type
                UVLM::Types::MapMatrixXint copy_flag_zeta_phantom(p_flag_zeta_phantom, 1, n_surf);
                uint N_row = copy_flag_zeta_phantom.rows();
                uint N_col = copy_flag_zeta_phantom.cols();
                flag_zeta_phantom.resize(N_row, N_col);
                for (uint i_row=0;i_row<N_row;i_row++)
                {
                    for (uint i_col=0;i_col<N_col;i_col++)
                    {
                        flag_zeta_phantom(i_row, i_col) = copy_flag_zeta_phantom(i_row, i_col);
                    }
                }
                phantom_cell_required = UVLM::Phantom::check_for_true_in_bool_vec_mat(flag_zeta_phantom);
                if (phantom_cell_required)
                {
                    UVLM::Phantom::create_phantom_zeta(zeta_lifting, zeta, flag_zeta_phantom);                     
                    UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
                    Ktotal = UVLM::Matrix::get_total_VecVecMat_size(zeta_col);                    
                    // Allocate phantom gamma                    
                    UVLM::Types::allocate_VecMat_from_VecVecMat(gamma, zeta, -1);
                    gamma_flat.resize(Ktotal);
                    // Allocate phantom gamma star 
                    // Allocate with all surfaces (Will cause error in Wing-tail combinations)                   
                    UVLM::Types::allocate_VecMat(gamma_star,
                                                 n_surf,
                                                 zeta_lifting_star[0][0].rows()-1,
                                                 zeta[0][0].cols()-1);
                    update_wake(zeta_lifting_star);   
                }
                else
                {
                    Ktotal = 0;
                }                
            }

            /**
             * @brief Get surface parameters such as surface vectors, and allocate necessary matrices.
             */
            void get_surface_parameters()
            {
                UVLM::Types::allocate_VecVecMat(normals, zeta, -1);   
                UVLM::Types::allocate_VecVecMat(longitudinals , zeta,-1 );
                UVLM::Types::allocate_VecVecMat(perpendiculars, zeta, -1); 
                UVLM::Geometry::generate_surface_vectors(zeta, normals, longitudinals, perpendiculars);                
            }

            /**
             * @brief Update the phantom wake based on the lifting wake lattice grid.
             *
             * @param zeta_lifting_star Mapped matrix with corner point coordinates of the lifting wake surface grid.
             */
            void update_wake(UVLM::Types::VecVecMapX zeta_lifting_star)
            {
                UVLM::Phantom::create_phantom_zeta_star(flag_zeta_phantom, zeta, zeta_lifting_star, zeta_star);   
            }

            /**
             * @brief Update gamma of phantom panels using interpolation of the circulation strength of the lifting junction panels.
             *
             * @param Ktotal_lifting Number of total collocation points of the lifting surfaces.
             * @param zeta_col_lifting  A matix containing collocation point coordinates of the lifting surfaces.
             * @param gamma_lifting A mapped matrix of circulation strength for each vortex ring panel of the lifting surfaces.
             */
            void update_gamma(const uint Ktotal_lifting,
                              UVLM::Types::VecVecMatrixX& zeta_col_lifting,
                              UVLM::Types::VecMapX& gamma_lifting)
            {
                // Allocate aic
                UVLM::Types::MatrixX aic = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal_lifting);
                // get aic phantom interp condition but without -1
                UVLM::Matrix::aic_phantom_interp_condition
                (
                    Ktotal_lifting,
                    Ktotal, //phantom
                    zeta_col_lifting,
                    zeta_col, //phantom
                    aic,
                    flag_zeta_phantom,
                    true
                );
                // get gamma phantom gamma = aic*gamma_lifting;
                gamma_flat.setZero();
                UVLM::Types::VectorX gamma_lifting_flat;
                UVLM::Matrix::reconstruct_vector_values_from_VecMatrixX(gamma_lifting,
                                                gamma_lifting_flat,
                                                zeta_col_lifting);
                                                
           
                for (uint i_row = 0; i_row < Ktotal; ++i_row)
                {
                    for (uint i_col = 0; i_col < Ktotal_lifting; ++i_col)
                    {
                        gamma_flat(i_row) += aic(i_row, i_col)* gamma_lifting_flat(i_col);
                    }
                    
                }
                UVLM::Matrix::reconstruct_VecMatrixX_values_from_vector(gamma_flat,
                                                                    gamma,
                                                                    zeta_col);
            }

            /**
             * @brief Update gamma_star based on interpolation of lifting surface gamma_star.
             *
             * @param zeta_lifting_star Mapped matrix with corner point coordinates of the lifting wake surface grid.
             * @param gamma_lifting A mapped matrix of circulation strength for each vortex ring panel of the lifting wake surfaces.
             */
            void update_gamma_wake(UVLM::Types::VecVecMapX zeta_lifting_star,
                                   UVLM::Types::VecMapX gamma_lifting_star)
            {
                double factor = 0.0;
                uint N_spanwise_panels, N_streamwise_panels;
                UVLM::Types::VecVecMatrixX zeta_col_lifting_star;
                UVLM::Geometry::generate_colocationMesh(zeta_lifting_star,zeta_col_lifting_star);
                for(uint i_surf=0; i_surf<n_surf; ++i_surf)
                {                    
                    N_spanwise_panels = gamma_star[i_surf].cols();
                    N_streamwise_panels = gamma_star[i_surf].rows();
                          
                    for (uint i_row=0; i_row<N_streamwise_panels; ++i_row)
                    {                
                        factor = (abs(gamma_lifting_star[1](i_row, 0)) - abs(gamma_lifting_star[0](i_row, 0)))
                                 /(zeta_lifting_star[1][1](i_row, 0) - zeta_lifting_star[0][1](i_row, 0));
                        for (uint i_col=0; i_col<N_spanwise_panels; ++i_col) 
                        {
                            gamma_star[i_surf](i_row, i_col) = abs(gamma_lifting_star[0](i_row, 0))
                                                                + (zeta_star[i_surf][1](i_row, i_col) - zeta_lifting_star[0][1](i_row, 0))
                                                                * factor;
                            if (gamma_lifting_star[i_surf](i_row, 0) < 0.0)
                            {
                                gamma_star[i_surf](i_row, i_col) *= -1;
                            }
                        
                        }              
                     }
                }
            }
        };
        
    /**
     * @struct nonlifting_body
     * @brief Structure representing the surface of the non-lifting bodies in the UVLM solver.
     */
    struct nonlifting_body : surface
    {      
    // Nonlifting body-specific attributes
    UVLM::Types::VecMapX sigma;                // Vector of sigma values.
    UVLM::Types::VecMapX pressure_coefficients; // Vector of pressure coefficients.
    UVLM::Types::VecVecMatrixX u_induced_col;   // Matrix of induced velocities at collocation points.

    // Aerodynamic Solver Inputs
    UVLM::Types::MatrixX aic_sources_x; // Matrix for AIC (Aerodynamic Influence Coefficient) in the x-direction.
    UVLM::Types::MatrixX aic_sources_y; // Matrix for AIC (Aerodynamic Influence Coefficient) in the y-direction.
    UVLM::Types::MatrixX aic_sources_z; // Matrix for AIC (Aerodynamic Influence Coefficient) in the z-direction.
    UVLM::Types::VectorX sigma_flat;     // Flattened sigma values.

    /**
     * @brief Constructor for the nonlifting_body struct.
     *
     * @param n_surfaces Number of non-lifting bodies.
     * @param p_dimensions Pointer to an array of dimensions (chordwise and spanwise) for each non-lifting body.
     * @param p_zeta Pointer to an array of coordinates for each non-lifting body.
     * @param p_u_ext Pointer to an array of external velocities for each non-lifting body.
     * @param p_forces Pointer to an array of forces for each non-lifting body.
     * @param p_pressure_coefficient_nonlifting Pointer to an array of pressure coefficients for each non-lifting body.
     * @param p_sigma Pointer to an array of sigma values for each non-lifting body.
     */
    nonlifting_body(
        uint n_surfaces,
        unsigned int** p_dimensions,
        double** p_zeta,
        double** p_u_ext,
        double** p_forces,
        double** p_pressure_coefficient_nonlifting,
        double** p_sigma
    ) : surface{n_surfaces, p_dimensions, p_zeta, p_u_ext, p_forces}
    {   
        UVLM::Mapping::map_VecMat(dimensions, p_sigma, sigma, 0);
        UVLM::Mapping::map_VecMat(dimensions, p_pressure_coefficient_nonlifting, pressure_coefficients, 0); 
    }

    /**
     * @brief Initialize surface parameters and allocate memory for induced velocities.
     *
     * @param phantom_wing_test Flag indicating whether this is a phantom wing test.
     */
    void get_surface_parameters(bool phantom_wing_test = false)
    {
        surface::get_surface_parameters();
        
        if (phantom_wing_test)
        {
            Ktotal = 0;
        }
        else
        {
            Ktotal = UVLM::Matrix::get_total_VecVecMat_size(uext_col);
        }

        UVLM::Types::allocate_VecVecMat(u_induced_col, uext_col);
    }

    /**
     * @brief Calculate and set the aerodynamic solver inputs including AIC and RHS.
     *
     * @param phantom_wing_test Flag indicating whether this is a phantom wing test.
     */
    void get_aerodynamic_solver_inputs(bool phantom_wing_test = false)
    {
        aic_sources_x = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
        aic_sources_y = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
        aic_sources_z = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
        rhs.resize(Ktotal);
        if (!phantom_wing_test)
        {
            UVLM::Matrix::RHS_nonlifting_body(uext_col,
                        normals,
                        rhs,
                    Ktotal,
                    n_surf);
        UVLM::Matrix::AIC_sources(zeta,
                                zeta_col,
                                longitudinals,
                                perpendiculars,
                                normals,
                                longitudinals,
                                perpendiculars,
                                normals,
                                aic_sources_x,
                                aic_sources_y,
                                aic_sources_z);
        }                
    }

    };
    }
}