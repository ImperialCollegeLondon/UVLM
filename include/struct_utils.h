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
        struct surface
        {
            // Parameters
            uint n_surf, Ktotal;
            UVLM::Types::VecDimensions dimensions;
            UVLM::Types::VecVecMapX zeta, u_ext, forces; 
                   
            UVLM::Types::VecVecMatrixX zeta_col, uext_col, uext_total, uext_total_col;
            UVLM::Types::VecVecMatrixX normals, longitudinals, perpendiculars;
            UVLM::Types::VectorX rhs;
            // Constructor
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

            //Functions 
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

        struct lifting_surface : surface
        {
            UVLM::Types::VecDimensions dimensions_star;
            UVLM::Types::VecMapX gamma, gamma_star;
            UVLM::Types::VecVecMapX zeta_dot, zeta_star;
            UVLM::Types::MatrixX aic;
            UVLM::Types::VectorX rbm_vel_g;
            
            // Constructor
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
                double* p_rbm_vel_g,
                unsigned int** p_dimensions_star
            ):surface{n_surfaces, p_dimensions, p_zeta, p_u_ext, p_forces}
            {    
                UVLM::Mapping::transform_dimensions(n_surf,p_dimensions_star,dimensions_star);                
                UVLM::Mapping::map_VecVecMat(dimensions_star,p_zeta_star,zeta_star,1);
                UVLM::Mapping::map_VecVecMat(dimensions, p_zeta_dot, zeta_dot, 1);
                UVLM::Mapping::map_VecMat(dimensions, p_gamma, gamma, 0);
                UVLM::Mapping::map_VecMat(dimensions_star, p_gamma_star, gamma_star, 0);
                UVLM::Mapping::map_VecX(2*UVLM::Constants::NDIM, p_rbm_vel_g, rbm_vel_g);
            }

            //Functions
            void get_surface_parameters()
            {
                surface::get_surface_parameters();
            }
            
            void get_aerodynamic_solver_inputs(const UVLM::Types::VMopts& options)
            {
                rhs.resize(Ktotal); 
                UVLM::Matrix::RHS(zeta_col,
                                zeta_star,
                                uext_col,
                                gamma_star,
                                normals,
                                options,
                                rhs,
                                Ktotal);
                aic = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
                UVLM::Matrix::AIC(zeta,
                                zeta_col,
                                zeta_star,
                                uext_col,
                                normals,
                                options,
                                false,
                                aic);          }
        };

        struct lifting_surface_unsteady : lifting_surface
        {
            UVLM::Types::VecVecMapX dynamic_forces, uext_star;
            UVLM::Types::VecMapX dist_to_orig;  
            UVLM::Types::VecVecMatrixX solid_vel; //
            // Constructor
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
                double* p_rbm_vel_g,
                unsigned int** p_dimensions_star,
                double**p_dist_to_orig,
                double**p_dynamic_forces,
                double** p_uext_star
            ):lifting_surface{n_surfaces, p_dimensions, p_zeta, p_u_ext, p_forces, p_zeta_star, p_zeta_dot, p_gamma, p_gamma_star, p_rbm_vel_g, p_dimensions_star}
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
            }
            void get_surface_parameters()
            {
                lifting_surface::get_surface_parameters();
                UVLM::Types::allocate_VecVecMat(solid_vel, u_ext); //
            }
        };
        struct phantom_surface
        {
            bool phantom_cell_required;
            uint Ktotal, n_surf;
            UVLM::Types::VecVecMatrixX zeta, zeta_col, zeta_star, normals, longitudinals, perpendiculars;

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
                UVLM::Types::MapMatrixXint flag_zeta_phantom(p_flag_zeta_phantom, dimensions_lifting[0].second, n_surf);
                phantom_cell_required = UVLM::Phantom::check_for_true_in_bool_vec_mat(flag_zeta_phantom);
                if (phantom_cell_required)
                {
                    UVLM::Phantom::create_phantom_zeta(zeta_lifting, zeta, flag_zeta_phantom);                     

                }
                else
                {
                    // TODO: Check if code adaptions needed for this case
                    Ktotal = 0;
                }
                
            }

            void get_surface_parameters()
            {
                UVLM::Geometry::generate_colocationMesh(zeta, zeta_col);
                Ktotal = UVLM::Matrix::get_total_VecVecMat_size(zeta_col);
                UVLM::Types::allocate_VecVecMat(normals, zeta, -1);   
                UVLM::Types::allocate_VecVecMat(longitudinals , zeta,-1 );
                UVLM::Types::allocate_VecVecMat(perpendiculars, zeta, -1); 
                UVLM::Geometry::generate_surface_vectors(zeta, normals, longitudinals, perpendiculars);                
            }

            void update_wake(UVLM::Types::VecVecMapX zeta_lifting_star)
            {
                UVLM::Phantom::create_phantom_zeta_star(zeta, zeta_lifting_star, zeta_star);   
            }
        };

        struct nonlifting_body : surface
        {
            // Nonlifting body specifix
            UVLM::Types::VecMapX sigma;
            UVLM::Types::VecVecMatrixX u_induced_col;

            // Aerodynamic Solver Inputs
            UVLM::Types::MatrixX aic_sources_x, aic_sources_y,aic_sources_z;            
            UVLM::Types::VectorX sigma_flat;
            nonlifting_body
            (
                uint n_surfaces,
                unsigned int** p_dimensions,
                double** p_zeta,
                double** p_u_ext,
                double** p_forces,
                double** p_sigma
            ):surface{n_surfaces, p_dimensions, p_zeta, p_u_ext, p_forces}
            {   
                UVLM::Mapping::map_VecMat(dimensions, p_sigma, sigma, 0);
            }
            void get_surface_parameters()
            {
                surface::get_surface_parameters();
                UVLM::Types::allocate_VecVecMat(u_induced_col, uext_col);
            }
            void get_aerodynamic_solver_inputs()
            {
                rhs.resize(Ktotal); 
                UVLM::Matrix::RHS_nonlifting_body(uext_col,
                            normals,
                            rhs,
                            Ktotal,
                            n_surf);
                            
                aic_sources_x = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
                aic_sources_y = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
                aic_sources_z = UVLM::Types::MatrixX::Zero(Ktotal, Ktotal);
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
        };
    }
}