/**
 * @file sources.h
 * @brief Header file containing functions related to source panel calculations in the UVLM framework.
 */

#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "debugutils.h"

#include <cmath>

/**
 * @namespace UVLM
 * @brief Namespace for the UVLM (Unsteady Vortex Lattice Method) framework.
 */
namespace UVLM
{
    /**
     * @namespace Sources
     * @brief Namespace for functions related to source panel calculations.
     */
    namespace Sources
    {
        /**
         * @brief Calculate influence coefficients for panels on the aerodynamic surface.
         *
         * This function computes the aerodynamic influence coefficients (AIC) for each panel on the aerodynamic nonlifting 
		 * surface (e.g. fuselage). These AICs represent the effect of velocity induced by a quadrilateral linear source panel 
		 * on an arbitary point (e.g. collocation points of a surface panel).
		 * 
		 * The equations used for the computations can be found Chapter 10.4.1 of:
		 * Katz, J., and Plotkin, A., Low-speed aerodynamics, 2nd, Cambridge University Press, 2001.
         *
         * @param zeta The coordinates of the aerodynamic nonlifting surface.
         * @param target_surface The collocation points of the target surface to compute influence coefficients for.
         * @param aic_x Output array for influence coefficients in the X-direction.
         * @param aic_y Output array for influence coefficients in the Y-direction.
         * @param aic_z Output array for influence coefficients in the Z-direction.
         * @param longitudinal_panel Longitudinal surface vectors of the panels in zeta.
         * @param perpendicular_panel Perpendicular surface vectors of the panels in zeta.
         * @param normal_panel Normal surface vectors of the panels in zeta.
         * @param longitudinal_col Longitudinal vectors of surface associated with the collocations included in target surface.
         * @param perpendicular_col Perpendicular vectors of surface associated with the collocations included in target surface.
         * @param normal_col Normal vectors of surface associated with the collocations included in target surface.
         * @param same_surface Flag indicating whether the target and source surfaces are the same.
         */
        template <typename t_zeta,
                  typename t_tsurface,
                  typename t_aic,
                  typename t_surf_vec_panel,
                  typename t_surf_vec_col>
        void get_influence_coefficient(
            const t_zeta& zeta,
            const t_tsurface& target_surface,
            t_aic& aic_x,
            t_aic& aic_y,
            t_aic& aic_z,
            const t_surf_vec_panel& longitudinal_panel = NULL,
            const t_surf_vec_panel& perpendicular_panel = NULL,
            const t_surf_vec_panel& normal_panel = NULL,
            const t_surf_vec_col& longitudinal_col = NULL,
            const t_surf_vec_col& perpendicular_col = NULL,
            const t_surf_vec_col& normal_col = NULL,
            const bool& same_surface = false
        );

        /**
         * @brief Calculate the Q-vector for induced velocity calculations.
         *
         * This function computes the Q-vector used in the calculation of induced velocities with
		 * \f$Q_{ij} = \ln{\frac{r_i + r_j - d_{ij}}{r_i + r_j + d_{ij}}}\f$
		 * found in Eq. 10.95 and 10.96 (Katz & Plotzkin, 2001).
         *
         * @param r Vector of radii for each panel.
         * @param d Vector of distances between panels.
         * @param Q Output array for the Q-vector.
         */
        void get_Q(
            UVLM::Types::VectorX& r,
            UVLM::Types::VectorX& d,
            UVLM::Types::VectorX& Q
        );

        /**
         * @brief Calculate the summation used in the computation of the vertical induced velocity..
         *
         * This function calculates the summation used in the computation of vertical induced velocit
         * for cases where the collocation points is not associated with the present source  panel. 
		 * Please refer to Eq. 10.97 (Katz & Plotzkin, 2001).
		 * 
         * @param r Vector of r Eq. 10.92 (Katz & Plotzkin, 2001).
         * @param diff_epsilon Vector representing differences in panel epsilon coordinates.
         * @param diff_eta Vector representing differences in panel eta coordinates.
         * @param delta_epsilon_x_vec Delta epsilon/x coordinates of panel with collocation point.
         * @param delta_eta_y_vec Delta eta/y coordinates of panel with collocation point.
         * @param z Z-coordinate for collocation points.
         * @return The computed summation value.
         */
        double get_j_sum(
            const UVLM::Types::VectorX& r,
            const UVLM::Types::VectorX& diff_epsilon,
            const UVLM::Types::VectorX& diff_eta,
            const UVLM::Types::VectorX& delta_epsilon_x_vec,
            const UVLM::Types::VectorX& delta_eta_y_vec,
            const double z
        );

        /**
         * @brief Calculate the dot product between two vectors.
         *
         * This function computes the dot product between two input vectors.
         *
         * @tparam vector Type of the input vectors.
         * @param vec1 The first vector.
         * @param vec2 The second vector.
         * @return The dot product between the two vectors.
         */
        template <typename vector>
        double dot_product
        (
             const vector& vec1,
             const vector& vec2
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
    }
}


template <typename t_zeta,
		  typename t_tsurface,
		  typename t_aic,
		  typename t_surf_vec_panel,
		  typename t_surf_vec_col>

void UVLM::Sources::get_influence_coefficient
    (
    const t_zeta&       zeta,
    const t_tsurface&   target_surface,
    t_aic&             aic_x,
    t_aic&             aic_y,
    t_aic&             aic_z,
	const t_surf_vec_panel&    longitudinal_panel,
	const t_surf_vec_panel&    perpendicular_panel,
	const t_surf_vec_panel&    normal_panel,
	const t_surf_vec_col&    longitudinal_col,
	const t_surf_vec_col&    perpendicular_col,
	const t_surf_vec_col&    normal_col,
	const bool& same_surface
    )
    {

    const unsigned int rows_collocation = target_surface[0].rows();
    const unsigned int cols_collocation = target_surface[0].cols();
    const unsigned int rows_panel = zeta[0].rows()-1;
    const unsigned int cols_panel = zeta[0].cols()-1;
	//keep
    UVLM::Types::Vector3 normal_panel_vec;
    UVLM::Types::Vector3 longitudinal_panel_vec;
    UVLM::Types::Vector3 perpendicular_panel_vec;

	//new
	UVLM::Types::VectorX panel_coordinates_epsilon;
	UVLM::Types::VectorX panel_coordinates_eta;
	UVLM::Types::VectorX panel_coordinates_z;
		
    UVLM::Types::VectorX delta_epsilon;
    UVLM::Types::VectorX delta_eta_vec;
	
	UVLM::Types::Vector3 collocation_point_transf;
	UVLM::Types::Vector3 induced_velocity_vec;
	UVLM::Types::VectorX delta_epsilon_x_vec;
	UVLM::Types::VectorX delta_eta_y_vec;
	UVLM::Types::VectorX r;
	UVLM::Types::VectorX d;
	UVLM::Types::VectorX S_vec;
	UVLM::Types::VectorX C_vec;
	UVLM::Types::VectorX Q;
	UVLM::Types::VectorX J_vec;

    UVLM::Types::Vector3 normal_col_vec;
    UVLM::Types::Vector3 longitudinal_col_vec;
    UVLM::Types::Vector3 perpendicular_col_vec;

	bool flag_triangle = false;
	bool flag_if_col_on_panel = false;
	int ignore_index;
    //PANELS
    uint panel_id = 0;
    uint collocation_id = 0;
    for (unsigned int i_panel=0; i_panel<rows_panel; ++i_panel)
    {
        for (unsigned int j_panel=0; j_panel<cols_panel; ++j_panel)
        {
            collocation_id = 0;
			// Convert Source Panel Corner Points from Global to Panel frame of reference
            longitudinal_panel_vec = UVLM::Types::Vector3(longitudinal_panel[0](i_panel, j_panel), longitudinal_panel[1](i_panel, j_panel), longitudinal_panel[2](i_panel, j_panel));
            perpendicular_panel_vec = UVLM::Types::Vector3(perpendicular_panel[0](i_panel, j_panel), perpendicular_panel[1](i_panel, j_panel), perpendicular_panel[2](i_panel, j_panel));
            normal_panel_vec = UVLM::Types::Vector3(normal_panel[0](i_panel, j_panel), normal_panel[1](i_panel, j_panel), normal_panel[2](i_panel, j_panel));
            UVLM::Geometry::convert_to_panel_coordinate_system(zeta[0].template block<2,2>(i_panel, j_panel),
                                                        zeta[1].template block<2,2>(i_panel, j_panel),
                                                        zeta[2].template block<2,2>(i_panel, j_panel),
                                                        longitudinal_panel_vec,
                                                        perpendicular_panel_vec,
                                                        normal_panel_vec,
                                                        panel_coordinates_epsilon,
                                                        panel_coordinates_eta,
                                                        panel_coordinates_z
                                                        );
			// Compute parts of the equations (10.95
            delta_epsilon = UVLM::Geometry::get_vector_diff(panel_coordinates_epsilon);
            delta_eta_vec = UVLM::Geometry::get_vector_diff(panel_coordinates_eta);

			// Checks if the source panel has a triangular shape (e.g. nose or tail panels)
            UVLM::Geometry::check_for_quadrilateral_panel(delta_epsilon, delta_eta_vec, flag_triangle, ignore_index);
			if (flag_triangle)
			{
				UVLM::Types::remove_row_from_VectorX(panel_coordinates_epsilon, ignore_index);
				delta_epsilon = UVLM::Geometry::get_vector_diff(panel_coordinates_epsilon);

				UVLM::Types::remove_row_from_VectorX(panel_coordinates_eta, ignore_index);
				delta_eta_vec = UVLM::Geometry::get_vector_diff(panel_coordinates_eta);
			}
			
			//  Compute Part of Eq. 10.95 and 10.96 (Katz & Plotzkin, 2001)
			d = (delta_epsilon.array().pow(2) + delta_eta_vec.array().pow(2)).sqrt();
			S_vec = delta_eta_vec.array()/d.array();
			C_vec = delta_epsilon.array()/d.array();
			for (unsigned int i_col=0; i_col<rows_collocation; ++i_col)
			{
				for (unsigned int j_col=0; j_col<cols_collocation; ++j_col)
				{

					if((same_surface) && (i_col==i_panel)&&(j_col==j_panel))
					{
						induced_velocity_vec[2] = 0.5;
						induced_velocity_vec[0] = 0.;
						induced_velocity_vec[1] = 0.;
					}
					else
					{
						UVLM::Geometry::convert_to_panel_coordinate_system(target_surface[0](i_col, j_col),
																target_surface[1](i_col, j_col),
																target_surface[2](i_col, j_col),
																longitudinal_panel_vec,
																perpendicular_panel_vec,
																normal_panel_vec,
																collocation_point_transf);
						collocation_point_transf[2] -= panel_coordinates_z[0];
					
						// Computes Eq. 10.92 (Katz & Plotzkin, 2001)
						delta_epsilon_x_vec = panel_coordinates_epsilon.array() - collocation_point_transf[0];
						delta_eta_y_vec = panel_coordinates_eta.array()- collocation_point_transf[1];
						r = (delta_epsilon_x_vec.array().pow(2) + delta_eta_y_vec.array().pow(2)+collocation_point_transf[2]*collocation_point_transf[2]).sqrt();
					
						//  Compute Part of Eq. 10.95 and 10.96 (Katz & Plotzkin, 2001)
						UVLM::Sources::get_Q(r,
											d,
											Q);

						//  Compute Eq. 10.95 (Katz & Plotzkin, 2001)
						induced_velocity_vec[0] = UVLM::Sources::dot_product(S_vec, Q)*UVLM::Constants::INV_PI4;
						//  Compute Eq. 10.96 (Katz & Plotzkin, 2001)
						induced_velocity_vec[1] = -UVLM::Sources::dot_product(C_vec, Q)*UVLM::Constants::INV_PI4;	
						// Compute Eq. 10.97 (Katz & Plotzkin, 2001)
						if (abs(collocation_point_transf[2])!= 0.0)
						{
							induced_velocity_vec[2] = UVLM::Sources::get_j_sum(r,
																				delta_epsilon,
																				delta_eta_vec,
																				delta_epsilon_x_vec,
																				delta_eta_y_vec,
																				collocation_point_transf[2]);
							induced_velocity_vec[2] *= UVLM::Constants::INV_PI4;
				
						}
						else
						{					
							induced_velocity_vec[2] = 0.0;
						}

						// convert induced panel velocity from panel coordinate system to collocation point coorindate system 
						longitudinal_col_vec = UVLM::Types::Vector3(longitudinal_col[0](i_col, j_col), longitudinal_col[1](i_col, j_col), longitudinal_col[2](i_col, j_col));
						perpendicular_col_vec = UVLM::Types::Vector3(perpendicular_col[0](i_col, j_col), perpendicular_col[1](i_col, j_col), perpendicular_col[2](i_col, j_col));
						normal_col_vec = UVLM::Types::Vector3(normal_col[0](i_col, j_col), normal_col[1](i_col, j_col), normal_col[2](i_col, j_col));

						UVLM::Geometry::convert_from_panel_A_to_panel_B_coordinate_system(induced_velocity_vec,
																						  longitudinal_panel_vec,
																						  perpendicular_panel_vec,
																						  normal_panel_vec,
																						  longitudinal_col_vec,
																						  perpendicular_col_vec,
																						  normal_col_vec
																						  );
					}
					// Store AIC coefficients in AIC matrices
					aic_x(collocation_id,panel_id) = induced_velocity_vec[0];
					aic_y(collocation_id,panel_id) = induced_velocity_vec[1];
					aic_z(collocation_id,panel_id) = induced_velocity_vec[2];
					collocation_id += 1;
				}
			}
        panel_id++;
        }
    }
}

void UVLM::Sources::get_Q(
    UVLM::Types::VectorX& r,
    UVLM::Types::VectorX& d,
    UVLM::Types::VectorX& Q
)
{
	const uint size_vector = r.rows();
	Q.resize(size_vector);
    for (uint i_row=0; i_row < size_vector; ++i_row)
    {
		if(i_row != size_vector-1)
		{
			Q[i_row] = log((r[i_row]+r[i_row+1]-d[i_row])
							  /(r[i_row]+r[i_row+1]+d[i_row]));
		}
		else
		{
			Q[i_row] = log((r[i_row]+r[0]-d[i_row])
						/(r[i_row]+r[0]+d[i_row]));
		}
    }
}

double UVLM::Sources::get_j_sum
(
    const UVLM::Types::VectorX& r,
    const UVLM::Types::VectorX& diff_epsilon,
    const UVLM::Types::VectorX& diff_eta,
    const UVLM::Types::VectorX& delta_epsilon_x_vec,
    const UVLM::Types::VectorX& delta_eta_y_vec,
    const double z
)
{
	double j_sum = 0;
	const uint vector_size = r.size();
    UVLM::Types::VectorX m_vec = diff_eta.array()/diff_epsilon.array();
    UVLM::Types::VectorX e_vec = delta_epsilon_x_vec.array().pow(2)+z*z;
    UVLM::Types::VectorX h_vec = delta_epsilon_x_vec.array()* delta_eta_y_vec.array();

	UVLM::Types::VectorX e_2_vec = UVLM::Types::reorder_vector_by_pushback(e_vec, 1);
	UVLM::Types::VectorX h_2_vec = UVLM::Types::reorder_vector_by_pushback(h_vec, 1);
	UVLM::Types::VectorX r_2_vec = UVLM::Types::reorder_vector_by_pushback(r, 1);

    for (uint i_point = 0; i_point < vector_size; i_point++)
	{
		j_sum += atan((m_vec[i_point]*e_vec[i_point]-h_vec[i_point])
		                      /(z*r[i_point])
							 ) -atan((m_vec[i_point]*e_2_vec[i_point]-h_2_vec[i_point])
		                      /(z*r_2_vec[i_point]));
	        
	}
	return j_sum;
}

template <typename vector>
double UVLM::Sources::dot_product
(
     const vector& vec1,
     const vector& vec2
 )
 {
    double dot_product = 0.0;
    for (uint i=0; i<vec1.size(); ++i)
    {
        dot_product += vec1[i]*vec2[i];
    }
    return dot_product;
}

