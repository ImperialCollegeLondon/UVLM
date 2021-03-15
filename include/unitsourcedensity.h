#pragma once

#include "EigenInclude.h"
#include "types.h"
#include "debugutils.h"

#include <math.h>
#include <cmath>


// Declaration for parallel computing
//#pragma omp declare reduction (sum_Vector3 : UVLM::Types::Vector3 : omp_out += omp_in) initializer(omp_priv = UVLM::Types::zeroVector3())

namespace UVLM
{
    namespace UnitSourceDensity
    {
        // DECLARATIONS
        template <typename t_zeta,
                  typename t_tsurface,
                  typename t_aic,
                  typename t_surf_vec_panel,
                  typename t_surf_vec_col>
        void get_influence_coefficient
        (
            const t_zeta&       zeta,
            const t_tsurface&   target_surface,
            t_aic&             aic_x,
            t_aic&             aic_y,
            t_aic&             aic_z,
            const t_surf_vec_panel&    longitudinal_panel = NULL,
            const t_surf_vec_panel&    perpendicular_panel = NULL,
            const t_surf_vec_panel&    normal_panel = NULL,
            const t_surf_vec_col&    longitudinal_col = NULL,
            const t_surf_vec_col&    perpendicular_col = NULL,
            const t_surf_vec_col&    normal_col = NULL,
			const bool& same_surface = false
        );

        void get_q_vec
        (
            UVLM::Types::VectorX& radius_vec,
            UVLM::Types::VectorX& d_vec,
            UVLM::Types::VectorX& Q_vec
        );

		double get_j_sum
		(
			const UVLM::Types::VectorX& radius_vec,
			const UVLM::Types::VectorX& diff_epsilon,
			const UVLM::Types::VectorX& diff_eta,
			const UVLM::Types::VectorX& delta_epsilon_x_vec,
			const UVLM::Types::VectorX& delta_eta_y_vec,
			const double z
		);

        template <typename vector>
        double dot_product
        (
             const vector& vec1,
             const vector& vec2
        );

		template <typename t_S_vec,
			  typename t_C_vec,
			  typename t_delta_eps_x_vec,
			  typename t_delta_eta_y_vec>
		bool check_if_col_on_panel
		(
			const t_S_vec& S_vec,
			const t_C_vec& C_vec,
			const t_delta_eps_x_vec& delta_epsilon_x_vec,
			const t_delta_eta_y_vec& delta_eta_y_vec
		);
    }
}
template <typename t_zeta,
		  typename t_tsurface,
		  typename t_aic,
		  typename t_surf_vec_panel,
		  typename t_surf_vec_col>

void UVLM::UnitSourceDensity::get_influence_coefficient
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
		
    UVLM::Types::VectorX delta_epsilon_vec;
    UVLM::Types::VectorX delta_eta_vec;
	
	UVLM::Types::Vector3 collocation_point_transf;
	UVLM::Types::Vector3 induced_velocity_vec;
	UVLM::Types::VectorX delta_epsilon_x_vec;
	UVLM::Types::VectorX delta_eta_y_vec;
	UVLM::Types::VectorX radius_vec;
	UVLM::Types::VectorX d_vec;
	UVLM::Types::VectorX S_vec;
	UVLM::Types::VectorX C_vec;
	UVLM::Types::VectorX Q_vec;
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
            delta_epsilon_vec = UVLM::Geometry::get_vector_diff(panel_coordinates_epsilon);
            delta_eta_vec = UVLM::Geometry::get_vector_diff(panel_coordinates_eta);
            UVLM::Geometry::check_for_quadrilateral_panel(delta_epsilon_vec, delta_eta_vec, flag_triangle, ignore_index);
			if (flag_triangle)
			{
				UVLM::Types::remove_row_from_VectorX(panel_coordinates_epsilon, ignore_index);
				delta_epsilon_vec = UVLM::Geometry::get_vector_diff(panel_coordinates_epsilon);

				UVLM::Types::remove_row_from_VectorX(panel_coordinates_eta, ignore_index);
				delta_eta_vec = UVLM::Geometry::get_vector_diff(panel_coordinates_eta);
			}
			d_vec = (delta_epsilon_vec.array().pow(2) + delta_eta_vec.array().pow(2)).sqrt();
			S_vec = delta_eta_vec.array()/d_vec.array();
			C_vec = delta_epsilon_vec.array()/d_vec.array();
			for (unsigned int i_col=0; i_col<rows_collocation; ++i_col)
			{
				for (unsigned int j_col=0; j_col<cols_collocation; ++j_col)
				{
					UVLM::Geometry::convert_to_panel_coordinate_system(target_surface[0](i_col, j_col),
															target_surface[1](i_col, j_col),
															target_surface[2](i_col, j_col),
															longitudinal_panel_vec,
															perpendicular_panel_vec,
															normal_panel_vec,
															collocation_point_transf);
					collocation_point_transf[2] -= panel_coordinates_z[0];
				
					delta_epsilon_x_vec = panel_coordinates_epsilon.array() - collocation_point_transf[0];
					delta_eta_y_vec = panel_coordinates_eta.array()- collocation_point_transf[1];
					radius_vec = (delta_epsilon_x_vec.array().pow(2) + delta_eta_y_vec.array().pow(2)+collocation_point_transf[2]*collocation_point_transf[2]).sqrt();
				
					UVLM::UnitSourceDensity::get_q_vec(radius_vec,
													   d_vec,
													   Q_vec);

					if((same_surface) && (i_col==i_panel)&&(j_col==j_panel))
					{
						induced_velocity_vec[0] = 0;
						induced_velocity_vec[1] = 0; 
						induced_velocity_vec[2] = 0.5;
					}
					else
					{
						
						induced_velocity_vec[0] = UVLM::UnitSourceDensity::dot_product(S_vec, Q_vec)*UVLM::Constants::INV_PI4;
						induced_velocity_vec[1] = -UVLM::UnitSourceDensity::dot_product(C_vec, Q_vec)*UVLM::Constants::INV_PI4;
						induced_velocity_vec[2] = 0;
						if (abs(collocation_point_transf[2])!= 0.0)
						{
							induced_velocity_vec[2] = UVLM::UnitSourceDensity::get_j_sum(radius_vec,
																						delta_epsilon_vec,
																						delta_eta_vec,
																						delta_epsilon_x_vec,
																						delta_eta_y_vec,
																						collocation_point_transf[2]);
							induced_velocity_vec[2] *= UVLM::Constants::INV_PI4;
				
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
					aic_x(collocation_id,panel_id) = induced_velocity_vec[0];//(4.0*UVLM::Constants::PI);
					aic_y(collocation_id,panel_id) = induced_velocity_vec[1];//(4.0*UVLM::Constants::PI);
					aic_z(collocation_id,panel_id) = induced_velocity_vec[2];//(4.0*UVLM::Constants::PI);
					collocation_id += 1;
				}
			}
        panel_id++;
        }
    }
}

void UVLM::UnitSourceDensity::get_q_vec(
    UVLM::Types::VectorX& radius_vec,
    UVLM::Types::VectorX& d_vec,
    UVLM::Types::VectorX& Q_vec
)
{
	const uint size_vector = radius_vec.rows();
	Q_vec.resize(size_vector);
    for (uint i_row=0; i_row < size_vector; ++i_row)
    {
		if(i_row != size_vector-1)
		{
			Q_vec[i_row] = log((radius_vec[i_row]+radius_vec[i_row+1]-d_vec[i_row])
							  /(radius_vec[i_row]+radius_vec[i_row+1]+d_vec[i_row]));
		}
		else
		{
			Q_vec[i_row] = log((radius_vec[i_row]+radius_vec[0]-d_vec[i_row])
						/(radius_vec[i_row]+radius_vec[0]+d_vec[i_row]));
		}
    }
}

double UVLM::UnitSourceDensity::get_j_sum
(
    const UVLM::Types::VectorX& radius_vec,
    const UVLM::Types::VectorX& diff_epsilon,
    const UVLM::Types::VectorX& diff_eta,
    const UVLM::Types::VectorX& delta_epsilon_x_vec,
    const UVLM::Types::VectorX& delta_eta_y_vec,
    const double z
)
{
	double j_sum = 0;
	const uint vector_size = radius_vec.size();
    UVLM::Types::VectorX m_vec = diff_eta.array()/diff_epsilon.array();
    UVLM::Types::VectorX e_vec = delta_epsilon_x_vec.array().pow(2)+z*z;
    UVLM::Types::VectorX h_vec = delta_epsilon_x_vec.array()* delta_eta_y_vec.array();

	UVLM::Types::VectorX e_2_vec = UVLM::Types::reorder_vector_by_pushback(e_vec, 1);
	UVLM::Types::VectorX h_2_vec = UVLM::Types::reorder_vector_by_pushback(h_vec, 1);
	UVLM::Types::VectorX r_2_vec = UVLM::Types::reorder_vector_by_pushback(radius_vec, 1);

    for (uint i_point = 0; i_point < vector_size; i_point++)
	{
		j_sum += atan((m_vec[i_point]*e_vec[i_point]-h_vec[i_point])
		                      /(z*radius_vec[i_point])
							 ) -atan((m_vec[i_point]*e_2_vec[i_point]-h_2_vec[i_point])
		                      /(z*r_2_vec[i_point]));
	        
	}
	return j_sum;
}

template <typename vector>
double UVLM::UnitSourceDensity::dot_product
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

template <typename t_S_vec,
		  typename t_C_vec,
		  typename t_delta_eps_x_vec,
		  typename t_delta_eta_y_vec>
bool UVLM::UnitSourceDensity::check_if_col_on_panel
(
	const t_S_vec& S_vec,
	const t_C_vec& C_vec,
	const t_delta_eps_x_vec& delta_epsilon_x_vec,
	const t_delta_eta_y_vec& delta_eta_y_vec
)
{
	// Check Hess and Smith (1967) p. 54
	if (((-delta_epsilon_x_vec.array()*S_vec.array() + delta_eta_y_vec.array()*C_vec.array()).array() <= 0.0).any())
	{
		return 1;
	}
	else
	{
		return 0;
	}
}