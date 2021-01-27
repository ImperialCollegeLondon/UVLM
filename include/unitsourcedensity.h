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
                  typename t_u_induced_col,
                  typename t_uout,
                  typename t_surf_vec_panel,
                  typename t_surf_vec_col>
        void get_influence_coefficient
        (
            const t_zeta&       zeta,
            const t_tsurface&   target_surface,
			t_u_induced_col&   u_induced_col_surface_x,
			t_u_induced_col&   u_induced_col_surface_y,
			t_u_induced_col&   u_induced_col_surface_z,
            t_uout&             uout,
            const t_surf_vec_panel&    longitudinal_panel = NULL,
            const t_surf_vec_panel&    perpendicular_panel = NULL,
            const t_surf_vec_panel&    normal_panel = NULL,
            const t_surf_vec_col&    longitudinal_col = NULL,
            const t_surf_vec_col&    perpendicular_col = NULL,
            const t_surf_vec_col&    normal_col = NULL
        );

        void get_q_vec
        (
            const UVLM::Types::Vector4& radius_vec,
            const UVLM::Types::Vector4& d_vec,
            UVLM::Types::Vector4& Q_vec
        );

        void get_q_vec
        (
            const UVLM::Types::Vector3& radius_vec,
            const UVLM::Types::Vector3& d_vec,
            UVLM::Types::Vector3& Q_vec
        );

        void get_j_vec
        (
            const UVLM::Types::Vector4& radius_vec,
            const UVLM::Types::Vector4& S_vec,
            const UVLM::Types::Vector4& C_vec,
            const UVLM::Types::Vector4& delta_epsilon_x_vec,
            const UVLM::Types::Vector4& delta_eta_y_vec,
            const double z,
            UVLM::Types::Vector4& Q_vec
        );
        void get_j_vec
        (
            const UVLM::Types::Vector3& radius_vec,
            const UVLM::Types::Vector3& S_vec,
            const UVLM::Types::Vector3& C_vec,
            const UVLM::Types::Vector3& delta_epsilon_x_vec,
            const UVLM::Types::Vector3& delta_eta_y_vec,
            const double z,
            UVLM::Types::Vector3& Q_vec
        );


		void get_j_vec_katz
		(
			const UVLM::Types::Vector3& radius_vec,
			const UVLM::Types::Vector3& diff_epsilon,
			const UVLM::Types::Vector3& diff_eta,
			const UVLM::Types::Vector3& delta_epsilon_x_vec,
			const UVLM::Types::Vector3& delta_eta_y_vec,
			const double z,
			UVLM::Types::Vector3& J_vec
		);

		void get_j_vec_katz
		(
			const UVLM::Types::Vector4& radius_vec,
			const UVLM::Types::Vector4& diff_epsilon,
			const UVLM::Types::Vector4& diff_eta,
			const UVLM::Types::Vector4& delta_epsilon_x_vec,
			const UVLM::Types::Vector4& delta_eta_y_vec,
			const double z,
			UVLM::Types::Vector4& J_vec
		);
		
        template <typename vector>
        double dot_product
        (
             const vector& vec1,
             const vector& vec2
        );
		template <typename t_panel_z,
				  typename vector_in,
				  typename t_tsurface,
				  typename t_u_induced_col,
				  typename t_uout,
                  typename t_longitudinals,
                  typename t_perpendiculars,
                  typename t_normals>
		void get_Aij_quadrilateral
		(
			const vector_in& panel_coordinates_epsilon,
			const vector_in& panel_coordinates_eta,
			const t_panel_z& panel_coordinate_z,
			const uint rows_collocation,
			const uint cols_collocation,
			const t_tsurface & target_surface,
			const UVLM::Types::Vector3&    longitudinal_panel,
			const UVLM::Types::Vector3&    perpendicular_panel,
			const UVLM::Types::Vector3&    normal_panel,
			const UVLM::Types::Vector4& delta_epsilon_vec,
			const UVLM::Types::Vector4& delta_eta_vec,
			t_u_induced_col&   u_induced_col_surface_x,
			t_u_induced_col&   u_induced_col_surface_y,
			t_u_induced_col&   u_induced_col_surface_z,
			t_uout&             uout,
			const uint panel_id,
			uint& collocation_id,
			const unsigned int i_panel,
			const unsigned int j_panel,
			const t_longitudinals&    longitudinal,
			const t_perpendiculars&    perpendicular,
			const t_normals&    normal
		);
		template <typename t_panel_z,
				  typename vector_in,
				  typename t_tsurface,
				  typename t_u_induced_col,
				  typename t_uout,
                  typename t_longitudinals,
                  typename t_perpendiculars,
                  typename t_normals>
		void get_Aij_triangle
		(
			const vector_in& panel_coordinates_epsilon,
			const vector_in& panel_coordinates_eta,
			const t_panel_z& panel_coordinate_z,
			const uint rows_collocation,
			const uint cols_collocation,
			const t_tsurface & target_surface,
			const UVLM::Types::Vector3&    longitudinal_panel,
			const UVLM::Types::Vector3&    perpendicular_panel,
			const UVLM::Types::Vector3&    normal_panel,
			t_u_induced_col&   u_induced_col_surface_x,
			t_u_induced_col&   u_induced_col_surface_y,
			t_u_induced_col&   u_induced_col_surface_z,
			t_uout&             uout,
			const uint panel_id,
			uint& collocation_id,
			const unsigned int i_panel,
			const unsigned int j_panel,
            const t_longitudinals&    longitudinal = NULL,
            const t_perpendiculars&    perpendicular = NULL,
            const t_normals&    normal = NULL
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
		  typename t_u_induced_col,
		  typename t_uout,
		  typename t_surf_vec_panel,
		  typename t_surf_vec_col>

void UVLM::UnitSourceDensity::get_influence_coefficient
    (
    const t_zeta&       zeta,
    const t_tsurface&   target_surface,
    t_u_induced_col&   u_induced_col_surface_x,
    t_u_induced_col&   u_induced_col_surface_y,
    t_u_induced_col&   u_induced_col_surface_z,
    t_uout&             uout,
	const t_surf_vec_panel&    longitudinal_panel,
	const t_surf_vec_panel&    perpendicular_panel,
	const t_surf_vec_panel&    normal_panel,
	const t_surf_vec_col&    longitudinal_col,
	const t_surf_vec_col&    perpendicular_col,
	const t_surf_vec_col&    normal_col
    )
    {

    const unsigned int rows_collocation = target_surface[0].rows();
    const unsigned int cols_collocation = target_surface[0].cols();
    const unsigned int rows_panel = zeta[0].rows()-1;
    const unsigned int cols_panel = zeta[0].cols()-1;

	UVLM::Types::Vector4 panel_coordinates_epsilon;
	UVLM::Types::Vector4 panel_coordinates_eta;
	UVLM::Types::Vector4 panel_coordinates_z;
    UVLM::Types::Vector4 delta_epsilon_vec;
    UVLM::Types::Vector4 delta_eta_vec;
    UVLM::Types::Vector3 normal_panel_vec;
    UVLM::Types::Vector3 longitudinal_panel_vec;
    UVLM::Types::Vector3 perpendicular_panel_vec;

	bool flag_triangle = false;
	int ignore_index;
    //PANELS
    uint panel_id = 0;
    uint collocation_id = 0;
    for (unsigned int i_panel=0; i_panel<rows_panel; ++i_panel)
    {
        for (unsigned int j_panel=0; j_panel<cols_panel; ++j_panel)
        {
            collocation_id = 0;
			UVLM::Types::Vector3 collocation_point_i = UVLM::Types::Vector3(target_surface[0](i_panel, j_panel),
													   target_surface[1](i_panel, j_panel),
													   target_surface[2](i_panel, j_panel));
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

            UVLM::Geometry::get_vector_diff(panel_coordinates_epsilon,
                                            delta_epsilon_vec);
            UVLM::Geometry::get_vector_diff(panel_coordinates_eta,
                                            delta_eta_vec);
            UVLM::Geometry::check_for_quadrilateral_panel(delta_epsilon_vec, delta_eta_vec, flag_triangle, ignore_index);
			if (flag_triangle)
			{
				UVLM::Types::Vector3 panel_coordinates_epsilon_triangle =UVLM::Types::remove_row(panel_coordinates_epsilon, ignore_index);
				UVLM::Types::Vector3 panel_coordinates_eta_triangle =UVLM::Types::remove_row(panel_coordinates_eta, ignore_index);
				UVLM::UnitSourceDensity::get_Aij_triangle
				(
					panel_coordinates_epsilon_triangle,
					panel_coordinates_eta_triangle,
					panel_coordinates_z[0],
					rows_collocation,
					cols_collocation,
					target_surface,
					longitudinal_panel_vec,
					perpendicular_panel_vec,
					normal_panel_vec,
					u_induced_col_surface_x,
					u_induced_col_surface_y,
					u_induced_col_surface_z,
					uout,
					panel_id,
					collocation_id,
					i_panel,
					j_panel,
					longitudinal_col,
					perpendicular_col,
					normal_col
				);
			}
			else
			{
				UVLM::UnitSourceDensity::get_Aij_quadrilateral
				(
					panel_coordinates_epsilon,
					panel_coordinates_eta,
					panel_coordinates_z[0],
					rows_collocation,
					cols_collocation,
					target_surface,
					longitudinal_panel_vec,
					perpendicular_panel_vec,
					normal_panel_vec,
					delta_epsilon_vec,
					delta_eta_vec,
					u_induced_col_surface_x,
					u_induced_col_surface_y,
					u_induced_col_surface_z,
					uout,
					panel_id,
					collocation_id,
					i_panel,
					j_panel,
					longitudinal_col,
					perpendicular_col,
					normal_col
				);
			}
        panel_id++;
        }
    }
}


template <typename t_panel_z,
		  typename vector_in,
		  typename t_tsurface,
		  typename t_u_induced_col,
		  typename t_uout,
			  typename t_longitudinals,
			  typename t_perpendiculars,
			  typename t_normals>
void UVLM::UnitSourceDensity::get_Aij_quadrilateral
(
	const vector_in& panel_coordinates_epsilon,
	const vector_in& panel_coordinates_eta,
	const t_panel_z& panel_coordinate_z,
	const uint rows_collocation,
	const uint cols_collocation,
	const t_tsurface & target_surface,
	const UVLM::Types::Vector3&    longitudinal_panel,
    const UVLM::Types::Vector3&    perpendicular_panel,
    const UVLM::Types::Vector3&    normal_panel,
	const UVLM::Types::Vector4& delta_epsilon_vec,
	const UVLM::Types::Vector4& delta_eta_vec,
	t_u_induced_col&   u_induced_col_surface_x,
    t_u_induced_col&   u_induced_col_surface_y,
    t_u_induced_col&   u_induced_col_surface_z,
    t_uout&             uout,
	const uint panel_id,
	uint& collocation_id,
	const unsigned int i_panel,
	const unsigned int j_panel,
    const t_longitudinals&    longitudinal,
    const t_perpendiculars&    perpendicular,
    const t_normals&    normal
)
{
	UVLM::Types::Vector4 d_vec = (delta_epsilon_vec.array().pow(2) + delta_eta_vec.array().pow(2)).sqrt();
	UVLM::Types::Vector4 S_vec = delta_eta_vec.array()/d_vec.array();
	UVLM::Types::Vector4 C_vec = delta_epsilon_vec.array()/d_vec.array();
	UVLM::Types::Vector4 Q_vec;
	UVLM::Types::Vector4 J_vec;
	UVLM::Types::Vector4 radius_vec;
	UVLM::Types::Vector4 delta_epsilon_x_vec;
	UVLM::Types::Vector4 delta_eta_y_vec;

    UVLM::Types::Vector3 collocation_point_transf;
    UVLM::Types::Vector3 induced_velocity_vec;
	
	bool flag_if_col_on_panel;
	for (unsigned int i_col=0; i_col<rows_collocation; ++i_col)
	{
		for (unsigned int j_col=0; j_col<cols_collocation; ++j_col)
		{
			UVLM::Geometry::convert_to_panel_coordinate_system(target_surface[0](i_col, j_col),
															  target_surface[1](i_col, j_col),
															  target_surface[2](i_col, j_col),
															  longitudinal_panel,
															  perpendicular_panel,
															  normal_panel,
															  collocation_point_transf);
			collocation_point_transf[2] -= panel_coordinate_z;
			
			delta_epsilon_x_vec = panel_coordinates_epsilon.array() - collocation_point_transf[0];
			delta_eta_y_vec = panel_coordinates_eta.array()- collocation_point_transf[1];
			radius_vec = (delta_epsilon_x_vec.array().pow(2) + delta_eta_y_vec.array().pow(2)+  collocation_point_transf[2]*collocation_point_transf[2]).sqrt();
			UVLM::UnitSourceDensity::get_q_vec(radius_vec,
											   d_vec,
											   Q_vec);

			flag_if_col_on_panel = false;
			if (abs(collocation_point_transf[2])< 0.00000000001)
			{
				flag_if_col_on_panel =UVLM::UnitSourceDensity::check_if_col_on_panel(S_vec,
																					C_vec,
																					delta_epsilon_x_vec,
																					delta_eta_y_vec); 
			    if(flag_if_col_on_panel)
					{
					induced_velocity_vec[0] = 0;
					induced_velocity_vec[1] = 0; 
					induced_velocity_vec[2] = 2.0*UVLM::Constants::PI;
				}
			}
			if (!flag_if_col_on_panel)
			{
				induced_velocity_vec[0] = UVLM::UnitSourceDensity::dot_product(S_vec, Q_vec);
				induced_velocity_vec[1] = -UVLM::UnitSourceDensity::dot_product(C_vec, Q_vec);;
				induced_velocity_vec[2] = 0;
				if (abs(collocation_point_transf[2])!=0.0)
				{
				UVLM::UnitSourceDensity::get_j_vec_katz
					(
						radius_vec,
						delta_epsilon_vec,
						delta_eta_vec,
						delta_epsilon_x_vec,
						delta_eta_y_vec,
						collocation_point_transf[2],
						J_vec
					);
					induced_velocity_vec[2] = (J_vec[0]+J_vec[1]+J_vec[2]+J_vec[3]);
				}
				// convert induced panel velocity from panel coordinate system to collocation point coorindate system 
				UVLM::Types::Vector3 longitudinal_col = UVLM::Types::Vector3(longitudinal[0](i_col, j_col), longitudinal[1](i_col, j_col), longitudinal[2](i_col, j_col));
				UVLM::Types::Vector3 perpendicular_col = UVLM::Types::Vector3(perpendicular[0](i_col, j_col), perpendicular[1](i_col, j_col), perpendicular[2](i_col, j_col));
				UVLM::Types::Vector3 normal_col = UVLM::Types::Vector3(normal[0](i_col, j_col), normal[1](i_col, j_col), normal[2](i_col, j_col));
				UVLM::Geometry::convert_from_panel_A_to_panel_B_coordinate_system
				(
					induced_velocity_vec,
					longitudinal_panel,
					perpendicular_panel,
					normal_panel,
					longitudinal_col,
					perpendicular_col,
					normal_col
				);
			}
			uout(collocation_id,panel_id) = induced_velocity_vec[2];
			u_induced_col_surface_x(collocation_id,panel_id) = induced_velocity_vec[0];
			u_induced_col_surface_y(collocation_id,panel_id) = induced_velocity_vec[1];
			u_induced_col_surface_z(collocation_id,panel_id) = induced_velocity_vec[2];
			collocation_id += 1;
		}
	}
}

template <typename t_panel_z,
		  typename vector_in,
		  typename t_tsurface,
		  typename t_u_induced_col,
          typename t_uout,
			  typename t_longitudinals,
			  typename t_perpendiculars,
			  typename t_normals>
void UVLM::UnitSourceDensity::get_Aij_triangle
(
	const vector_in& panel_coordinates_epsilon,
	const vector_in& panel_coordinates_eta,
	const t_panel_z& panel_coordinate_z,
	const uint rows_collocation,
	const uint cols_collocation,
	const t_tsurface & target_surface,
	const UVLM::Types::Vector3&    longitudinal_panel,
    const UVLM::Types::Vector3&    perpendicular_panel,
    const UVLM::Types::Vector3&    normal_panel,
	t_u_induced_col&   u_induced_col_surface_x,
    t_u_induced_col&   u_induced_col_surface_y,
    t_u_induced_col&   u_induced_col_surface_z,
    t_uout&             uout,
	const uint panel_id,
	uint& collocation_id,
	const unsigned int i_panel,
	const unsigned int j_panel,
    const t_longitudinals&    longitudinal,
    const t_perpendiculars&    perpendicular,
    const t_normals&    normal
)
{
	UVLM::Types::Vector3 delta_epsilon_vec;
	UVLM::Types::Vector3 delta_eta_vec;
	UVLM::Geometry::get_vector_diff(panel_coordinates_epsilon, delta_epsilon_vec);
	UVLM::Geometry::get_vector_diff(panel_coordinates_eta, delta_eta_vec);

	bool flag_if_col_on_panel;
	for (uint i_point = 0; i_point < 4; i_point++)
	{
		if ((delta_epsilon_vec[i_point]< 0) && (delta_epsilon_vec[i_point]> -0.00000001))
		{
			delta_epsilon_vec[i_point] = 0.0;
		}
	}
	UVLM::Types::Vector3 d_vec = (delta_epsilon_vec.array().pow(2) + delta_eta_vec.array().pow(2)).sqrt();
	UVLM::Types::Vector3 S_vec = delta_eta_vec.array()/d_vec.array();
	UVLM::Types::Vector3 C_vec = delta_epsilon_vec.array()/d_vec.array();
	// std::cout << "\n delta eps: \n" << delta_epsilon_vec << "\n delta eta: \n" << delta_eta_vec << std::endl;
	// std::cout << "\n d_vec: \n" << d_vec << std::endl;
	// std::cout << "\n S_vec: \n" << S_vec << std::endl;
	// std::cout << "\n C_vec: \n" << C_vec << std::endl;
	
	UVLM::Types::Vector3 Q_vec;
	UVLM::Types::Vector3 J_vec;
	UVLM::Types::Vector3 radius_vec;
    UVLM::Types::Vector3 induced_velocity_vec;
	
    UVLM::Types::Vector3 collocation_point_transf;
	UVLM::Types::Vector3 delta_epsilon_x_vec;
	UVLM::Types::Vector3 delta_eta_y_vec;
	for (unsigned int i_col=0; i_col<rows_collocation; ++i_col)
	{
		for (unsigned int j_col=0; j_col<cols_collocation; ++j_col)
		{
			UVLM::Geometry::convert_to_panel_coordinate_system(target_surface[0](i_col, j_col),
															  target_surface[1](i_col, j_col),
															  target_surface[2](i_col, j_col),
															  longitudinal_panel,
															  perpendicular_panel,
															  normal_panel,
															  collocation_point_transf);
			collocation_point_transf[2] -= panel_coordinate_z;
			delta_epsilon_x_vec = panel_coordinates_epsilon.array() - collocation_point_transf[0];
			delta_eta_y_vec = panel_coordinates_eta.array()- collocation_point_transf[1];
			radius_vec = (delta_epsilon_x_vec.array().pow(2) + delta_eta_y_vec.array().pow(2)+  collocation_point_transf[2]*collocation_point_transf[2]).sqrt();
			UVLM::UnitSourceDensity::get_q_vec(radius_vec,
											   d_vec,
											   Q_vec);
		

			flag_if_col_on_panel = false;
			if (abs(collocation_point_transf[2])< 0.00000000001)
			{
				flag_if_col_on_panel =UVLM::UnitSourceDensity::check_if_col_on_panel(S_vec,
																					C_vec,
																					delta_epsilon_x_vec,
																					delta_eta_y_vec); 
			    if(flag_if_col_on_panel)
					{
					induced_velocity_vec[0] = 0;
					induced_velocity_vec[1] = 0; 
					induced_velocity_vec[2] = 2.0*UVLM::Constants::PI;
				}
			}
			if (!flag_if_col_on_panel)
			{
				induced_velocity_vec[0] = UVLM::UnitSourceDensity::dot_product(S_vec, Q_vec);
				induced_velocity_vec[1] = -UVLM::UnitSourceDensity::dot_product(C_vec, Q_vec);
				induced_velocity_vec[2] = 0;
				if (abs(collocation_point_transf[2])!= 0.0)
				{
				UVLM::UnitSourceDensity::get_j_vec_katz
				(
					radius_vec,
					delta_epsilon_vec,
					delta_eta_vec,
					delta_epsilon_x_vec,
					delta_eta_y_vec,
					collocation_point_transf[2],
					J_vec
				);
				induced_velocity_vec[2] = (J_vec[0]+J_vec[1]+J_vec[2]);
				
				}
				// convert induced panel velocity from panel coordinate system to collocation point coorindate system 
				UVLM::Types::Vector3 longitudinal_col = UVLM::Types::Vector3(longitudinal[0](i_col, j_col), longitudinal[1](i_col, j_col), longitudinal[2](i_col, j_col));
				UVLM::Types::Vector3 perpendicular_col = UVLM::Types::Vector3(perpendicular[0](i_col, j_col), perpendicular[1](i_col, j_col), perpendicular[2](i_col, j_col));
				UVLM::Types::Vector3 normal_col = UVLM::Types::Vector3(normal[0](i_col, j_col), normal[1](i_col, j_col), normal[2](i_col, j_col));

				UVLM::Geometry::convert_from_panel_A_to_panel_B_coordinate_system(induced_velocity_vec,
																  longitudinal_panel,
																  perpendicular_panel,
																  normal_panel,
																  longitudinal_col,
																  perpendicular_col,
																  normal_col
																  );
			}
			uout(collocation_id,panel_id) = induced_velocity_vec[2];
			u_induced_col_surface_x(collocation_id,panel_id) = induced_velocity_vec[0];
			u_induced_col_surface_y(collocation_id,panel_id) = induced_velocity_vec[1];
			u_induced_col_surface_z(collocation_id,panel_id) = induced_velocity_vec[2];
			collocation_id++;
		}
	}
}

void UVLM::UnitSourceDensity::get_q_vec(
    const UVLM::Types::Vector4& radius_vec,
    const UVLM::Types::Vector4& d_vec,
    UVLM::Types::Vector4& Q_vec
)
{
    std::vector<int> combinations = {0, 1, 2, 3, 0};
    for (uint i_row=0; i_row < 4; ++i_row)
    {
        Q_vec[i_row] = log((radius_vec[combinations[i_row]]+radius_vec[combinations[i_row+1]]-d_vec[i_row])
						  /(radius_vec[combinations[i_row]]+radius_vec[combinations[i_row+1]]+d_vec[i_row]));
    }
}

void UVLM::UnitSourceDensity::get_q_vec(
    const UVLM::Types::Vector3& radius_vec,
    const UVLM::Types::Vector3& d_vec,
    UVLM::Types::Vector3& Q_vec
)
{
    std::vector<int> combinations = {0, 1, 2, 0};
    for (uint i_row=0; i_row < 3; ++i_row)
    {
        Q_vec[i_row] = log((radius_vec[combinations[i_row]]+radius_vec[combinations[i_row+1]]-d_vec[i_row])
						  /(radius_vec[combinations[i_row]]+radius_vec[combinations[i_row+1]]+d_vec[i_row]));
    }
}

void UVLM::UnitSourceDensity::get_j_vec
(
    const UVLM::Types::Vector4& radius_vec,
    const UVLM::Types::Vector4& S_vec,
    const UVLM::Types::Vector4& C_vec,
    const UVLM::Types::Vector4& delta_epsilon_x_vec,
    const UVLM::Types::Vector4& delta_eta_y_vec,
    const double z,
    UVLM::Types::Vector4& J_vec
)
{
    UVLM::Types::Vector4 s_1_vec = delta_epsilon_x_vec.array()*C_vec.array() + delta_eta_y_vec.array()*S_vec.array();
    UVLM::Types::Vector4 s_2_vec = UVLM::Types::Vector4(delta_epsilon_x_vec[1], delta_epsilon_x_vec[2], delta_epsilon_x_vec[3], delta_epsilon_x_vec[0]).array()*C_vec.array()
                                 + UVLM::Types::Vector4(delta_eta_y_vec[1], delta_eta_y_vec[2], delta_eta_y_vec[3], delta_eta_y_vec[0]).array()*S_vec.array();
    UVLM::Types::Vector4 R_vec = - delta_epsilon_x_vec.array()*S_vec.array() + delta_eta_y_vec.array()*C_vec.array();

    UVLM::Types::Vector4 radius_2_vec =UVLM::Types::Vector4(radius_vec[1], radius_vec[2], radius_vec[3], radius_vec[0]);
    J_vec = atan(R_vec.array()*abs(z)*(radius_vec.array()*s_2_vec.array()-radius_2_vec.array()*s_1_vec.array())
          /(radius_vec.array()*radius_2_vec.array()*R_vec.array()*R_vec.array() + z*z*s_2_vec.array()*s_1_vec.array()));
}

void UVLM::UnitSourceDensity::get_j_vec
(
    const UVLM::Types::Vector3& radius_vec,
    const UVLM::Types::Vector3& S_vec,
    const UVLM::Types::Vector3& C_vec,
    const UVLM::Types::Vector3& delta_epsilon_x_vec,
    const UVLM::Types::Vector3& delta_eta_y_vec,
    const double z,
    UVLM::Types::Vector3& J_vec
)
{
    UVLM::Types::Vector3 s_1_vec = delta_epsilon_x_vec.array()*C_vec.array() + delta_eta_y_vec.array()*S_vec.array();
    UVLM::Types::Vector3 s_2_vec = UVLM::Types::Vector3(delta_epsilon_x_vec[1], delta_epsilon_x_vec[2], delta_epsilon_x_vec[0]).array()*C_vec.array()
                                 + UVLM::Types::Vector3(delta_eta_y_vec[1], delta_eta_y_vec[2], delta_eta_y_vec[0]).array()*S_vec.array();
    UVLM::Types::Vector3 R_vec = - delta_epsilon_x_vec.array()*S_vec.array() + delta_eta_y_vec.array()*C_vec.array();

    UVLM::Types::Vector3 radius_2_vec =UVLM::Types::Vector3(radius_vec[1], radius_vec[2], radius_vec[0]);
    J_vec = atan(R_vec.array()*abs(z)*(radius_vec.array()*s_2_vec.array()-radius_2_vec.array()*s_1_vec.array())
          /(radius_vec.array()*radius_2_vec.array()*R_vec.array()*R_vec.array() + z*z*s_2_vec.array()*s_1_vec.array()));
}



void UVLM::UnitSourceDensity::get_j_vec_katz
(
    const UVLM::Types::Vector4& radius_vec,
    const UVLM::Types::Vector4& diff_epsilon,
    const UVLM::Types::Vector4& diff_eta,
    const UVLM::Types::Vector4& delta_epsilon_x_vec,
    const UVLM::Types::Vector4& delta_eta_y_vec,
    const double z,
    UVLM::Types::Vector4& J_vec
)
{

    UVLM::Types::Vector4 m_vec = diff_eta.array()/diff_epsilon.array();
	//std::cout << "M_VEC =  \n" << m_vec << std::endl;
    UVLM::Types::Vector4 e_vec = delta_epsilon_x_vec.array().pow(2)+z*z;
    UVLM::Types::Vector4 h_vec = delta_epsilon_x_vec.array()* delta_eta_y_vec.array();
    //UVLM::Types::Vector4 s_1_vec = delta_epsilon_x_vec.array()*C_vec.array() + delta_eta_y_vec.array()*S_vec.array();
    //UVLM::Types::Vector4 s_2_vec = UVLM::Types::Vector4(delta_epsilon_x_vec[1], delta_epsilon_x_vec[2], delta_epsilon_x_vec[3], delta_epsilon_x_vec[0]).array()*C_vec.array()
    //                             + UVLM::Types::Vector4(delta_eta_y_vec[1], delta_eta_y_vec[2], delta_eta_y_vec[3], delta_eta_y_vec[0]).array()*S_vec.array();
    //UVLM::Types::Vector4 R_vec = - delta_epsilon_x_vec.array()*S_vec.array() + delta_eta_y_vec.array()*C_vec.array();
    UVLM::Types::Vector4 e_2_vec =UVLM::Types::Vector4(e_vec[1], e_vec[2], e_vec[3], e_vec[0]);
    UVLM::Types::Vector4 h_2_vec =UVLM::Types::Vector4(h_vec[1], h_vec[2], h_vec[3], h_vec[0]);
    UVLM::Types::Vector4 r_2_vec =UVLM::Types::Vector4(radius_vec[1], radius_vec[2], radius_vec[3], radius_vec[0]);

    for (uint i_point = 0; i_point < 4; i_point++)
	{
		J_vec[i_point] = atan((m_vec[i_point]*e_vec[i_point]-h_vec[i_point])
		                      /(z*radius_vec[i_point])
							 ) -atan((m_vec[i_point]*e_2_vec[i_point]-h_2_vec[i_point])
		                      /(z*r_2_vec[i_point]));
	        
	}
//	std::cout << "\nPoint z = " << z << ",  |z| = " << abs(z) << std::endl;
//	J_vec = atan(R_vec.array()*abs(z)*(radius_vec.array()*s_2_vec.array()-radius_2_vec.array()*s_1_vec.array())
//          /(radius_vec.array()*radius_2_vec.array()*R_vec.array()*R_vec.array() + z*z*s_2_vec.array()*s_1_vec.array()));
}

void UVLM::UnitSourceDensity::get_j_vec_katz
(
    const UVLM::Types::Vector3& radius_vec,
    const UVLM::Types::Vector3& diff_epsilon,
    const UVLM::Types::Vector3& diff_eta,
    const UVLM::Types::Vector3& delta_epsilon_x_vec,
    const UVLM::Types::Vector3& delta_eta_y_vec,
    const double z,
    UVLM::Types::Vector3& J_vec
)
{
    UVLM::Types::Vector3 m_vec = diff_eta.array()/diff_epsilon.array();
    UVLM::Types::Vector3 e_vec = delta_epsilon_x_vec.array().pow(2)+z*z;
    UVLM::Types::Vector3 h_vec = delta_epsilon_x_vec.array()* delta_eta_y_vec.array();

    //UVLM::Types::Vector4 s_1_vec = delta_epsilon_x_vec.array()*C_vec.array() + delta_eta_y_vec.array()*S_vec.array();
    //UVLM::Types::Vector4 s_2_vec = UVLM::Types::Vector4(delta_epsilon_x_vec[1], delta_epsilon_x_vec[2], delta_epsilon_x_vec[3], delta_epsilon_x_vec[0]).array()*C_vec.array()
    //                             + UVLM::Types::Vector4(delta_eta_y_vec[1], delta_eta_y_vec[2], delta_eta_y_vec[3], delta_eta_y_vec[0]).array()*S_vec.array();
    //UVLM::Types::Vector4 R_vec = - delta_epsilon_x_vec.array()*S_vec.array() + delta_eta_y_vec.array()*C_vec.array();
    UVLM::Types::Vector3 e_2_vec =UVLM::Types::Vector3(e_vec[1], e_vec[2], e_vec[0]);
    UVLM::Types::Vector3 h_2_vec =UVLM::Types::Vector3(h_vec[1], h_vec[2], h_vec[0]);
    UVLM::Types::Vector3 r_2_vec =UVLM::Types::Vector3(radius_vec[1], radius_vec[2], radius_vec[0]);

    for (uint i_point = 0; i_point < 3; i_point++)
	{
		J_vec[i_point] = atan((m_vec[i_point]*e_vec[i_point]-h_vec[i_point])
		                      /(z*radius_vec[i_point])
							 ) -atan((m_vec[i_point]*e_2_vec[i_point]-h_2_vec[i_point])
		                      /(z*r_2_vec[i_point]));
	        
	}
//	std::cout << "\nPoint z = " << z << ",  |z| = " << abs(z) << std::endl;
//	J_vec = atan(R_vec.array()*abs(z)*(radius_vec.array()*s_2_vec.array()-radius_2_vec.array()*s_1_vec.array())
//          /(radius_vec.array()*radius_2_vec.array()*R_vec.array()*R_vec.array() + z*z*s_2_vec.array()*s_1_vec.array()));
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