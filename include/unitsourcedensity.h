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
                  typename t_longitudinals,
                  typename t_perpendiculars,
                  typename t_normals>
        void get_influence_coefficient
        (
            const t_zeta&       zeta,
            const t_tsurface&   target_surface,
			t_u_induced_col&   u_induced_col_surface_x,
			t_u_induced_col&   u_induced_col_surface_y,
			t_u_induced_col&   u_induced_col_surface_z,
            t_uout&             uout,
            const t_longitudinals&    longitudinal = NULL,
            const t_perpendiculars&    perpendicular = NULL,
            const t_normals&    normal = NULL
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
				  typename t_uout>
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
			const unsigned int j_panel
		);
		template <typename t_panel_z,
				  typename vector_in,
				  typename t_tsurface,
				  typename t_u_induced_col,
				  typename t_uout>
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
			const unsigned int j_panel
		);
    }
}

template <typename t_zeta,
          typename t_tsurface,
          typename t_u_induced_col,
          typename t_uout,
          typename t_longitudinals,
          typename t_perpendiculars,
          typename t_normals>
void UVLM::UnitSourceDensity::get_influence_coefficient
    (
    const t_zeta&       zeta,
    const t_tsurface&   target_surface,
    t_u_induced_col&   u_induced_col_surface_x,
    t_u_induced_col&   u_induced_col_surface_y,
    t_u_induced_col&   u_induced_col_surface_z,
    t_uout&             uout,
    const t_longitudinals&    longitudinal,
    const t_perpendiculars&    perpendicular,
    const t_normals&    normal
    )
    {

    const unsigned int rows_collocation = target_surface[0].rows();
    const unsigned int cols_collocation = target_surface[0].cols();

	UVLM::Types::Vector4 panel_coordinates_epsilon;
	UVLM::Types::Vector4 panel_coordinates_eta;
	UVLM::Types::Vector4 panel_coordinates_z;
    UVLM::Types::Vector4 delta_epsilon_vec;
    UVLM::Types::Vector4 delta_eta_vec;
    UVLM::Types::Vector3 normal_panel;
    UVLM::Types::Vector3 longitudinal_panel;
    UVLM::Types::Vector3 perpendicular_panel;

	bool flag_triangle = false;
	int ignore_index;
    //PANELS
    uint panel_id = 0;
    uint collocation_id = 0;
    for (unsigned int i_panel=0; i_panel<rows_collocation; ++i_panel)
    {
        for (unsigned int j_panel=0; j_panel<cols_collocation; ++j_panel)
        {
            collocation_id = 0;
            longitudinal_panel = UVLM::Types::Vector3(longitudinal[0](i_panel, j_panel), longitudinal[1](i_panel, j_panel), longitudinal[2](i_panel, j_panel));
            perpendicular_panel = UVLM::Types::Vector3(perpendicular[0](i_panel, j_panel), perpendicular[1](i_panel, j_panel), perpendicular[2](i_panel, j_panel));
            normal_panel = UVLM::Types::Vector3(normal[0](i_panel, j_panel), normal[1](i_panel, j_panel), normal[2](i_panel, j_panel));
            UVLM::Geometry::convert_to_panel_coordinate_system(zeta[0].template block<2,2>(i_panel, j_panel),
                                                        zeta[1].template block<2,2>(i_panel, j_panel),
                                                        zeta[2].template block<2,2>(i_panel, j_panel),
                                                        longitudinal_panel,
                                                        perpendicular_panel,
                                                        normal_panel,
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
					longitudinal_panel,
					perpendicular_panel,
					normal_panel,
					u_induced_col_surface_x,
					u_induced_col_surface_y,
					u_induced_col_surface_z,
					uout,
					panel_id,
					collocation_id,
					i_panel,
					j_panel
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
					longitudinal_panel,
					perpendicular_panel,
					normal_panel,
					delta_epsilon_vec,
					delta_eta_vec,
					u_induced_col_surface_x,
					u_induced_col_surface_y,
					u_induced_col_surface_z,
					uout,
					panel_id,
					collocation_id,
					i_panel,
					j_panel
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
		  typename t_uout>
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
	const unsigned int j_panel
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

			induced_velocity_vec[0] = -UVLM::UnitSourceDensity::dot_product(S_vec, Q_vec);
			induced_velocity_vec[1] = UVLM::UnitSourceDensity::dot_product(C_vec, Q_vec);

			if (i_panel == i_col && j_panel == j_col)
			{
				induced_velocity_vec[2] = 2*UVLM::Constants::PI;
			}
			else
			{
				induced_velocity_vec[2] = 0;
			}
			if (abs(collocation_point_transf[2]) > 0.00001)
			{
				UVLM::UnitSourceDensity::get_j_vec(radius_vec,
					S_vec,
					C_vec,
					delta_epsilon_x_vec,
					delta_eta_y_vec,
					collocation_point_transf[2],
					J_vec);
				induced_velocity_vec[2] -= J_vec[0]+J_vec[1]+J_vec[2]+J_vec[3];
				if (collocation_point_transf[2] < 0.0)
				{
					induced_velocity_vec[2] *= -1;
				}
			}
			uout(panel_id, collocation_id) = UVLM::UnitSourceDensity::dot_product(normal_panel, induced_velocity_vec);
			u_induced_col_surface_x(panel_id, collocation_id) = induced_velocity_vec[0];
			u_induced_col_surface_y(panel_id, collocation_id) = induced_velocity_vec[1];
			u_induced_col_surface_z(panel_id, collocation_id) = induced_velocity_vec[2];
			collocation_id += 1;
		}
	}
}

template <typename t_panel_z,
		  typename vector_in,
		  typename t_tsurface,
		  typename t_u_induced_col,
          typename t_uout>
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
	const unsigned int j_panel
)
{
	UVLM::Types::Vector3 delta_epsilon_vec;
	UVLM::Types::Vector3 delta_eta_vec;
	UVLM::Geometry::get_vector_diff(panel_coordinates_epsilon, delta_epsilon_vec);
	UVLM::Geometry::get_vector_diff(panel_coordinates_eta, delta_eta_vec);
	UVLM::Types::Vector3 d_vec = (delta_epsilon_vec.array().pow(2) + delta_eta_vec.array().pow(2)).sqrt();
	UVLM::Types::Vector3 S_vec = delta_eta_vec.array()/d_vec.array();
	UVLM::Types::Vector3 C_vec = delta_epsilon_vec.array()/d_vec.array();
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

			induced_velocity_vec[0] = -UVLM::UnitSourceDensity::dot_product(S_vec, Q_vec);
			induced_velocity_vec[1] = UVLM::UnitSourceDensity::dot_product(C_vec, Q_vec);

			if (i_panel == i_col && j_panel == j_col)
			{
				induced_velocity_vec[2] = 2*UVLM::Constants::PI;
			}
			else
			{
				induced_velocity_vec[2] = 0;
			}
			if (abs(collocation_point_transf[2]) > 0.00001)
			{
				UVLM::UnitSourceDensity::get_j_vec(radius_vec,
					S_vec,
					C_vec,
					delta_epsilon_x_vec,
					delta_eta_y_vec,
					collocation_point_transf[2],
					J_vec);
				induced_velocity_vec[2] -= J_vec[0]+J_vec[1]+J_vec[2];
				if (collocation_point_transf[2] < 0.0)
				{
					induced_velocity_vec[2] *= -1;
				}
			}
			uout(panel_id, collocation_id) = UVLM::UnitSourceDensity::dot_product(normal_panel, induced_velocity_vec);
			u_induced_col_surface_x(panel_id, collocation_id) = induced_velocity_vec[0];
			u_induced_col_surface_y(panel_id, collocation_id) = induced_velocity_vec[1];
			u_induced_col_surface_z(panel_id, collocation_id) = induced_velocity_vec[2];
			collocation_id += 1;
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
        Q_vec[i_row] = radius_vec[combinations[i_row]]+radius_vec[combinations[i_row+1]];
        Q_vec[i_row] = log((Q_vec[i_row]+ d_vec[i_row])/(Q_vec[i_row] - d_vec[i_row]));
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
        Q_vec[i_row] = radius_vec[combinations[i_row]]+radius_vec[combinations[i_row+1]];
        Q_vec[i_row] = log((Q_vec[i_row]+ d_vec[i_row])/(Q_vec[i_row] - d_vec[i_row]));
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