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

        }
    }
}
