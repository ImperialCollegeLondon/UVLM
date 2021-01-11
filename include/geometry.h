#pragma once

#include "EigenInclude.h"
#include <unsupported/Eigen/Splines>
#include "types.h"
#include "mapping.h"

#include <iostream>
#include <cmath>

namespace UVLM
{
    namespace Geometry
    {
        // calculates the area of a triangle given its
        // side lengths: a, b, and c.
        // It uses Heron's fomula:
        // s = 0.5*(a + b + c)
        // A = sqrt(s*(s-a)*(s-b)*(s-c))
        UVLM::Types::Real triangle_area
        (
            const UVLM::Types::Real& a,
            const UVLM::Types::Real& b,
            const UVLM::Types::Real& c
        )
        {
            UVLM::Types::Real s = 0.5*(a + b + c);
            return std::sqrt(s*(s - a)*(s - b)*(s - c));
        }



        // Calculates the area of a quadrilateral in 3D
        // The method used is:
        // 1) divide the quad with a diagonal from 0 to 2
        // 2) calculate area of resulting triangles
        // 3) divide the quad with a diagonal from 1 to 3
        // 4) calculate area of resulting triangles
        // 5) average the two areas
        template <typename t_block>
        UVLM::Types::Real panel_area
        (
            const t_block& x,
            const t_block& y,
            const t_block& z
        )
        {
            UVLM::Types::Real area = 0;
            // calculate side length
            UVLM::Types::VectorX sides;
            sides.resize(4);
            uint i_side = 0;
            for (uint i_side=0; i_side<4; ++i_side)
            {
                uint i_first = UVLM::Mapping::vortex_indices(i_side, 0);
                uint j_first = UVLM::Mapping::vortex_indices(i_side, 1);
                uint i_second = UVLM::Mapping::vortex_indices((i_side + 1) % 4, 0);
                uint j_second = UVLM::Mapping::vortex_indices((i_side + 1) % 4, 1);
                sides(i_side) = std::sqrt(
                    (x(i_second,j_second) - x(i_first,j_first))*(x(i_second,j_second) - x(i_first,j_first)) +
                    (y(i_second,j_second) - y(i_first,j_first))*(y(i_second,j_second) - y(i_first,j_first)) +
                    (z(i_second,j_second) - z(i_first,j_first))*(z(i_second,j_second) - z(i_first,j_first)));
            }

            // diagonal from 0 to 2
            UVLM::Types::Real diagonal = 0;
            diagonal = std::sqrt(
                    (x(1,1) - x(0,0))*(x(1,1) - x(0,0)) +
                    (y(1,1) - y(0,0))*(y(1,1) - y(0,0)) +
                    (z(1,1) - z(0,0))*(z(1,1) - z(0,0)));

            area += triangle_area(sides(0), sides(1), diagonal);
            area += triangle_area(sides(2), sides(3), diagonal);

            // diagonal from 1 to 3
            diagonal = std::sqrt(
                    (x(1,0) - x(0,1))*(x(1,0) - x(0,1)) +
                    (y(1,0) - y(0,1))*(y(1,0) - y(0,1)) +
                    (z(1,0) - z(0,1))*(z(1,0) - z(0,1)));

            area += triangle_area(sides(1), sides(2), diagonal);
            area += triangle_area(sides(0), sides(3), diagonal);
            area *= 0.5;

            return area;
        }

        template <typename type>
        void panel_longitudinal_vector(type& x,
                                       type& y,
                                       type& z,
                                       UVLM::Types::Vector3& longitudinal_vec
                                       )
        {
            longitudinal_vec = UVLM::Types::Vector3(x(1,1) - x(0,0),
                                                    y(1,1) - y(0,0),
                                                    z(1,1) - z(0,0));
            longitudinal_vec.normalize();
        }

        template <typename type>
        void panel_tangential_vector(type& x,
                                     type& y,
                                     type& z,
                                     UVLM::Types::Vector3& tangential_vec
                                     )
        {
            tangential_vec = UVLM::Types::Vector3(x(1,0) - x(0,1),
                                                  y(1,0) - y(0,1),
                                                  z(1,0) - z(0,1));
            tangential_vec.normalize();
        }



        template <typename type>
        void panel_normal(type& x,
                          type& y,
                          type& z,
                          UVLM::Types::Vector3& normal
                          )
        {
            // correction for left-oriented panels
            UVLM::Types::Vector3 v_01(x(0,1) - x(0,0),
                                      y(0,1) - y(0,0),
                                      z(0,1) - z(0,0));
            UVLM::Types::Vector3 v_03(x(1,0) - x(0,0),
                                      y(1,0) - y(0,0),
                                      z(1,0) - z(0,0));
            UVLM::Types::Vector3 diff = v_01.cross(v_03);

            UVLM::Types::Vector3 A(x(1,1) - x(0,0),
                                   y(1,1) - y(0,0),
                                   z(1,1) - z(0,0));

            UVLM::Types::Vector3 B(x(1,0) - x(0,1),
                                   y(1,0) - y(0,1),
                                   z(1,0) - z(0,1));

            // if (diff(2) < 0.0)
            // {
                normal = B.cross(A);
            // } else
            // {
                // normal = A.cross(B);
            // }
            normal.normalize();
        }

        template <typename type>
        void panel_normal(type& x,
                          type& y,
                          type& z,
                          UVLM::Types::Real& xnormal,
                          UVLM::Types::Real& ynormal,
                          UVLM::Types::Real& znormal
                          )
        {
            UVLM::Types::Vector3 A;
            panel_normal(x, y, z, A);
            xnormal = A(0);
            ynormal = A(1);
            znormal = A(2);
        }

        template <typename type_in,
                  typename type_out>
        void generate_surfaceNormal(const type_in& zeta,
                                    type_out& normal)
        {
            for (unsigned int i_surf=0; i_surf<zeta.size(); ++i_surf)
            {
                for (unsigned int i_dim=0; i_dim<zeta[i_surf].size(); i_dim++)
                {
                    const unsigned int M = zeta[i_surf][i_dim].rows() - 1;
                    const unsigned int N = zeta[i_surf][i_dim].cols() - 1;

                    for (unsigned int iM=0; iM<M; ++iM)
                    {
                        for (unsigned int jN=0; jN<N; ++jN)
                        {
                            UVLM::Types::Vector3 temp_normal;
                            panel_normal(zeta[i_surf][0].template block<2,2>(iM,jN),
                                         zeta[i_surf][1].template block<2,2>(iM,jN),
                                         zeta[i_surf][2].template block<2,2>(iM,jN),
                                         temp_normal);

                            normal[i_surf][0](iM,jN) = temp_normal[0];
                            normal[i_surf][1](iM,jN) = temp_normal[1];
                            normal[i_surf][2](iM,jN) = temp_normal[2];
                        }
                    }
                }
            }
        }
		template <typename type_in>
		void check_for_quadrilateral_panel(const type_in& delta_vec_1,
										   const type_in& delta_vec_2,
										   bool& flag_triangle,
										   int& ignored_index)
		{
			for (int i=0; i<delta_vec_1.size(); ++i)
			{
				flag_triangle = false;
				if ((abs(delta_vec_1[i]) < 0.00001) && (abs(delta_vec_2[i]) < 0.00001))
				{
					flag_triangle = true;
					ignored_index = i;
					break;
				}
			}
		}
        template <typename type_in,
                  typename type_out>
        void generate_surface_vectors(const type_in& zeta,
                                      type_out& normal,
                                      type_out& long_vec,
                                      type_out& perpendicular_vec)
        {
            for (unsigned int i_surf=0; i_surf<zeta.size(); ++i_surf)
            {
                for (unsigned int i_dim=0; i_dim<zeta[i_surf].size(); i_dim++)
                {
                    const unsigned int M = zeta[i_surf][i_dim].rows() - 1;
                    const unsigned int N = zeta[i_surf][i_dim].cols() - 1;

                    for (unsigned int iM=0; iM<M; ++iM)
                    {
                        for (unsigned int jN=0; jN<N; ++jN)
                        {
                            UVLM::Types::Vector3 temp_long_vec;
                            panel_longitudinal_vector(zeta[i_surf][0].template block<2,2>(iM,jN),
                                                      zeta[i_surf][1].template block<2,2>(iM,jN),
                                                      zeta[i_surf][2].template block<2,2>(iM,jN),
                                                      temp_long_vec);

                            long_vec[i_surf][0](iM,jN) = temp_long_vec[0];
                            long_vec[i_surf][1](iM,jN) = temp_long_vec[1];
                            long_vec[i_surf][2](iM,jN) = temp_long_vec[2];

                            UVLM::Types::Vector3 temp_tangential_vec;
                            panel_tangential_vector(zeta[i_surf][0].template block<2,2>(iM,jN),
                                                    zeta[i_surf][1].template block<2,2>(iM,jN),
                                                    zeta[i_surf][2].template block<2,2>(iM,jN),
                                                    temp_tangential_vec);

                            UVLM::Types::Vector3 temp_normal_vec;
                            temp_normal_vec = temp_tangential_vec.cross(temp_long_vec);
                            temp_normal_vec.normalize();
                            normal[i_surf][0](iM,jN) = temp_normal_vec[0];
                            normal[i_surf][1](iM,jN) = temp_normal_vec[1];
                            normal[i_surf][2](iM,jN) = temp_normal_vec[2];

                            UVLM::Types::Vector3 temp_perpendicular_vec;
                            temp_perpendicular_vec = temp_long_vec.cross(temp_normal_vec);
                            temp_perpendicular_vec.normalize();
                            perpendicular_vec[i_surf][0](iM,jN) = temp_perpendicular_vec[0];
                            perpendicular_vec[i_surf][1](iM,jN) = temp_perpendicular_vec[1];
                            perpendicular_vec[i_surf][2](iM,jN) = temp_perpendicular_vec[2];
                        }
                }
            }
        }

        template <typename type>
        void convert_to_panel_coordinate_system(const type& x_G,
                                                const type& y_G,
                                                const type& z_G,
                                                const UVLM::Types::Vector3& chordwise_vec,
                                                const UVLM::Types::Vector3& tangential_vec,
                                                const UVLM::Types::Vector3& normal_vec,
                                                UVLM::Types::Vector4& x_transf,
                                                UVLM::Types::Vector4& y_transf,
                                                UVLM::Types::Vector4& z_transf
                                                )
        {
            UVLM::Types::Vector4 x = UVLM::Types::Vector4(x_G(0,0), x_G(1, 0), x_G(1, 1), x_G(0, 1));
            UVLM::Types::Vector4 y = UVLM::Types::Vector4(y_G(0,0), y_G(1, 0), y_G(1, 1), y_G(0, 1));
            UVLM::Types::Vector4 z = UVLM::Types::Vector4(z_G(0,0), z_G(1, 0), z_G(1, 1), z_G(0, 1));
            x_transf = x *chordwise_vec[0] + y *chordwise_vec[1] + z *chordwise_vec[2];
            y_transf = x *tangential_vec[0] + y *tangential_vec[1] + z *tangential_vec[2];
            z_transf = x *normal_vec[0] + y *normal_vec[1] + z * normal_vec[2];
        }

        template <typename type>
        void convert_to_panel_coordinate_system(const type& x_G,
                                                const type& y_G,
                                                const type& z_G,
                                                const UVLM::Types::Vector3& chordwise_vec,
                                                const UVLM::Types::Vector3& tangential_vec,
                                                const UVLM::Types::Vector3& normal_vec,
                                                UVLM::Types::Vector3& point_transf
                                                )
        {
            point_transf[0] = x_G *chordwise_vec[0] + y_G *chordwise_vec[1] + z_G *chordwise_vec[2];
            point_transf[1] = x_G *tangential_vec[0] + y_G *tangential_vec[1] + z_G *tangential_vec[2];
            point_transf[2] = x_G *normal_vec[0] + y_G *normal_vec[1] + z_G * normal_vec[2];
        }
		


        template <typename t_in,
                  typename t_out>
        void generate_colocationMesh
        (
            t_in& vortex_mesh,
            t_out& collocation_mesh
        )
        {
            // Size of surfaces contained in a vector of tuples
            UVLM::Types::VecDimensions dimensions;
            dimensions.resize(vortex_mesh.size());
            for (unsigned int i_surf=0; i_surf<dimensions.size(); ++i_surf)
            {
                dimensions[i_surf] = UVLM::Types::IntPair(
                                                    vortex_mesh[i_surf][0].rows(),
                                                    vortex_mesh[i_surf][0].cols());
            }

            if (collocation_mesh.empty())
            {
                UVLM::Types::allocate_VecVecMat(collocation_mesh,
                                                UVLM::Constants::NDIM,
                                                dimensions,
                                                -1);
            }
            for (unsigned int i_surf=0; i_surf<dimensions.size(); ++i_surf)
            {
                UVLM::Mapping::BilinearMapping(vortex_mesh[i_surf],
                                               collocation_mesh[i_surf]);
            }
        }
        void get_vector_diff(const UVLM::Types::Vector4& vec,
                             UVLM::Types::Vector4& diff_vec)
        {
            // Calcualtes difference between adjascent vector scalars
            diff_vec = UVLM::Types::Vector4(vec[1]-vec[0],
                                            vec[2]-vec[1],
                                            vec[3]-vec[2],
                                            vec[0]-vec[3]);
        }
    } // geometry

    namespace Interpolation
    {
        template <typename t_dist,
                  typename t_dist_conv,
                  typename t_coord,
                  typename t_coord_conv>
        void linear
        (
            uint M,
            const t_dist& dist_to_orig,
            const t_dist_conv& dist_to_orig_conv,
            const t_coord_conv& coord0,
            const t_coord_conv& coord1,
            const t_coord_conv& coord2,
            t_coord& new_coord0,
            t_coord& new_coord1,
            t_coord& new_coord2
        )
        {
            UVLM::Types::Real to_prev, to_next, prev_to_next;
            uint i_conv=0;
            for (unsigned int i_m=0; i_m<M; ++i_m)
            {
                while ((dist_to_orig_conv(i_conv) <= dist_to_orig(i_m)) and (i_conv < M))
                {i_conv++;}

                to_prev = dist_to_orig(i_m) - dist_to_orig_conv(i_conv - 1);
                to_next = dist_to_orig_conv(i_conv) - dist_to_orig(i_m);
                prev_to_next = dist_to_orig_conv(i_conv) - dist_to_orig_conv(i_conv - 1);
                new_coord0(i_m) = (to_prev*coord0(i_conv) + to_next*coord0(i_conv - 1))/prev_to_next;
                new_coord1(i_m) = (to_prev*coord1(i_conv) + to_next*coord1(i_conv - 1))/prev_to_next;
                new_coord2(i_m) = (to_prev*coord2(i_conv) + to_next*coord2(i_conv - 1))/prev_to_next;
            }
        } // linear

        template <typename t_dist,
                  typename t_dist_conv,
                  typename t_coord,
                  typename t_coord_conv>
        void parabolic
        (
            uint M,
            const t_dist& dist_to_orig,
            const t_dist_conv& dist_to_orig_conv,
            const t_coord_conv& coord0,
            const t_coord_conv& coord1,
            const t_coord_conv& coord2,
            t_coord& new_coord0,
            t_coord& new_coord1,
            t_coord& new_coord2
        )
        {
            UVLM::Types::Vector3 b, abc;
            Eigen::Matrix<UVLM::Types::Real, 3, 3> Amat, Ainv;
            uint i_conv=0;
            for (unsigned int i_m=0; i_m<M; ++i_m)
            {
                while ((dist_to_orig_conv(i_conv) <= dist_to_orig(i_m)) and (i_conv < M))
                {i_conv++;}

                if (i_conv == 1)
                {
                    Amat << dist_to_orig_conv(i_conv - 1)*dist_to_orig_conv(i_conv - 1), dist_to_orig_conv(i_conv - 1), 1.,
                            dist_to_orig_conv(i_conv    )*dist_to_orig_conv(i_conv    ), dist_to_orig_conv(i_conv    ), 1.,
                            dist_to_orig_conv(i_conv + 1)*dist_to_orig_conv(i_conv + 1), dist_to_orig_conv(i_conv + 1), 1.;
                    Ainv = Amat.inverse();

                    b << coord0(i_conv - 1), coord0(i_conv), coord0(i_conv + 1);
                    abc = Ainv*b;
                    new_coord0(i_m) = abc(0)*dist_to_orig(i_m)*dist_to_orig(i_m) +
                                      abc(1)*dist_to_orig(i_m) +
                                      abc(2);

                    b << coord1(i_conv - 1), coord1(i_conv), coord1(i_conv + 1);
                    abc = Ainv*b;
                    new_coord1(i_m) = abc(0)*dist_to_orig(i_m)*dist_to_orig(i_m) +
                                      abc(1)*dist_to_orig(i_m) +
                                      abc(2);

                    b << coord2(i_conv - 1), coord2(i_conv), coord2(i_conv + 1);
                    abc = Ainv*b;
                    new_coord2(i_m) = abc(0)*dist_to_orig(i_m)*dist_to_orig(i_m) +
                                      abc(1)*dist_to_orig(i_m) +
                                      abc(2);
                } else {
                    Amat << dist_to_orig_conv(i_conv - 2)*dist_to_orig_conv(i_conv - 2), dist_to_orig_conv(i_conv - 2), 1.,
                            dist_to_orig_conv(i_conv - 1)*dist_to_orig_conv(i_conv - 1), dist_to_orig_conv(i_conv - 1), 1.,
                            dist_to_orig_conv(i_conv    )*dist_to_orig_conv(i_conv    ), dist_to_orig_conv(i_conv    ), 1.;
                    Ainv = Amat.inverse();

                    b << coord0(i_conv - 2), coord0(i_conv - 1), coord0(i_conv);
                    abc = Ainv*b;
                    new_coord0(i_m) = abc(0)*dist_to_orig(i_m)*dist_to_orig(i_m) +
                                      abc(1)*dist_to_orig(i_m) +
                                      abc(2);

                    b << coord1(i_conv - 2), coord1(i_conv - 1), coord1(i_conv);
                    abc = Ainv*b;
                    new_coord1(i_m) = abc(0)*dist_to_orig(i_m)*dist_to_orig(i_m) +
                                      abc(1)*dist_to_orig(i_m) +
                                      abc(2);

                    b << coord2(i_conv - 2), coord2(i_conv - 1), coord2(i_conv);
                    abc = Ainv*b;
                    new_coord2(i_m) = abc(0)*dist_to_orig(i_m)*dist_to_orig(i_m) +
                                      abc(1)*dist_to_orig(i_m) +
                                      abc(2);
                }
            }
        } // parabolic

        template <typename t_dist,
                  typename t_dist_conv,
                  typename t_coord,
                  typename t_coord_conv>
        void splines
        (
            uint M,
            const t_dist& dist_to_orig,
            const t_dist_conv& dist_to_orig_conv,
            const t_coord_conv& coord0,
            const t_coord_conv& coord1,
            const t_coord_conv& coord2,
            t_coord& new_coord0,
            t_coord& new_coord1,
            t_coord& new_coord2
        )
        {

            const unsigned int splines_degree=4;

            Eigen::Spline<UVLM::Types::Real, 1, splines_degree> spline0 = Eigen::SplineFitting<Eigen::Spline<UVLM::Types::Real, 1, splines_degree>>::Interpolate(coord0, splines_degree, dist_to_orig_conv);

            Eigen::Spline<UVLM::Types::Real, 1, splines_degree> spline1 = Eigen::SplineFitting<Eigen::Spline<UVLM::Types::Real, 1, splines_degree>>::Interpolate(coord1, splines_degree, dist_to_orig_conv);

            Eigen::Spline<UVLM::Types::Real, 1, splines_degree> spline2 = Eigen::SplineFitting<Eigen::Spline<UVLM::Types::Real, 1, splines_degree>>::Interpolate(coord2, splines_degree, dist_to_orig_conv);

            for (uint i_m=0; i_m<M; ++i_m)
            {
                new_coord0(i_m) = spline0(dist_to_orig(i_m))(0);
                new_coord1(i_m) = spline1(dist_to_orig(i_m))(0);
                new_coord2(i_m) = spline2(dist_to_orig(i_m))(0);
            }
        } // splines

        template <typename t_dist,
                  typename t_dist_conv,
                  typename t_coord,
                  typename t_coord_conv>
        void slerp_z
        (
            uint M,
            const t_dist& dist_to_orig,
            const t_dist_conv& dist_to_orig_conv,
            const t_coord_conv& coord0,
            const t_coord_conv& coord1,
            const t_coord_conv& coord2,
            t_coord& new_coord0,
            t_coord& new_coord1,
            t_coord& new_coord2
        )
        {
            // https://en.wikipedia.org/wiki/Slerp
            UVLM::Types::Real to_prev, to_next, prev_to_next, omega, coef_prev, coef_next, mod_next, mod_prev;
            uint i_conv=0;
            for (unsigned int i_m=0; i_m<M; ++i_m)
            {
                new_coord0(i_m) = coord0(i_m);
                new_coord1(i_m) = coord1(i_m);
                new_coord2(i_m) = coord2(i_m);

                while ((dist_to_orig_conv(i_conv) <= dist_to_orig(i_m)) and (i_conv < M))
                {i_conv++;}

                to_prev = dist_to_orig(i_m) - dist_to_orig_conv(i_conv - 1);
                to_next = dist_to_orig_conv(i_conv) - dist_to_orig(i_m);
                prev_to_next = dist_to_orig_conv(i_conv) - dist_to_orig_conv(i_conv - 1);

                mod_prev = std::sqrt(coord0(i_conv - 1)*coord0(i_conv - 1) +
                                     coord1(i_conv - 1)*coord1(i_conv - 1));
                                     // coord2(i_conv - 1)*coord2(i_conv - 1));
                mod_next = std::sqrt(coord0(i_conv)*coord0(i_conv) +
                                     coord1(i_conv)*coord1(i_conv));
                                     // coord2(i_conv)*coord2(i_conv));

                omega = std::acos((coord0(i_conv - 1)*coord0(i_conv) +
                                  coord1(i_conv - 1)*coord1(i_conv))/mod_prev/mod_next);
                                  // coord2(i_conv - 1)*coord2(i_conv))/mod_prev/mod_next);


                if (std::abs(std::sin(omega)) > 1e-6)
                {
                    coef_prev = std::sin(to_next*omega/prev_to_next)/std::sin(omega);
                    coef_next = std::sin(to_prev*omega/prev_to_next)/std::sin(omega);

                    new_coord0(i_m) = (coef_next*coord0(i_conv) + coef_prev*coord0(i_conv - 1));
                    new_coord1(i_m) = (coef_next*coord1(i_conv) + coef_prev*coord1(i_conv - 1));
                } else {
                    new_coord0(i_m) = (to_prev*coord0(i_conv) + to_next*coord0(i_conv - 1))/prev_to_next;
                    new_coord1(i_m) = (to_prev*coord1(i_conv) + to_next*coord1(i_conv - 1))/prev_to_next;
                }

                new_coord2(i_m) = (to_prev*coord2(i_conv) + to_next*coord2(i_conv - 1))/prev_to_next;
            }
        } // slerp_z

        template <typename t_dist,
                  typename t_dist_conv,
                  typename t_coord,
                  typename t_coord_conv>
        void slerp_yaw
        (
            uint M,
            const UVLM::Types::Real yaw,
            const t_dist& dist_to_orig,
            const t_dist_conv& dist_to_orig_conv,
            const t_coord_conv& coord0,
            const t_coord_conv& coord1,
            const t_coord_conv& coord2,
            t_coord& new_coord0,
            t_coord& new_coord1,
            t_coord& new_coord2
        )
        {
            // https://en.wikipedia.org/wiki/Slerp
            // This function computes the slerp interpolation around an axis on the y-z plane rotated 'yaw' degrees around x
            UVLM::Types::Real to_prev, to_next, prev_to_next, omega, coef_prev, coef_next, mod_next, mod_prev;
            UVLM::Types::VectorX aux_coord0, aux_coord1, aux_coord2;
            UVLM::Types::VectorX aux_new_coord0, aux_new_coord1, aux_new_coord2;

            aux_coord0.resize(M + 1);
            aux_coord1.resize(M + 1);
            aux_coord2.resize(M + 1);

            aux_new_coord0.resize(M);
            aux_new_coord1.resize(M);
            aux_new_coord2.resize(M);

            // Transform the coordinates
            for (unsigned int i_m=0; i_m<M+1; ++i_m)
            {
                aux_coord0(i_m) = coord0(i_m);
                aux_coord1(i_m) = coord1(i_m)*cos(yaw) + coord2(i_m)*sin(yaw);
                aux_coord2(i_m) = -1.0*coord1(i_m)*sin(yaw) + coord2(i_m)*cos(yaw);
            }

            // Compute the new coordinates in the yaw FoR
            uint i_conv=0;
            for (unsigned int i_m=0; i_m<M; ++i_m)
            {
                while ((dist_to_orig_conv(i_conv) <= dist_to_orig(i_m)) and (i_conv < M))
                {i_conv++;}

                to_prev = dist_to_orig(i_m) - dist_to_orig_conv(i_conv - 1);
                to_next = dist_to_orig_conv(i_conv) - dist_to_orig(i_m);
                prev_to_next = dist_to_orig_conv(i_conv) - dist_to_orig_conv(i_conv - 1);

                mod_prev = std::sqrt(aux_coord0(i_conv - 1)*aux_coord0(i_conv - 1) +
                                     aux_coord1(i_conv - 1)*aux_coord1(i_conv - 1));
                mod_next = std::sqrt(aux_coord0(i_conv)*aux_coord0(i_conv) +
                                     aux_coord1(i_conv)*aux_coord1(i_conv));

                omega = std::acos((aux_coord0(i_conv - 1)*aux_coord0(i_conv) +
                                  aux_coord1(i_conv - 1)*aux_coord1(i_conv))/mod_prev/mod_next);

                if (std::abs(std::sin(omega)) > 1e-6)
                {
                    coef_prev = std::sin(to_next*omega/prev_to_next)/std::sin(omega);
                    coef_next = std::sin(to_prev*omega/prev_to_next)/std::sin(omega);

                    aux_new_coord0(i_m) = (coef_next*aux_coord0(i_conv) + coef_prev*aux_coord0(i_conv - 1));
                    aux_new_coord1(i_m) = (coef_next*aux_coord1(i_conv) + coef_prev*aux_coord1(i_conv - 1));
                } else {
                    aux_new_coord0(i_m) = (to_prev*aux_coord0(i_conv) + to_next*aux_coord0(i_conv - 1))/prev_to_next;
                    aux_new_coord1(i_m) = (to_prev*aux_coord1(i_conv) + to_next*aux_coord1(i_conv - 1))/prev_to_next;
                }

                aux_new_coord2(i_m) = (to_prev*aux_coord2(i_conv) + to_next*aux_coord2(i_conv - 1))/prev_to_next;
            }

            // Transform back the coordinates
            for (unsigned int i_m=0; i_m<M; ++i_m)
            {
                new_coord0(i_m) = aux_new_coord0(i_m);
                new_coord1(i_m) = aux_new_coord1(i_m)*cos(yaw) - aux_new_coord2(i_m)*sin(yaw);
                new_coord2(i_m) = aux_new_coord1(i_m)*sin(yaw) + aux_new_coord2(i_m)*cos(yaw);
            }

        } // slerp_yaw
    } // Interpolation
    //
    namespace Filters
    {
        template <typename t_coord>
        void moving_average
        (
            uint M,
            const unsigned int window,
            const t_coord& x,
            t_coord& coord0,
            t_coord& coord1,
            t_coord& coord2
        )
        {
            unsigned int sp;
            UVLM::Types::Real aux_coord0[M], aux_coord1[M], aux_coord2[M];

            if (window % 2 == 1)
            {
                // window has to be odd
                sp = int((window - 1)/2);
            } else
            {
                std::cerr << "window has to be odd" << std::endl;
            }

            // Copy values
            for (uint i_m = 0; i_m < M; i_m++)
            {
                aux_coord0[i_m] = coord0(i_m);
                aux_coord1[i_m] = coord1(i_m);
                aux_coord2[i_m] = coord2(i_m);
            }

            for (uint i_m = 1; i_m < M; i_m++)
            {
                if (i_m < sp)
                {
                    // First points
                    for (uint i_avg = 0; i_avg < i_m; i_avg++)
                    {
                        // At coord0(i_m) I already have the value at that point
                        coord0(i_m) += aux_coord0[i_m - i_avg - 1] + coord0[i_m + i_avg + 1];
                        coord1(i_m) += aux_coord1[i_m - i_avg - 1] + coord1[i_m + i_avg + 1];
                        coord2(i_m) += aux_coord2[i_m - i_avg - 1] + coord2[i_m + i_avg + 1];
                    }
                    coord0(i_m) /= (2*i_m + 1);
                    coord1(i_m) /= (2*i_m + 1);
                    coord2(i_m) /= (2*i_m + 1);

                } else if (i_m < M - sp)
                {
                    // Intermediate points
                    for (uint i_avg = 0; i_avg < sp; i_avg++)
                    {
                        // At coord0(i_m) I already have the value at that point
                        coord0(i_m) += aux_coord0[i_m - i_avg - 1] + coord0[i_m + i_avg + 1];
                        coord1(i_m) += aux_coord1[i_m - i_avg - 1] + coord1[i_m + i_avg + 1];
                        coord2(i_m) += aux_coord2[i_m - i_avg - 1] + coord2[i_m + i_avg + 1];
                    }
                    coord0(i_m) /= window;
                    coord1(i_m) /= window;
                    coord2(i_m) /= window;
                } else
                {
                    // Last points
                    for (uint i_avg = 0; i_avg < (M - 1 - i_m); i_avg++)
                    {
                        // At coord0(i_m) I already have the value at that point
                        coord0(i_m) += aux_coord0[i_m - i_avg - 1] + coord0[i_m + i_avg + 1];
                        coord1(i_m) += aux_coord1[i_m - i_avg - 1] + coord1[i_m + i_avg + 1];
                        coord2(i_m) += aux_coord2[i_m - i_avg - 1] + coord2[i_m + i_avg + 1];
                    }
                    coord0(i_m) /= (2*(M - 1 - i_m) + 1);
                    coord1(i_m) /= (2*(M - 1 - i_m) + 1);
                    coord2(i_m) /= (2*(M - 1 - i_m) + 1);
                }
            }
        } // moving_average
    } // Filters

} // UVLM
