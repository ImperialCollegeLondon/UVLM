#pragma once

#include "EigenInclude.h"
#include "types.h"

#include <iostream>

namespace UVLM
{
    namespace Geometry
    {
        template <typename type>
        void panel_normal(type& x,
                          type& y,
                          type& z,
                          UVLM::Types::Vector3& normal
                          )
        {
            UVLM::Types::Vector3 A(x(1,1) - x(0,0),
                                   y(1,1) - y(0,0),
                                   z(1,1) - z(0,0));

            UVLM::Types::Vector3 B(x(1,0) - x(0,1),
                                   y(1,0) - y(0,1),
                                   z(1,0) - z(0,1));

            normal = A.cross(B);
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
                    const unsigned int M = normal[i_surf][i_dim].rows();
                    const unsigned int N = normal[i_surf][i_dim].cols();

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



        
    }
}
