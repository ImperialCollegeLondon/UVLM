#pragma once
#include "EigenInclude.h"
#include "types.h"
#include "constants.h"
#include "geometry.h"
#include "steady.h"
#include "unsteady.h"

#include <iostream>
//#ifdef _OPENMP
    //#include "omp.h"
//#endif

#define DLLEXPORT extern "C"
// #pragma warning disable 1017

namespace UVLM
{
    namespace CppInterface
    {
        void map_VecVecMat(const UVLM::Types::VecDimensions& dimensions,
                           double** in,
                           UVLM::Types::VecVecMapX& map,
                           const int& correction=0,
                           const unsigned int& n_dim=UVLM::Constants::NDIM)
        {
            const unsigned int n_surf = dimensions.size();
            map.resize(n_surf);
            unsigned int counter = 0;
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                for (uint i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    map[i_surf].push_back(UVLM::Types::MapMatrixX (in[counter],
                                                                   dimensions[i_surf].first + correction,
                                                                   dimensions[i_surf].second + correction));
                    counter++;
                }
            }
        }

        void map_VecMat(const UVLM::Types::VecDimensions& dimensions,
                        double** in,
                        UVLM::Types::VecMapX& map,
                        const int& correction=0)
        {
            const unsigned int n_surf = dimensions.size();
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                map.push_back(UVLM::Types::MapMatrixX (in[i_surf],
                                                       dimensions[i_surf].first + correction,
                                                       dimensions[i_surf].second + correction));
            }
        }

        void map_VecVec1(const UVLM::Types::VecDimensions& dimensions,
                        double** in,
                        UVLM::Types::VecMapVX& map,
                        //UVLM::Types::MapVectorX& map,
                        const int& correction=0)
        {
            // Generates a variable that will be indexed as map[i_surf](i_m)
            const unsigned int n_surf = dimensions.size();
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                map.push_back(UVLM::Types::MapVectorX (in[i_surf],
                                                       dimensions[i_surf].first + correction));
                                                      // dimensions[i_surf].second + correction));
            }
        }

        void transform_dimensions(unsigned int& n_surf,
                                  unsigned int** dimensions_in,
                                  UVLM::Types::VecDimensions& dimensions)
        {
            dimensions.resize(n_surf);
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                dimensions[i_surf].first = dimensions_in[i_surf][0];
                dimensions[i_surf].second= dimensions_in[i_surf][1];
            }
        }
    }
}


namespace UVLMlin{

  extern "C" void call_der_biot_panel(double p_DerP[9],
                    double p_DerVertices[36],
                    double p_zetaP[3],
                    double p_ZetaPanel[12],
                    const double& gamma,
                    double& vortex_radius);


  extern "C" void call_biot_panel(double p_vel[3],
                  double p_zetaP[3],
                  double p_ZetaPanel[12],
                  const double& gamma,
                  double& vortex_radius);
}
