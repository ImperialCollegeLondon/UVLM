#include "EigenInclude.h"
#include "constants.h"
#include "types.h"
#include "geometry.h"

#include <iostream>

#define DLLEXPORT extern "C"

namespace UVLM
{
    namespace CppInterface
    {
        void map_VecVecMat(const UVLM::Types::VecDimensions& dimensions,
                           double** in,
                           UVLM::Types::VecVecMapX& map,
                           const int& correction=0,
                           const int& n_dim=UVLM::Constants::NDIM)
        {
            const unsigned int n_surf = dimensions.size();
            map.resize(n_surf);
            unsigned int counter = 0;
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                for (unsigned int i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    map[i_surf].push_back(UVLM::Types::MapMatrixX (in[counter],
                                                                   dimensions[i_surf].first + correction,
                                                                   dimensions[i_surf].second + correction));
                    counter++;
                }
            }
        }

        void transform_dimensions(unsigned int& n_surf,
                                  int** dimensions_in,
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
