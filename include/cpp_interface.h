#include "EigenInclude.h"
#include "constants.h"
#include "types.h"

#include <iostream>


namespace UVLM
{
    namespace CppInterface
    {
        void map_VecVecMat(const UVLM::Types::VecDimensions& dimensions,
                           double** in,
                           UVLM::Types::VecVecMapX& map)
        {
            const unsigned int n_surf = dimensions.size();
            map.resize(n_surf);
            unsigned int counter = -1;
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                // map[i_surf].resize(UVLM::Constants::NDIM);
                for (unsigned int i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
                {
                    map[i_surf].push_back(UVLM::Types::MapMatrixX (in[++counter],
                                                                   dimensions[i_surf].first,
                                                                   dimensions[i_surf].second));
                }
            }
        }

        void transform_dimensions(int& n_surf,
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
