#include "cpp_interface.h"


DLLEXPORT void run_VLM
(
    const UVLM::Types::VMopts& options,
    const int& n_dim,
    int** dimensions_in,
    int** dimensions_star_in,
    double** zeta_in,
    double** zeta_star_in,
    double** normals_in
    // double** u_ext,
    // double** zeta_star,
    // double** gamma,
    // double** gamma_star,
    // UVLM::Types::VMopts& VMOPTS
)
{
    unsigned int n_surf;
    n_surf = options.NumSurfaces;
    UVLM::Types::VecDimensions dimensions;
    UVLM::CppInterface::transform_dimensions(n_surf,
                                             dimensions_in,
                                             dimensions);
    UVLM::Types::VecDimensions dimensions_star;
    UVLM::CppInterface::transform_dimensions(n_surf,
                                             dimensions_star_in,
                                             dimensions_star);

    UVLM::Types::VecVecMapX zeta;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      zeta_in,
                                      zeta,
                                      1);

    UVLM::Types::VecVecMapX zeta_star;
    UVLM::CppInterface::map_VecVecMat(dimensions_star,
                                      zeta_star_in,
                                      zeta_star,
                                      1);

    UVLM::Types::VecVecMapX normals;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      normals_in,
                                      normals);

    UVLM::Geometry::generate_surfaceNormal(zeta, normals);
}
