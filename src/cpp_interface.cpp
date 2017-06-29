#include "cpp_interface.h"
#include <fenv.h>


DLLEXPORT void run_VLM
(
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    unsigned int** p_dimensions,
    unsigned int** p_dimensions_star,
    double** p_zeta,
    double** p_zeta_star,
    double** p_u_ext,
    double** p_gamma,
    double** p_gamma_star,
    double** p_normals,
    double** p_forces
)
{
    // feenableexcept(FE_INVALID | FE_OVERFLOW);
    Eigen::setNbThreads(options.NumCores);
    unsigned int n_surf;
    n_surf = options.NumSurfaces;
    UVLM::Types::VecDimensions dimensions;
    UVLM::CppInterface::transform_dimensions(n_surf,
                                             p_dimensions,
                                             dimensions);
    UVLM::Types::VecDimensions dimensions_star;
    UVLM::CppInterface::transform_dimensions(n_surf,
                                             p_dimensions_star,
                                             dimensions_star);

    UVLM::Types::VecVecMapX zeta;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      p_zeta,
                                      zeta,
                                      1);

    UVLM::Types::VecVecMapX zeta_star;
    UVLM::CppInterface::map_VecVecMat(dimensions_star,
                                      p_zeta_star,
                                      zeta_star,
                                      1);

    UVLM::Types::VecVecMapX u_ext;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      p_u_ext,
                                      u_ext,
                                      1);

    UVLM::Types::VecMapX gamma;
    UVLM::CppInterface::map_VecMat(dimensions,
                                   p_gamma,
                                   gamma,
                                   0);

    UVLM::Types::VecMapX gamma_star;
    UVLM::CppInterface::map_VecMat(dimensions_star,
                                   p_gamma_star,
                                   gamma_star,
                                   0);

    UVLM::Types::VecVecMapX normals;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      p_normals,
                                      normals,
                                      0);

    UVLM::Types::VecVecMapX forces;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      p_forces,
                                      forces,
                                      1,
                                      2*UVLM::Constants::NDIM);

    // zeta_dot is zero for VLM simulations
    UVLM::Types::VecVecMatrixX zeta_dot;
    UVLM::Types::allocate_VecVecMat(zeta_dot,
                                    zeta);

    UVLM::Solver::solve(zeta,
                        zeta_dot,
                        u_ext,
                        zeta_star,
                        gamma,
                        gamma_star,
                        normals,
                        forces,
                        options,
                        flightconditions);
}
