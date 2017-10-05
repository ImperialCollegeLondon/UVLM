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

    UVLM::Steady::solver(zeta,
                         zeta_dot,
                         u_ext,
                         zeta_star,
                         gamma,
                         gamma_star,
                         forces,
                         options,
                         flightconditions);
}


DLLEXPORT void init_UVLM
(
    const UVLM::Types::VMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    unsigned int** p_dimensions,
    unsigned int** p_dimensions_star,
    double** p_uext,
    double** p_zeta,
    double** p_zeta_star,
    double** p_zeta_dot,
    double** p_zeta_star_dot,
    double*  p_rbm_vel,
    double** p_gamma,
    double** p_gamma_star,
    double** p_normals,
    double** p_forces
)
{
    // feenableexcept(FE_INVALID | FE_OVERFLOW);
    Eigen::setNbThreads(options.NumCores);
    uint n_surf = options.NumSurfaces;

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

    UVLM::Types::VecVecMapX zeta_dot;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      p_zeta_dot,
                                      zeta_dot,
                                      1);

    UVLM::Types::VecVecMapX zeta_star_dot;
    UVLM::CppInterface::map_VecVecMat(dimensions_star,
                                      p_zeta_star_dot,
                                      zeta_star_dot,
                                      1);

    UVLM::Types::MapVectorX rbm_velocity (p_rbm_vel, 2*UVLM::Constants::NDIM);

    UVLM::Types::VecVecMapX uext;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      p_uext,
                                      uext,
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

    UVLM::Types::VecVecMapX forces;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      p_forces,
                                      forces,
                                      1,
                                      2*UVLM::Constants::NDIM);


    // UVLM::Types::VMopts steady_options = UVLM::Types::UVMopts2VMopts(options);
    UVLM::Unsteady::initialise(
        zeta,
        zeta_dot,
        zeta_star,
        zeta_star_dot,
        uext,
        gamma,
        gamma_star,
        rbm_velocity,
        forces,
        options,
        flightconditions
    );
}


DLLEXPORT void run_UVLM
(
    const UVLM::Types::UVMopts& options,
    const UVLM::Types::FlightConditions& flightconditions,
    unsigned int** p_dimensions,
    unsigned int** p_dimensions_star,
    unsigned int i_iter,
    double** p_uext,
    double** p_uext_star,
    double** p_zeta,
    double** p_zeta_star,
    double** p_zeta_dot,
    double** p_zeta_star_dot,
    double*  p_rbm_vel,
    double** p_gamma,
    double** p_gamma_star,
    double** p_previous_gamma,
    double** p_normals,
    double** p_forces,
    double** p_dynamic_forces
)
{
    Eigen::setNbThreads(options.NumCores);
    uint n_surf = options.NumSurfaces;

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

    UVLM::Types::VecVecMapX zeta_dot;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      p_zeta_dot,
                                      zeta_dot,
                                      1);

    UVLM::Types::VecVecMapX zeta_star_dot;
    UVLM::CppInterface::map_VecVecMat(dimensions_star,
                                      p_zeta_star_dot,
                                      zeta_star_dot,
                                      1);

    UVLM::Types::MapVectorX rbm_velocity (p_rbm_vel, 2*UVLM::Constants::NDIM);

    UVLM::Types::VecVecMapX uext;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      p_uext,
                                      uext,
                                      1);

    UVLM::Types::VecVecMapX uext_star;
    UVLM::CppInterface::map_VecVecMat(dimensions_star,
                                      p_uext_star,
                                      uext_star,
                                      1);

    UVLM::Types::VecMapX gamma;
    UVLM::CppInterface::map_VecMat(dimensions,
                                   p_gamma,
                                   gamma,
                                   0);

    UVLM::Types::VecMapX previous_gamma;
    UVLM::CppInterface::map_VecMat(dimensions,
                                   p_previous_gamma,
                                   previous_gamma,
                                   0);

    UVLM::Types::VecMapX gamma_star;
    UVLM::CppInterface::map_VecMat(dimensions_star,
                                   p_gamma_star,
                                   gamma_star,
                                   0);

    UVLM::Types::VecVecMapX forces;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      p_forces,
                                      forces,
                                      1,
                                      2*UVLM::Constants::NDIM);

    UVLM::Types::VecVecMapX dynamic_forces;
    UVLM::CppInterface::map_VecVecMat(dimensions,
                                      p_dynamic_forces,
                                      dynamic_forces,
                                      1,
                                      2*UVLM::Constants::NDIM);

    UVLM::Unsteady::solver
    (
        i_iter,
        zeta,
        zeta_dot,
        uext,
        uext_star,
        zeta_star,
        gamma,
        gamma_star,
        previous_gamma,
        rbm_velocity,
        forces,
        dynamic_forces,
        options,
        flightconditions
    );
}
