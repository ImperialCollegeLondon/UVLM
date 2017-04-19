#include "EigenInclude.h"
#include "solver.h"
#include "constants.h"
#include "types.h"
#include "timing.h"

#include <iostream>
#include <vector>
#include <math.h>

int main()
{
    Eigen::initParallel();
    // create input data
    UVLM::Types::VMopts VMOPTS;
    VMOPTS.ImageMethod = false;
    VMOPTS.Mstar = 1;
    VMOPTS.NumSurfaces = 2;
    VMOPTS.Steady = false;

    int Mstar = VMOPTS.Mstar;
    int n_surf = VMOPTS.NumSurfaces;

    unsigned int M = 2;
    unsigned int N = 3;

    UVLM::Types::VecDimensions dimensions;
    dimensions.resize(n_surf);
    dimensions[0] = UVLM::Types::IntPair(M, N);
    dimensions[1] = UVLM::Types::IntPair(M, N);
    UVLM::Types::VecDimensions dimensions_star;
    dimensions_star.resize(n_surf);
    dimensions_star[0] = UVLM::Types::IntPair(VMOPTS.Mstar, N);
    dimensions_star[1] = UVLM::Types::IntPair(VMOPTS.Mstar, N);

    UVLM::Types::Real u = 1.0;
    UVLM::Types::Real alpha = 5; // degrees

    UVLM::Types::VecVecMatrixX zeta;
    UVLM::Types::allocate_VecVecMat(zeta,
                                    UVLM::Constants::NDIM,
                                    dimensions,
                                    1);
    UVLM::Types::VecVecMatrixX zeta_dot;
    UVLM::Types::allocate_VecVecMat(zeta_dot,
                                    UVLM::Constants::NDIM,
                                    dimensions,
                                    1);
    UVLM::Types::VecVecMatrixX u_ext;
    UVLM::Types::allocate_VecVecMat(u_ext,
                                    UVLM::Constants::NDIM,
                                    dimensions,
                                    1);
    UVLM::Types::VecVecMatrixX zeta_star;
    UVLM::Types::allocate_VecVecMat(zeta_star,
                                    UVLM::Constants::NDIM,
                                    dimensions_star,
                                    1);
    UVLM::Types::VecMatrixX gamma;
    UVLM::Types::allocate_VecMat(gamma,
                                 dimensions);
    UVLM::Types::VecMatrixX gamma_star;
    UVLM::Types::allocate_VecMat(gamma_star,
                                 dimensions_star,
                                 0,
                                 1.0);

    // surface 0------------------------------------------------------
    unsigned int i_surf = 0;
    Eigen::VectorXd temp1, temp2;
    temp1.setLinSpaced(dimensions[0].first + 1, 0.0, 1.0);
    temp2.setLinSpaced(dimensions[0].second + 1, 1.0, 1.0);
    zeta[i_surf][0].noalias() = temp1*temp2.transpose();
    // std::cout << zeta[i_surf][0] << std::endl;

    temp1.setLinSpaced(dimensions[0].first + 1, 1.0, 1.0);
    temp2.setLinSpaced(dimensions[0].second + 1, 0.0, 2.0);
    zeta[i_surf][1].noalias() = temp1*temp2.transpose();
    // std::cout << zeta[i_surf][1] << std::endl;

    temp1.setLinSpaced(dimensions_star[0].first + 1, 1.0, 2.0);
    temp2.setLinSpaced(dimensions_star[0].second + 1, 1.0, 1.0);
    zeta_star[i_surf][0].noalias() = temp1*temp2.transpose();
    // std::cout << zeta_star[i_surf][0] << std::endl;

    temp1.setLinSpaced(dimensions_star[0].first + 1, 1.0, 1.0);
    temp2.setLinSpaced(dimensions_star[0].second + 1, 0.0, 2.0);
    zeta_star[i_surf][1].noalias() = temp1*temp2.transpose();
    // std::cout << zeta_star[i_surf][1] << std::endl;

    u_ext[i_surf][0].setConstant(dimensions[0].first + 1, dimensions[0].second + 1,
                                 u*std::cos(alpha*UVLM::Constants::DEGREES2RAD));
    u_ext[i_surf][2].setConstant(dimensions[0].first + 1, dimensions[0].second + 1,
                                 u*std::sin(alpha*UVLM::Constants::DEGREES2RAD));

    // surface 1------------------------------------------------------
    i_surf = 1;
    temp1.setLinSpaced(dimensions[1].first + 1, 0.0, 1.0);
    temp2.setLinSpaced(dimensions[1].second + 1, 1.0, 1.0);
    zeta[i_surf][0].noalias() = temp1*temp2.transpose();
    // std::cout << zeta[i_surf][0] << std::endl;

    temp1.setLinSpaced(dimensions[1].first + 1, 1.0, 1.0);
    temp2.setLinSpaced(dimensions[1].second + 1, -2.0, 0.0);
    zeta[i_surf][1].noalias() = temp1*temp2.transpose();
    // std::cout << zeta[i_surf][1] << std::endl;

    temp1.setLinSpaced(dimensions_star[1].first + 1, 1.0, 2.0);
    temp2.setLinSpaced(dimensions_star[1].second + 1, 1.0, 1.0);
    zeta_star[i_surf][0].noalias() = temp1*temp2.transpose();
    // std::cout << zeta_star[i_surf][1] << std::endl;

    temp1.setLinSpaced(dimensions_star[1].first + 1, 1.0, 1.0);
    temp2.setLinSpaced(dimensions_star[1].second + 1, -2.0, 0.0);
    zeta_star[i_surf][1].noalias() = temp1*temp2.transpose();
    // std::cout << zeta_star[i_surf][1] << std::endl;

    u_ext[i_surf][0].setConstant(dimensions[1].first + 1, dimensions[1].second + 1,
                                 u*std::cos(alpha*UVLM::Constants::DEGREES2RAD));
    u_ext[i_surf][2].setConstant(dimensions[1].first + 1, dimensions[1].second + 1,
                                 u*std::sin(alpha*UVLM::Constants::DEGREES2RAD));

    std::cout << "calling solver" << std::endl;
    UVLM::Timing::Timer time;
    time.tic();
    UVLM::Solver::solve(zeta, zeta_dot, u_ext, zeta_star, gamma, gamma_star, VMOPTS);
    time.toc();

}
