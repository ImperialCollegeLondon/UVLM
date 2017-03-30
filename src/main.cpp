#include "EigenInclude.h"
#include "solver.h"
#include "constants.h"

#include <iostream>
#include <vector>

int main()
{
    // create input data
    //TODO M, N and Mstar independant for different surfaces
    UVLM::Types::VMopts VMOPTS;
    VMOPTS.N = 2;
    VMOPTS.M = 2;
    VMOPTS.ImageMethod = false;
    VMOPTS.Mstar = 1;
    VMOPTS.NumSurfaces = 1;


    int M = VMOPTS.M;
    int N = VMOPTS.N;
    int Mstar = VMOPTS.Mstar;
    int n_surf = VMOPTS.NumSurfaces;

    UVLM::Types::Real u = 1.0;

    UVLM::Types::VecVecMatrixX zeta;
    UVLM::Types::allocate_VecVecMat(zeta,
                                    n_surf,
                                    UVLM::Constants::NDIM,
                                    M+1, N+1);
    UVLM::Types::VecVecMatrixX zeta_dot;
    UVLM::Types::allocate_VecVecMat(zeta_dot,
                                    n_surf,
                                    UVLM::Constants::NDIM,
                                    M+1, N+1);
    UVLM::Types::VecVecMatrixX u_ext;
    UVLM::Types::allocate_VecVecMat(u_ext,
                                    n_surf,
                                    UVLM::Constants::NDIM,
                                    M+1, N+1);
    UVLM::Types::VecVecMatrixX zeta_star;
    UVLM::Types::allocate_VecVecMat(zeta_star,
                                    n_surf,
                                    UVLM::Constants::NDIM,
                                    Mstar+1, N+1);

    for (int i_surf=0; i_surf<n_surf; ++i_surf)
    {
        for (int i_dim=0; i_dim<UVLM::Constants::NDIM; ++i_dim)
        {
            if (i_dim == 0)
            {
                for (int iM=0; iM<M+1; ++iM)
                {
                    u_ext[i_surf][i_dim].row(iM) << u, u, u;
                    zeta[i_surf][i_dim].row(iM) << 0.0, 0.5, 1.0;
                }
                for (int iMstar=0; iMstar<Mstar+1; ++iMstar)
                {
                    zeta_star[i_surf][i_dim].row(iMstar) << 0.0, 0.5, 1.0;
                }
            } else if (i_dim == 1)
            {
                for (int iN=0; iN<N+1; ++iN)
                {

                    zeta[i_surf][i_dim].col(iN) << 0.0, 0.5, 1.0;
                    zeta_star[i_surf][i_dim].col(iN) << 1.0, 1.5;
                }
            }
        }
        // std::cout << "zeta x" << std::endl;
        // std::cout << zeta[i_surf][0] << std::endl;
        // std::cout << "---" << std::endl;
        // std::cout << "zeta y" << std::endl;
        // std::cout << zeta[i_surf][1] << std::endl;
        // std::cout << "---" << std::endl;
        // std::cout << "zeta z" << std::endl;
        // std::cout << zeta[i_surf][2] << std::endl;
        // std::cout << "---" << std::endl;
    }

    std::cout << "calling solver" << std::endl;
    UVLM::Solver::solve(zeta, zeta_dot, u_ext, zeta_star, VMOPTS);


}
