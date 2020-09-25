#pragma once

#include "EigenInclude.h"
#include <vector>
#include <utility>
#include <iostream>

// convenience declarations
typedef unsigned int uint;

namespace UVLM
{
    namespace Types
    {
        // Working precision
        typedef double Real;

        // Eigen shortcuts
        // Matrices
        typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixX;
        typedef Eigen::Map<MatrixX> MapMatrixX;
        typedef Eigen::Map<const MatrixX> cMapMatrixX;
        typedef std::vector<MatrixX> VecMatrixX;
        typedef std::vector<VecMatrixX> VecVecMatrixX;
        typedef std::vector<VecVecMatrixX> VecVecVecMatrixX;
        typedef std::vector<MapMatrixX> VecMapX;
        typedef std::vector<VecMapX> VecVecMapX;
        typedef std::vector<VecVecMapX> VecVecVecMapX;

        typedef Eigen::DenseBase<Real> DenseBase;
        typedef Eigen::Block<MatrixX> Block;

        // Vectors
        typedef Eigen::Matrix<Real, 3, 1> Vector3;
        typedef Eigen::Matrix<Real, 6, 1> Vector6;
        typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorX;
        typedef Eigen::Map<VectorX> MapVectorX;
        typedef std::vector<MapVectorX> VecMapVX;

        // std custom containers
        typedef std::pair<unsigned int, unsigned int> IntPair;
        typedef std::vector<IntPair> VecDimensions;


        struct VMopts
        {
        	bool ImageMethod;
        	// unsigned int Mstar;
        	bool Steady;
            bool horseshoe;
        	bool KJMeth;
        	bool NewAIC;
        	double DelTime;
        	bool Rollup;
        	unsigned int NumCores;
        	unsigned int NumSurfaces;
            double dt;
            unsigned int n_rollup;
            double rollup_tolerance;
            unsigned int rollup_aic_refresh;
            bool iterative_solver;
            double iterative_tol;
            bool iterative_precond;
            bool cfl1;
            double vortex_radius;
            double vortex_radius_wake_ind;
        };

        struct UVMopts
        {
            double dt;
            uint NumCores;
            uint NumSurfaces;
            // uint steady_n_rollup;
            // uint steady_rollup_tolerance;
            // uint steady_rollup_aic_refresh;
            uint convection_scheme;
            // uint Mstar;
            bool ImageMethod;
            bool iterative_solver;
            double iterative_tol;
            bool iterative_precond;
            bool convect_wake;
            bool cfl1;
            double vortex_radius;
            double vortex_radius_wake_ind;
            uint interp_coords;
            uint filter_method;
            uint interp_method;
            double yaw_slerp;
            bool quasi_steady;
        };

        VMopts UVMopts2VMopts(const UVMopts& uvm)
        {
            VMopts vm;
            vm.dt = uvm.dt;
            vm.NumCores = uvm.NumCores;
            vm.NumSurfaces = uvm.NumSurfaces;
            // vm.Mstar = uvm.Mstar;
            vm.ImageMethod = uvm.ImageMethod;
            vm.iterative_solver = uvm.iterative_solver;
            vm.iterative_tol = uvm.iterative_tol;
            vm.iterative_precond = uvm.iterative_precond;
            vm.horseshoe = false;
            vm.vortex_radius = uvm.vortex_radius;
            vm.vortex_radius_wake_ind = uvm.vortex_radius_wake_ind;
            if (uvm.quasi_steady)
            {
                vm.Steady = true;
            } else {
                vm.Steady = false;
            }
            return vm;
        };

        struct FlightConditions
        {
            double uinf = 1.0;
            double uinf_direction[3];
            double rho = 1.225;
            double c_ref = 1.0;
        };

        template <typename t_mat>
        inline void generate_dimensions
        (
            const t_mat& mat,
            UVLM::Types::VecDimensions& dimensions,
            const int& correction = 0
        )
        {
            dimensions.resize(mat.size());
            for (unsigned int i_surf=0; i_surf<dimensions.size(); ++i_surf)
            {
                dimensions[i_surf] = UVLM::Types::IntPair
                                                    (
                                                        mat[i_surf][0].rows() + correction,
                                                        mat[i_surf][0].cols() + correction
                                                    );
            }
        }


        inline void allocate_VecMat
        (
            UVLM::Types::VecMatrixX& mat,
            const UVLM::Types::VecDimensions& dimensions,
            const int& correction = 0,
            const UVLM::Types::Real& initial_value = 0.0
        )
        {
            unsigned int n_surf = dimensions.size();
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                int M = dimensions[i_surf].first + correction;
                int N = dimensions[i_surf].second + correction;
                mat.push_back(UVLM::Types::MatrixX());
                if (initial_value == 0.0)
                {
                    mat[i_surf].setZero(M,N);
                } else if (initial_value == 1.0)
                {
                    mat[i_surf].setOnes(M,N);
                } else
                {
                    mat[i_surf].setConstant(M,N,initial_value);
                }


            }
        }

        template <typename t_dimensions_in>
        inline void allocate_VecMat
        (
            UVLM::Types::VecMatrixX& mat,
            const t_dimensions_in& dimensions_in,
            const int& correction = 0,
            const double& initial_value = 0.0
        )
        {
            const unsigned int n_mats = dimensions_in.size();
            mat.resize(n_mats);
            for (unsigned int i=0; i<n_mats; ++i)
            {
                mat[i].setConstant
                (
                    dimensions_in[i].rows(),
                    dimensions_in[i].cols(),
                    initial_value
                );
            }
        }

        inline void allocate_VecMat
        (
             UVLM::Types::VecMatrixX& mat,
             const unsigned int& n_surf,
             const unsigned int& M,
             const unsigned int& N
        )
        {
             mat.resize(n_surf);
             for (auto& surf: mat)
             {
                 surf.resize(M, N);
                 surf.setZero(M, N);
             }
        }

        inline void initialise_VecMat
        (
            UVLM::Types::VecMatrixX& mat,
            const double& value = 0.0
        )
        {
            const unsigned int n_mats = mat.size();
            for (unsigned int i=0; i<n_mats; ++i)
            {
                mat[i].setConstant(
                    mat[i].rows(),
                    mat[i].cols(),
                    value);
            }
        }

        template <typename t_mat>
        inline void initialise_VecVecMat
        (
            t_mat& mat,
            const double& value = 0.0
        )
        {
            const unsigned int n_surf = mat.size();
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const uint n_dim = mat[i_surf].size();
                for (unsigned int i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    mat[i_surf][i_dim].setConstant(value);
                }
            }
        }

        inline void allocate_VecVecMat
        (
            UVLM::Types::VecVecMatrixX& mat,
            const unsigned int& n_surf,
            const unsigned int& n_dim,
            const unsigned int& M,
            const unsigned int& N
        )
        {
            mat.resize(n_surf);
            for (auto& surf: mat)
            {
                surf.resize(n_dim);
                for (unsigned int i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    surf.push_back(UVLM::Types::MatrixX());
                    surf[i_dim].resize(M, N);
                    surf[i_dim].setZero(M, N);
                }
            }
        }

        inline void allocate_VecVecMat
        (
            UVLM::Types::VecVecMatrixX& mat,
            const unsigned int& n_dim,
            const UVLM::Types::VecDimensions& dimensions,
            const int& correction = 0
        )
        {
            unsigned int n_surf = dimensions.size();
            mat.resize(n_surf);
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                int M = dimensions[i_surf].first + correction;
                int N = dimensions[i_surf].second + correction;
                mat[i_surf].resize(n_dim);
                for (unsigned int i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    mat[i_surf].push_back(UVLM::Types::MatrixX());
                    mat[i_surf][i_dim].resize(M, N);
                    mat[i_surf][i_dim].setZero(M, N);
                }
            }
        }

        template <typename t_dimensions>
        inline void allocate_VecVecMat
        (
            UVLM::Types::VecVecMatrixX& mat,
            const t_dimensions& in_dimensions,
            const int& correction = 0
        )
        {
            unsigned int n_surf = in_dimensions.size();
            mat.resize(n_surf);
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                const unsigned int M = in_dimensions[i_surf][0].rows() + correction;
                const unsigned int N = in_dimensions[i_surf][0].cols() + correction;
                mat[i_surf].resize(in_dimensions[i_surf].size());
                for (unsigned int i_dim=0; i_dim<in_dimensions[i_surf].size(); ++i_dim)
                {
                    mat[i_surf][i_dim].resize(M, N);
                    mat[i_surf][i_dim].setZero(M, N);
                }
            }
        }

        template <typename t_in,
                  typename t_out>
        inline void copy_VecVecMat
        (
            const t_in& in,
            t_out& out
        )
        {
            uint M, N;
            uint n_surf = in.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                uint n_dim = in[i_surf].size();
                M = in[i_surf][0].rows();
                N = in[i_surf][0].cols();
                for (uint i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    for (uint i_m=0; i_m<M; ++i_m)
                    {
                        for (uint i_n=0; i_n<N; ++i_n)
                        {
                            out[i_surf][i_dim](i_m, i_n) = in[i_surf][i_dim](i_m, i_n);
                        }
                    }
                }
            }
        }

        template <typename t_mat>
        inline double norm_VecVec_mat
        (
            const t_mat& mat
        )
        {
            double norm = 0.0;
            uint n_surf = mat.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                uint n_dim = mat[i_surf].size();
                for (uint i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    norm += mat[i_surf][i_dim].norm();
                }
            }
            return norm;
        }

        template <typename t_mat>
        inline double max_VecVecMat
        (
            const t_mat& mat
        )
        {
            double max = 0.0;
            uint n_surf = mat.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                uint n_dim = mat[i_surf].size();
                for (uint i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    max = std::max(max, std::abs(mat[i_surf][i_dim].maxCoeff()));
                }
            }
            return max;
        }

        UVLM::Types::Vector3 zeroVector3()
        {
            UVLM::Types::Vector3 vec;
            vec.setZero();
            return vec;
        }

        // inline void allocate_VecVecVecMat
        // (
        //     UVLM::Types::VecVecVecMatrixX& mat,
        //     const int& n_tsteps,
        //     const UVLM::Types::VecDimensions& dimensions,
        //     const int& correction = 0
        // )
        // {
        //     mat.resize(n_tsteps);
        //     UVLM::Types::allocate_VecVecMat(mat,
        //                                     UVLM::Constants::NDIM,
        //                                     dimensions,
        //                                     correction);
        // }

        // template <typename t_in_dimensions>
        // inline void allocate_VecVecVecMat
        // (
        //     UVLM::Types::VecVecVecMatrixX& mat,
        //     const t_in_dimensions& in_dimensions,
        //     const int& correction = 0
        // )
        // {
        //     mat.resize(n_tsteps);
        //     UVLM::Types::allocate_VecVecMat(mat,
        //                                     in_dimensions,
        //                                     correction);
        // }
    }
}



/*
Define types for C++ interface on linear UVLM routines.
Matrices size is specified whenever possible to maximise speed.
*/

namespace UVLMlin{

    using namespace Eigen;

    typedef Matrix<double,4,3,RowMajor> Matrix4by3d;

    // map 1d arrays into Eigen Matrices (interface for 1D or 2D python arrays)
    typedef Map< Matrix<double,Dynamic,Dynamic,RowMajor> > map_Mat;
    typedef Map< Matrix<double,4,3,RowMajor> > map_Mat4by3;
    typedef Map< Matrix<double,3,3,RowMajor> > map_Mat3by3;
    typedef Map< Matrix<double,1,3> > map_RowVec3;

    // mapp 3D python arrays
    typedef std::vector<map_Mat3by3> Vec_map_Mat3by3;
    typedef std::vector<map_Mat> Vec_map_Mat;
}



#include "typeutils.h"
