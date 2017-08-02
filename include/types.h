#pragma once

#include "EigenInclude.h"
#include "constants.h"
#include <vector>
#include <utility>

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
        typedef std::vector<MapMatrixX> VecMapX;
        typedef std::vector<VecMapX> VecVecMapX;

        typedef Eigen::DenseBase<Real> DenseBase;
        typedef Eigen::Block<MatrixX> Block;

        // Vectors
        typedef Eigen::Matrix<Real, 3, 1> Vector3;
        typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorX;

        // std custom containers
        typedef std::pair<unsigned int, unsigned int> IntPair;
        typedef std::vector<IntPair> VecDimensions;


        struct VMopts {
        	bool ImageMethod;
        	unsigned int Mstar;
        	bool Steady;
            bool horseshoe;
        	bool KJMeth;
        	bool NewAIC;
        	double DelTime;
        	bool Rollup;
        	unsigned int NumCores;
        	unsigned int NumSurfaces;
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
                    mat[i_surf].push_back(UVLM::Types::MatrixX());
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
            uint n_surf = in.size();
            for (uint i_surf=0; i_surf<n_surf; ++i_surf)
            {
                uint n_dim = in[i_surf].size();
                for (uint i_dim=0; i_dim<n_dim; ++i_dim)
                {
                    out[i_surf][i_dim] = in[i_surf][i_dim];
                }
            }
        }
    }
}


#include "typeutils.h"
