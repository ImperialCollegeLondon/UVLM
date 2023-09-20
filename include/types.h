/**
 * @file types.h
 * @brief Header file containing type definitions and utility functions for the UVLM (Unsteady Vortex Lattice Method) framework.
 */

#pragma once

#include "EigenInclude.h"
#include <vector>
#include <utility>
#include <iostream>

// Convenience declarations
typedef unsigned int uint;

namespace UVLM
{
    namespace Types
    {        
        /**
         * @brief Working precision type.
         */
        typedef double Real;

        // Eigen shortcuts for matrices
        typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixX;
        typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXint;
        typedef Eigen::Map<MatrixX> MapMatrixX;
        typedef Eigen::Map<MatrixXint> MapMatrixXint;
        typedef Eigen::Map<const MatrixX> cMapMatrixX;
        typedef std::vector<MatrixX> VecMatrixX;;
        typedef std::vector<VecMatrixX> VecVecMatrixX;
        typedef std::vector<VecVecMatrixX> VecVecVecMatrixX;
        typedef std::vector<MapMatrixX> VecMapX;
        typedef std::vector<MapMatrixXint> VecMapXint;
        typedef std::vector<VecMapX> VecVecMapX;
        typedef std::vector<VecVecMapX> VecVecVecMapX;
        
        /**
         * @brief Type for dense Eigen matrices.
         */
        typedef Eigen::DenseBase<Real> DenseBase;

        /**
         * @brief Type for Eigen matrix block.
         */
        typedef Eigen::Block<MatrixX> Block;

        // Eigen shortcuts for vectors
        typedef Eigen::Matrix<Real, 3, 1> Vector3;
        typedef Eigen::Matrix<Real, 4, 1> Vector4;
        typedef Eigen::Matrix<Real, 5, 1> Vector5;
        typedef Eigen::Matrix<Real, 6, 1> Vector6;
        typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorX;
        typedef Eigen::Map<VectorX> MapVectorX;
        typedef std::vector<MapVectorX> VecMapVX;

        /**
         * @brief Type for integer pairs (e.g., dimensions).
         */
        typedef std::pair<unsigned int, unsigned int> IntPair;
        /**
         * @brief Type for a vector of integer pairs (e.g., dimensions for multiple surfaces).
         */
        typedef std::vector<IntPair> VecDimensions;

        /**
         * @brief Structure representing various options for the VM (Vortex Method) solver.
         */
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
            bool only_lifting;
            bool only_nonlifting;
            bool phantom_wing_test;
        	unsigned int NumCores;
        	unsigned int NumSurfaces;
        	unsigned int NumSurfacesNonlifting;
            double dt;
            unsigned int n_rollup;
            double rollup_tolerance;
            unsigned int rollup_aic_refresh;
            bool iterative_solver;
            double iterative_tol;
            bool iterative_precond;
            double vortex_radius;
            double vortex_radius_wake_ind;
            bool consider_u_ind_by_sources_for_lifting_forces;
            uint ignore_first_x_nodes_in_force_calculation;
        };

        /**
         * @brief Structure representing various options for the UVM (Unsteady Vortex Method) solver.
         */
        struct UVMopts
        {
            double dt;
            uint NumCores;
            uint NumSurfaces;
            uint NumSurfacesNonlifting;
            bool only_lifting;
            bool only_nonlifting;
            bool phantom_wing_test;
            uint convection_scheme;
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
            uint num_spanwise_panels_wo_induced_velocity;
            bool consider_u_ind_by_sources_for_lifting_forces;
            uint ignore_first_x_nodes_in_force_calculation;
        };
         
         /**
         * @brief This function creates a VMopts object from a UVMopts object.
         */
        VMopts UVMopts2VMopts(const UVMopts& uvm)
        {
            VMopts vm;
            vm.dt = uvm.dt;
            vm.NumCores = uvm.NumCores;
            vm.NumSurfaces = uvm.NumSurfaces;
            vm.NumSurfacesNonlifting = uvm.NumSurfacesNonlifting;
            vm.ImageMethod = uvm.ImageMethod;
            vm.iterative_solver = uvm.iterative_solver;
            vm.iterative_tol = uvm.iterative_tol;
            vm.iterative_precond = uvm.iterative_precond;
            vm.horseshoe = false;
            vm.vortex_radius = uvm.vortex_radius;
            vm.vortex_radius_wake_ind = uvm.vortex_radius_wake_ind;
            vm.only_lifting = uvm.only_lifting;
            vm.only_nonlifting = uvm.only_nonlifting;
            vm.Steady = uvm.quasi_steady;
            vm.phantom_wing_test = uvm.phantom_wing_test;
            vm.consider_u_ind_by_sources_for_lifting_forces = uvm.consider_u_ind_by_sources_for_lifting_forces;
            vm.ignore_first_x_nodes_in_force_calculation = uvm.ignore_first_x_nodes_in_force_calculation;
            return vm;
        };

         /**
         * @brief Structure representing flight conditions.
         * 
         * Structure includes flow speed and direction, as well as reference chord and density.
         */
        struct FlightConditions
        {
            double uinf = 1.0;
            double uinf_direction[3];
            double rho = 1.225;
            double c_ref = 1.0;
        };

        /**
         * @brief Correct the dimensions by applying a correction value.
         *
         * This function corrects dimensions by applying a correction value. It ensures that the dimensions
         * are not negative after the correction.
         *
         * @param correction The correction value to be applied.
         * @param dimension The original dimension.
         * @return The corrected dimension.
         */
        inline int correct_dimensions
        (
            const int correction,
            int dimension
        )
        {
            if (correction < 0 && dimension < abs(correction))
            {
                return 0;
            }
            else
            {
                return dimension + correction;
            }
        }

        /**
         * @brief Generate dimensions from a matrix.
         *
         * This function generates dimensions from a matrix and stores them in a vector of integer pairs.
         * It applies a correction to the dimensions if needed.
         *
         * @param mat The input matrix.
         * @param dimensions The vector of dimensions to be filled.
         * @param correction The correction value to be applied (default is 0).
         */
        template <typename t_mat>
        inline void generate_dimensions
        (
            const t_mat& mat,
            UVLM::Types::VecDimensions& dimensions,
            const int& correction = 0
        )
        {
            int M, N;
            dimensions.resize(mat.size());
            
            for (unsigned int i_surf=0; i_surf<dimensions.size(); ++i_surf)
            {
                M = correct_dimensions(correction, mat[i_surf][0].rows());
                N = correct_dimensions(correction, mat[i_surf][0].cols());
                dimensions[i_surf] = UVLM::Types::IntPair(M, N);
            }
        }
        
        /**
         * @brief Allocate a vector of matrices with specified dimensions.
         *
         * This function allocates a vector of matrices with specified dimensions and initializes them with
         * an initial value.
         *
         * @param mat The vector of matrices to be allocated.
         * @param dimensions The dimensions for each matrix.
         * @param correction The correction value to be applied to dimensions (default is 0).
         * @param initial_value The initial value to fill the matrices (default is 0.0).
         */
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

        inline void allocate_VecMat_from_VecVecMat
        (
            UVLM::Types::VecMatrixX& mat,
            const UVLM::Types::VecVecMatrixX& dimensions_in,
            const int& correction = 0,
            const double& initial_value = 0.0
        )
        {
            const unsigned int n_mats = dimensions_in.size();
            mat.resize(n_mats);
            for (unsigned int i=0; i<n_mats; ++i)
            {
                if ((dimensions_in[i][0].rows() == 0 )|| (dimensions_in[i][0].cols() == 0))
                {
                    mat[i].resize(0,0);
                }

                else
                {
                    mat[i].setConstant
                    (
                        dimensions_in[i][0].rows() + correction,
                        dimensions_in[i][0].cols() + correction,
                        initial_value
                    );
                }
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
        /**
         * @brief Initialize a vector of matrices with a specific value.
         *
         * This function initializes a vector of matrices with a specific value.
         *
         * @param mat The vector of matrices to be initialized.
         * @param value The value to initialize the matrices with (default is 0.0).
         */
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
        /**
         * @brief Allocate a vector of vector of matrices with specified dimensions.
         *
         * This function allocates a vector of vector of matrices with specified dimensions and initializes
         * them with an initial value.
         *
         * @param mat The vector of vector of matrices to be allocated.
         * @param n_surf The number of surfaces.
         * @param n_dim The number of dimensions for each matrix.
         * @param M The number of rows for each matrix.
         * @param N The number of columns for each matrix.
         */
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
            int M, N;
            mat.resize(n_surf);
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                M = correct_dimensions(correction, dimensions[i_surf].first);
                N = correct_dimensions(correction, dimensions[i_surf].second);
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
            int M, N;
            mat.resize(n_surf);
            for (unsigned int i_surf=0; i_surf<n_surf; ++i_surf)
            {
                
                M = correct_dimensions(correction, in_dimensions[i_surf][0].rows());
                N = correct_dimensions(correction, in_dimensions[i_surf][0].cols());  
                mat[i_surf].resize(in_dimensions[i_surf].size());
                for (unsigned int i_dim=0; i_dim<in_dimensions[i_surf].size(); ++i_dim)
                {
                    mat[i_surf][i_dim].resize(M, N);
                    mat[i_surf][i_dim].setZero(M, N);
                }
            }
        }
        /**
         * @brief Copy data from a vector of vector of matrices to another.
         *
         * This function copies data from one vector of vector of matrices to another.
         *
         * @tparam t_in The input vector of vector of matrices.
         * @tparam t_out The output vector of vector of matrices.
         * @param in The input vector of vector of matrices.
         * @param out The output vector of vector of matrices.
         */
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
        /**
         * @brief Copy data from a matrix to a block within another matrix.
         *
         * This function copies data from a matrix to a specified block within another matrix.
         *
         * @tparam mat_in The input matrix.
         * @param in The input matrix to copy from.
         * @param out The output matrix (block) to copy to.
         * @param i_start The row index to start copying.
         * @param j_start The column index to start copying.
         */
        template<typename mat_in>
        inline void copy_Mat_to_block
        (
            mat_in& in,
            UVLM::Types::MatrixX& out,
            uint i_start,
            uint j_start
        )
        {
            uint M = in.rows();
            uint N = in.cols();
            for (uint i_m=0; i_m<M; ++i_m)
            {
                for (uint i_n=0; i_n<N; ++i_n)
                {
                    out(i_m+i_start, i_n+j_start) = in(i_m, i_n);
                }
            }
        }
        /**
         * @brief Join two vectors into a single vector.
         *
         * This function concatenates two vectors into a single vector.
         *
         * @param vec1 The first vector.
         * @param vec2 The second vector.
         * @return The concatenated vector.
         */
        UVLM::Types::VectorX join_vectors
        (
            const UVLM::Types::VectorX& vec1,
            const UVLM::Types::VectorX& vec2
        )
        {
            UVLM::Types::VectorX vec_joined(vec1.size() + vec2.size());
            vec_joined << vec1, vec2;
            return vec_joined;
        }

        /**
         * @brief Calculate the norm of a vector of vector of matrices.
         *
         * This function calculates the norm of a vector of vector of matrices by summing the norms of individual matrices.
         *
         * @tparam t_mat The input vector of vector of matrices.
         * @param mat The input vector of vector of matrices.
         * @return The norm of the vector of vector of matrices.
         */        template <typename t_mat>
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
        /**
         * @brief Calculate the norm of a vector.
         *
         * This function calculates the norm of a 3D vector.
         *
         * @param value_1 The first component of the vector.
         * @param value_2 The second component of the vector.
         * @param value_3 The third component of the vector.
         * @return The norm of the vector.
         */

        inline double norm_Vec
        (
            const double value_1,
            const double value_2,
            const double value_3
        )
        {
            
            UVLM::Types::Vector3 vector;
            vector << value_1, value_2, value_3;
            return vector.norm();
        }

        /**
         * @brief Calculate the squared norm of a vector.
         *
         * This function calculates the squared norm of a 3D vector.
         *
         * @param value_1 The first component of the vector.
         * @param value_2 The second component of the vector.
         * @param value_3 The third component of the vector.
         * @return The squared norm of the vector.
         */
        inline double norm_Vec_squared
        (
            const double value_1,
            const double value_2,
            const double value_3
        )
        {
            
            return value_1 * value_1 + 
                   value_2 * value_2 + 
                   value_3 * value_3;
        }
        /**
         * @brief Calculate the maximum value in a vector of vector of matrices.
         *
         * This function calculates the maximum value in a vector of vector of matrices by finding the maximum
         * absolute value in all matrices.
         *
         * @tparam t_mat The input vector of vector of matrices.
         * @param mat The input vector of vector of matrices.
         * @return The maximum value in the vector of vector of matrices.
         */
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

        template<typename t_vecmatrix>
        inline void pass_3D_Vector_to_VecMatrix
        (
            t_vecmatrix& mat,
            UVLM::Types::Vector3& vec,
            const uint row,
            const uint col

        )
        {
            for (uint i_dim=0; i_dim<3; i_dim++)
            {
                mat[i_dim](row, col) = vec(i_dim);
            }
        }

        void remove_row_from_VectorX
		(
		UVLM::Types::VectorX & vector_in,
		const int rowToRemove
		)
		{
			int counter_idx = 0;
            UVLM::Types::VectorX intitial_vector_in = vector_in;
            uint vector_initial_size = vector_in.rows();
            vector_in.resize(vector_initial_size-1,1);
			for (uint i_row = 0;i_row<vector_initial_size;++i_row)
			{
				if (i_row!=rowToRemove)
				{
					vector_in[counter_idx] = intitial_vector_in[i_row];
					counter_idx++;
				}
			}
		}

        template<typename vec_in>
        UVLM::Types::VectorX reorder_vector_by_pushback
		(
		vec_in& vector_in,
		const int numb_of_push_backs
		)
		{
            uint vector_size = vector_in.rows();            
            UVLM::Types::VectorX vec_out(vector_size);
            uint counter = 0;
			for (uint i_row = numb_of_push_backs;i_row<vector_size;++i_row)
			{
                vec_out[counter] = vector_in[i_row];
                counter++;
			}
            for (uint i_row = 0 ;i_row<numb_of_push_backs;++i_row)
			{
                vec_out[counter] = vector_in[i_row];
                counter++;
			}
            return vec_out;
		}
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
