/*
 * aicMats.hpp
 *
 *  Created on: 13 May 2015
 *      Author: rjs10
 */

#ifndef AICMATS_HPP_
#define AICMATS_HPP_

#include<Eigen/Dense>
#include<datatypesx.hpp>
using namespace Eigen;

void genW(const double* zeta_,
		   const unsigned int M,
		   const unsigned int N,
		   double* W_);

void getNormals(const VectorXd& zeta,
		          const unsigned int M,
		          const unsigned int N,
		          const VectorXd& normals_);

void genAstar(const VectorXd& zetaSrc,
				const unsigned int mSrc,
				const unsigned int nSrc,
		        const VectorXd& zetaTgt,
		        const unsigned int mTgt,
		        const unsigned int nTgt,
		        const EigDynMatrixRM& Astar_);

double fGeom(const double* r0,
          const double* r1,
	 	  const double* r2,
	 	  const double* n);

void fGeom3(const double* r0,
          const double* r1,
	 	  const double* r2,
	 	  double* out);

void df_dgeom(const double* r0,
		        const double* r1,
			 	const double* r2,
			 	const double* n,
			 	const Vector3d& f_r0_,
			 	const Vector3d& f_r1_,
			 	const Vector3d& f_r2_,
			 	const Vector3d& f_n_);

Matrix3d dxHat_dx(const Vector3d& x);

Matrix3d duHat2_dr1(const Vector3d& r1,
				     const Vector3d& r2);

Matrix3d duHat2_dr2(const Vector3d& r1,
				     const Vector3d& r2);

Matrix3d skew(const Vector3d& x);

Matrix3d dn_dd(const Vector3d& d, const Vector3d& e);

Matrix3d dn_de(const Vector3d& d, const Vector3d& e);

void dAgamma0_dZeta(const double* zetaSrc_,
					 const unsigned int mSrc,
					 const unsigned int nSrc,
					 const double* gamma0_,
					 const double* zetaTgt_,
					 const unsigned int mTgt,
					 const unsigned int nTgt,
//<<<<<<< HEAD
//					 double* dX_);
//=======
					 const bool imageMeth,
					 double* dX_);

void dAgamma0_dZeta_num(const double* zetaSrc_,
					 	   const unsigned int mSrc,
					 	   const unsigned int nSrc,
					 	   const double* gamma0_,
					 	   const double* zetaTgt_,
					 	   const unsigned int mTgt,
					 	   const unsigned int nTgt,
					 	   const bool imageMeth,
					 	   double* dX_);

//>>>>>>> rob/master
void AIC(const double* zetaSrc_,
		  const unsigned int mSrc,
		  const unsigned int nSrc,
		  const double* zetaTgt_,
		  const unsigned int mTgt,
		  const unsigned int nTgt,
		  const bool imageMeth,
		  double* dX_);

void dWzetaPri0_dZeta(const double* zeta_,
					  const unsigned int m,
					  const unsigned int n,
					  const double* zetaPri_,
					  double* dX_);

void genH(const unsigned int m,
		   const unsigned int n,
		   double* H_);

void genXi(const unsigned int m,
		    const unsigned int n,
		    const double eta1,
		    const double eta2,
		    double* Xi_);

void AIC3(const double* zetaSrc_,
		  const unsigned int mSrc,
		  const unsigned int nSrc,
		  const double* zetaTgt_,
		  const unsigned int mTgt,
		  const unsigned int nTgt,
		  double* dX_);

void AIC3noTE(const double* zetaSrc_,
		  const unsigned int mSrc,
		  const unsigned int nSrc,
		  const double* zetaTgt_,
		  const unsigned int mTgt,
		  const unsigned int nTgt,
		  const bool wakeSrc,
		  double* dX_);

void AIC3s(const double* zetaSrc_,
		  	const unsigned int mSrc,
		  	const unsigned int nSrc,
		  	const double* zetaTgt_,
		  	const unsigned int mTgt,
		  	const unsigned int nTgt,
//<<<<<<< HEAD
//=======
		  	const bool imageMeth,
//>>>>>>> rob/master
		  	double* dX_);

void AIC3s_noTE(const double* zetaSrc_,
		  	const unsigned int mSrc,
		  	const unsigned int nSrc,
		  	const double* zetaTgt_,
		  	const unsigned int mTgt,
		  	const unsigned int nTgt,
		  	const bool wakeSrc,
		  	double* dX_);

void dA3gamma0_dZeta(const double* zetaSrc_,
					 const unsigned int mSrc,
					 const unsigned int nSrc,
					 const double* gamma0_,
					 const double* zetaTgt_,
					 const unsigned int mTgt,
					 const unsigned int nTgt,
					 double* dX_);

void Y1(const double* vM_,
		const double* zeta_,
		const unsigned int m,
		const unsigned int n,
		double* Y1_);

void Y2(const double* gamma_,
		const double* vM_,
		const unsigned int m,
		const unsigned int n,
		double* Y2_);

void Y3(const double* gamma_,
		const double* zeta_,
		const unsigned int m,
		const unsigned int n,
		double* Y3_);

void Y4(const double* zeta_,
		const unsigned int m,
		const unsigned int n,
		double* Y4_);

void Y5(const double* gammaPri_,
		const double* zeta_,
		const unsigned int m,
		const unsigned int n,
		double* Y5_);

void dAs3gam0_dZeta_numerical(const double* zetaSrc_,
								const unsigned int mSrc,
								const unsigned int nSrc,
								const double* gamma0_,
								const double* zetaTgt_,
								const unsigned int mTgt,
								const unsigned int nTgt,
//<<<<<<< HEAD
//=======
								const bool imageMeth,
//>>>>>>> rob/master
								double* dX_);

#endif /* AICMATS_HPP_ */
