/*
 * wrapper.cpp
 *
 *  Created on: 29 Jan 2013
 *      Author: rjs10
 */

#include "vorticity.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <VLM.hpp>
#include <stdio.h>
#include <datatypesx.hpp>
#include <aicMats.hpp>
#include <indices.hpp>

// forward declare some useful functions

Eigen::Map<Eigen::Vector3d> C2Eig_triad_map(double* x_triad) {
	return Eigen::Map<Eigen::Vector3d>(x_triad);
}

Eigen::Vector3d C2Eig_triad(const double* x_triad) {
	return Eigen::Vector3d(x_triad[0], x_triad[1], x_triad[2]);
}

void Eig2C_triad(const Eigen::Vector3d& Eig, double* x_triad) {
	x_triad[0] = Eig[0];
	x_triad[1] = Eig[1];
	x_triad[2] = Eig[2];
}

//wrapper interface must by in C-syntax
extern "C" {

void cpp_wrap_vorticity_biotsegment(double* xp_triad, double* x1_triad, \
									 double* x2_triad, double& Gamma, \
									 double* Uind_triad) {
	/** @brief wrapper for BiotSegment.
	 *  @details c-style arrays are converted into Eigen::Vector3d.
	 */

	//call BiotSegment
	Eigen::Vector3d Uind = BiotSegment(C2Eig_triad(xp_triad),\
									   C2Eig_triad(x1_triad),\
									   C2Eig_triad(x2_triad),\
									   Gamma);

	//overwrite Uind with result
	Eig2C_triad(Uind,Uind_triad);

	return;
}

void cpp_wrap_vorticity_biotsegment_map(double* xp_triad, double* x1_triad, \
									 double* x2_triad, double& Gamma, \
									 double* Uind_triad) {
	/** @brief wrapper for BiotSegment.
	 *  @details c-style arrays are converted into Eigen::Vector3d.
	 */

	//call BiotSegment
	BiotSegment_map(C2Eig_triad_map(xp_triad),\
					C2Eig_triad_map(x1_triad),\
					C2Eig_triad_map(x2_triad),\
					Gamma,\
					C2Eig_triad_map(Uind_triad) );

	//overwrite Uind with result
	//Eig2C_triad(Uind,Uind_triad);

	return;
}

void c_wrap_vorticity_biotsegment(double* xp_triad, double* x1_triad, \
									 double* x2_triad, double& Gamma, \
									 double* Uind_triad) {
	/** @brief wrapper for C_BiotSegment.
	  * @details c-style pointers only.
	  */
	C_BiotSegment(xp_triad, x1_triad, x2_triad,\
				  Gamma, Uind_triad);
}

void cpp_wrap_test_biotsegment(const int& NumTests) {
	Eigen::Vector3d xp(0.0,0.0,-1.0);
	Eigen::Vector3d x1(-0.5,0.0,0.0);
	Eigen::Vector3d x2(0.5,0.0,0.0);
	double gam  = 1.0;
	Eigen::Vector3d Result(0.0,0.0,0.0);

	//call function
	for (int i = 0; i < NumTests; i++) {
		Result = BiotSegment(xp,x1,x2,gam);
	}
	return;
}

void c_wrap_test_biotsegment(const int& NumTests) {
	double xp[] = {0.0,0.0,-1.0};
	double x1[] = {-0.5,0.0,0.0};
	double x2[] = {0.5,0.0,0.0};
	double gam  = 1.0;
	double Result[] = {0.0,0.0,0.0};

	//call function
	for (int i = 0; i < NumTests; i++) {
		C_BiotSegment(xp,x1,x2,gam,Result);
	}
	return;
}

void cpp_wrap_solver_vlm(const double* Zeta_Vec, const double* ZetaDot_Vec, \
		const double* Uext_Vec, \
		double* ZetaStar_Vec, \
		unsigned int& VMOPTS_M, \
		unsigned int& VMOPTS_N, \
		bool& VMOPTS_ImageMethod,\
		unsigned int& VMOPTS_Mstar,\
		bool& VMOPTS_Steady,\
		bool& VMOPTS_KJMeth,\
		bool& VMOPTS_NewAIC,\
		double& VMOPTS_DelTime,\
		bool& VMOPTS_Rollup,\
		unsigned int& VMOPTS_NumCores,\
		double* Forces_Vec, \
		double* Gamma_Vec, \
		double* GammaStar_Vec,\
		double* AIC_Vec,\
	    double* BIC_Vec) {
	/** @brief wrapper for cpp_solver_vlm
	 *  @details wrapped function takes c arguments only
	 */

	// Convert VMOPTS_* into class
	VMopts VMOPTS;
	VMOPTS.M = VMOPTS_M;
	VMOPTS.N = VMOPTS_N;
	VMOPTS.ImageMethod = VMOPTS_ImageMethod;
	VMOPTS.Mstar = VMOPTS_Mstar;
	VMOPTS.Steady = VMOPTS_Steady;
	VMOPTS.KJMeth = VMOPTS_KJMeth;
	VMOPTS.NewAIC = VMOPTS_NewAIC;
	VMOPTS.DelTime = VMOPTS_DelTime;
	VMOPTS.Rollup = VMOPTS_Rollup;
	VMOPTS.NumCores = VMOPTS_NumCores;

	//call solver
	cpp_solver_vlm(Zeta_Vec, ZetaDot_Vec, Uext_Vec, ZetaStar_Vec, \
			VMOPTS, Forces_Vec, Gamma_Vec, GammaStar_Vec, \
			AIC_Vec,\
			BIC_Vec);

}

void cpp_wrap_KJMethodForces(const double* Zeta_Vec,
					  const double* Gamma_Vec,
					  const double* ZetaStar_Vec,
					  const double* GammaStar_Vec,
					  const double* ZetaDot_Vec,
					  const double* Uext_Vec,
					  unsigned int& VMOPTS_M,
					  unsigned int& VMOPTS_N,
					  bool& VMOPTS_ImageMethod,
					  unsigned int& VMOPTS_Mstar,
					  bool& VMOPTS_Steady,
					  bool& VMOPTS_KJMeth,
					  bool& VMOPTS_NewAIC,
					  double& VMOPTS_DelTime,
					  bool& VMOPTS_Rollup,
					  unsigned int& VMOPTS_NumCores,
					  const double* Gamma_tm1_Vec,
					  double* Forces_Vec) {
	/** @brief wrapper for cpp_solver_vlm
	 */

	// Convert VMOPTS_* into class
	VMopts VMOPTS;
	VMOPTS.M = VMOPTS_M;
	VMOPTS.N = VMOPTS_N;
	VMOPTS.ImageMethod = VMOPTS_ImageMethod;
	VMOPTS.Mstar = VMOPTS_Mstar;
	VMOPTS.Steady = VMOPTS_Steady;
	VMOPTS.KJMeth = VMOPTS_KJMeth;
	VMOPTS.NewAIC = VMOPTS_NewAIC;
	VMOPTS.DelTime = VMOPTS_DelTime;
	VMOPTS.Rollup = VMOPTS_Rollup;
	VMOPTS.NumCores = VMOPTS_NumCores;

	KJMethodForces(Zeta_Vec, Gamma_Vec, ZetaStar_Vec, GammaStar_Vec,
				   ZetaDot_Vec, Uext_Vec, VMOPTS, Gamma_tm1_Vec, Forces_Vec);

}

void cpp_wrap_KJMethodForces_vC(const double* Zeta_Vec,
					  const double* Gamma_Vec,
					  const double* ZetaStar_Vec,
					  const double* GammaStar_Vec,
					  const double* ZetaDot_Vec,
					  const double* Uext_Vec,
					  unsigned int& VMOPTS_M,
					  unsigned int& VMOPTS_N,
					  bool& VMOPTS_ImageMethod,
					  unsigned int& VMOPTS_Mstar,
					  bool& VMOPTS_Steady,
					  bool& VMOPTS_KJMeth,
					  bool& VMOPTS_NewAIC,
					  double& VMOPTS_DelTime,
					  bool& VMOPTS_Rollup,
					  unsigned int& VMOPTS_NumCores,
					  const double* Gamma_tm1_Vec,
					  double* Forces_Vec) {
	/** @brief wrapper for cpp_solver_vlm
	 */

	// Convert VMOPTS_* into class
	VMopts VMOPTS;
	VMOPTS.M = VMOPTS_M;
	VMOPTS.N = VMOPTS_N;
	VMOPTS.ImageMethod = VMOPTS_ImageMethod;
	VMOPTS.Mstar = VMOPTS_Mstar;
	VMOPTS.Steady = VMOPTS_Steady;
	VMOPTS.KJMeth = VMOPTS_KJMeth;
	VMOPTS.NewAIC = VMOPTS_NewAIC;
	VMOPTS.DelTime = VMOPTS_DelTime;
	VMOPTS.Rollup = VMOPTS_Rollup;
	VMOPTS.NumCores = VMOPTS_NumCores;

	KJMethodForces_vC(Zeta_Vec, Gamma_Vec, ZetaStar_Vec, GammaStar_Vec,
				   ZetaDot_Vec, Uext_Vec, VMOPTS, Gamma_tm1_Vec, Forces_Vec);

}

void cpp_wrap_KJMethodForces_vC_mod(const double* Zeta_Vec,
					  const double* Gamma_Vec,
					  const double* ZetaStar_Vec,
					  const double* GammaStar_Vec,
					  const double* ZetaDot_Vec,
					  const double* Uext_Vec,
					  unsigned int& VMOPTS_M,
					  unsigned int& VMOPTS_N,
					  bool& VMOPTS_ImageMethod,
					  unsigned int& VMOPTS_Mstar,
					  bool& VMOPTS_Steady,
					  bool& VMOPTS_KJMeth,
					  bool& VMOPTS_NewAIC,
					  double& VMOPTS_DelTime,
					  bool& VMOPTS_Rollup,
					  unsigned int& VMOPTS_NumCores,
					  const double* Gamma_tm1_Vec,
					  double* Forces_Vec) {
	/** @brief wrapper for cpp_solver_vlm
	 */

	// Convert VMOPTS_* into class
	VMopts VMOPTS;
	VMOPTS.M = VMOPTS_M;
	VMOPTS.N = VMOPTS_N;
	VMOPTS.ImageMethod = VMOPTS_ImageMethod;
	VMOPTS.Mstar = VMOPTS_Mstar;
	VMOPTS.Steady = VMOPTS_Steady;
	VMOPTS.KJMeth = VMOPTS_KJMeth;
	VMOPTS.NewAIC = VMOPTS_NewAIC;
	VMOPTS.DelTime = VMOPTS_DelTime;
	VMOPTS.Rollup = VMOPTS_Rollup;
	VMOPTS.NumCores = VMOPTS_NumCores;

	KJMethodForces_vC_mod(Zeta_Vec, Gamma_Vec, ZetaStar_Vec, GammaStar_Vec,
				   ZetaDot_Vec, Uext_Vec, VMOPTS, Gamma_tm1_Vec, Forces_Vec);

}

void cpp_wrap_AIC(const double* zetaSrc_,
				  const unsigned int& mSrc,
				  const unsigned int& nSrc,
				  const double* zetaTgt_,
				  const unsigned int& mTgt,
				  const unsigned int& nTgt,
				  const bool& imageMeth,
				  double* Aic_) {
	//call AIC generator
	AIC(zetaSrc_,mSrc,nSrc,zetaTgt_,mTgt,nTgt,imageMeth,Aic_);
}

void cpp_wrap_dAgamma0_dZeta(const double* zetaSrc_,
								const unsigned int& mSrc,
								const unsigned int& nSrc,
								const double* gamma0_,
								const double* zetaTgt_,
								const unsigned int& mTgt,
								const unsigned int& nTgt,
//<<<<<<< HEAD
//								double* dAgam0_) {
//	// call function
//	dAgamma0_dZeta(zetaSrc_,mSrc,nSrc,gamma0_,zetaTgt_,mTgt,nTgt,dAgam0_);
//=======
								const bool& imageMeth,
								double* dAgam0_) {
	// call function
	dAgamma0_dZeta(zetaSrc_,mSrc,nSrc,gamma0_,zetaTgt_,mTgt,nTgt,imageMeth,dAgam0_);
//>>>>>>> rob/master
}

void cpp_wrap_dWzetaPri0_dZeta(const double* zeta_,
								  const unsigned int& m,
								  const unsigned int& n,
								  const double* zetaPri_,
								  double* dX_) {
	//call function
	dWzetaPri0_dZeta(zeta_, m, n, zetaPri_, dX_);
}

void cpp_wrap_genW(const double* zeta_,
					  const unsigned int& M,
					  const unsigned int& N,
					  double* W_){
	genW(zeta_,M,N,W_);
}

void cpp_wrap_genH(const unsigned int& m,
		   	   	     const unsigned int& n,
		   	   	     double* H_){
	genH(m,n,H_);
}

void cpp_wrap_genXi(const unsigned int& m,
		    		 const unsigned int& n,
		    		 const double& eta1,
		    		 const double& eta2,
		    		 double* Xi_) {
	genXi(m,n,eta1,eta2,Xi_);
}

void cpp_wrap_AIC3(const double* zetaSrc_,
		  const unsigned int& mSrc,
		  const unsigned int& nSrc,
		  const double* zetaTgt_,
		  const unsigned int& mTgt,
		  const unsigned int& nTgt,
		  double* dX_){
	AIC3(zetaSrc_,mSrc,nSrc,zetaTgt_,mTgt,nTgt,dX_);
}

void cpp_wrap_AIC3noTE(const double* zetaSrc_,
		  const unsigned int& mSrc,
		  const unsigned int& nSrc,
		  const double* zetaTgt_,
		  const unsigned int& mTgt,
		  const unsigned int& nTgt,
		  const bool& wakeSrc,
		  double* dX_){
	AIC3noTE(zetaSrc_,mSrc,nSrc,zetaTgt_,mTgt,nTgt,wakeSrc,dX_);
}

void cpp_wrap_dA3gamma0_dZeta(const double* zetaSrc_,
								const unsigned int& mSrc,
								const unsigned int& nSrc,
								const double* gamma0_,
								const double* zetaTgt_,
								const unsigned int& mTgt,
								const unsigned int& nTgt,
								double* dA3gam0_) {
	// call function
	dA3gamma0_dZeta(zetaSrc_,mSrc,nSrc,gamma0_,zetaTgt_,mTgt,nTgt,dA3gam0_);
}

void cpp_wrap_Y1(const double* vM_,
				   const double* zeta_,
				   const unsigned int& m,
				   const unsigned int& n,
				   double* Y1_) {
	Y1(vM_,zeta_,m,n,Y1_);
}

void cpp_wrap_Y2(const double* gamma_,
				   const double* vM_,
				   const unsigned int& m,
				   const unsigned int& n,
				   double* Y2_) {
	Y2(gamma_,vM_,m,n,Y2_);
}

void cpp_wrap_Y3(const double* gamma_,
				   const double* zeta_,
				   const unsigned int& m,
				   const unsigned int& n,
				   double* Y3_) {
	Y3(gamma_,zeta_,m,n,Y3_);
}

void cpp_wrap_Y4(const double* zeta_,
				  const unsigned int& m,
				  const unsigned int& n,
				  double* Y4_) {
	Y4(zeta_,m,n,Y4_);
}

void cpp_wrap_Y5(const double* gammaPri_,
				   const double* zeta_,
				   const unsigned int& m,
				   const unsigned int& n,
				   double* Y5_) {
	Y5(gammaPri_,zeta_,m,n,Y5_);
}

void cpp_wrap_AIC3s(const double* zetaSrc_,
		  const unsigned int& mSrc,
		  const unsigned int& nSrc,
		  const double* zetaTgt_,
		  const unsigned int& mTgt,
		  const unsigned int& nTgt,
//<<<<<<< HEAD
//		  double* dX_){
//	AIC3s(zetaSrc_,mSrc,nSrc,zetaTgt_,mTgt,nTgt,dX_);
//=======
		  const bool& imageMeth,
		  double* dX_){
	AIC3s(zetaSrc_,mSrc,nSrc,zetaTgt_,mTgt,nTgt,imageMeth,dX_);
//>>>>>>> rob/master
}

void cpp_wrap_AIC3s_noTE(const double* zetaSrc_,
		  const unsigned int& mSrc,
		  const unsigned int& nSrc,
		  const double* zetaTgt_,
		  const unsigned int& mTgt,
		  const unsigned int& nTgt,
		  const bool& wakeSrc,
		  double* dX_){
	AIC3s_noTE(zetaSrc_,mSrc,nSrc,zetaTgt_,mTgt,nTgt,wakeSrc,dX_);
}

unsigned int cpp_wrap_q_k(const unsigned int k_,
				   	         const unsigned int N,
				   	         const unsigned int no) {
	return(q_k(k_,N,no));
}

void cpp_wrap_dAs3gam0_dZeta_numerical(const double* zetaSrc_,
									  const unsigned int& mSrc,
									  const unsigned int& nSrc,
									  const double* gamma0_,
									  const double* zetaTgt_,
									  const unsigned int& mTgt,
									  const unsigned int& nTgt,
//<<<<<<< HEAD
//									  double* dX_) {
//	dAs3gam0_dZeta_numerical(zetaSrc_,mSrc,nSrc,gamma0_,zetaTgt_,mTgt,nTgt,dX_);
//=======
									  const bool& imageMeth,
									  double* dX_) {
	dAs3gam0_dZeta_numerical(zetaSrc_,mSrc,nSrc,gamma0_,zetaTgt_,mTgt,nTgt,imageMeth,dX_);
//>>>>>>> rob/master
}

} // END extern C