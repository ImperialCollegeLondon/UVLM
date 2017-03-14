/*
 * main.cpp
 *
 *  Created on: 28 Jan 2013
 *      Author: rjs10
 */

#include <aicMats.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <triads.hpp>
#include <indices.hpp>
using namespace std;
using namespace Eigen;

double fRand(double fMin, double fMax) {
	// generate random numbers
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void wFunc(const MatrixXd& A,
		    const VectorXd& gamma,
		    const VectorXd& w_) {
	VectorXd& w = const_cast<VectorXd&> (w_);
	w = A*gamma;
}

int main() {
	// Test UVLMLib functions during development.
	srand(time(NULL));

	// grid vars
	const int M=1;
//<<<<<<< HEAD
//	const int N=1;
//=======
	const int N=2;
//>>>>>>> rob/master
	const int K=M*N;
	const int K_zeta = (M+1)*(N+1);

	Eigen::VectorXd zeta(3*K_zeta); zeta.setZero();
	Eigen::MatrixXd W(K,3*K_zeta); W.setZero();
	Eigen::MatrixXd Astar(3*K_zeta,K); Astar.setZero();
	Eigen::VectorXd gam(K); gam.setOnes(); // init gamma as 1
//<<<<<<< HEAD
//	Eigen::MatrixXd dAgam0_dzeta(K,3*K_zeta); dAgam0_dzeta.setZero();
//
////	zeta << 0.0, 0.0, 0.0,
////			0.0, 1.0, 0.0,
////			0.0, 2.0, 0.0,
////			1.0, 0.0, 0.0,
////			1.0, 1.0, 0.0,
////			1.0, 2.0, 0.0;
//
//	zeta << 0.0,    -1000.0, 0.0,
//			0.0,     1000.0, 0.0,
//			0.1,    -1000.0, 0.0,
//			0.1,     1000.0, 0.0;
//=======
	EigDynMatrixRM dAgam0_dzeta(K,3*K_zeta); dAgam0_dzeta.setZero();

//	zeta << 0.0, 0.0, 0.0,
//			0.0, 1.0, 0.0,
//			1.0, 0.0, 0.0,
//			1.0, 1.0, 0.0,
//			2.0, 0.0, 0.0,
//			2.0, 1.0, 0.0;

	zeta << 0.0, 1.0, 0.0,
			0.0, 2.0, 0.0,
			0.0, 3.0, 0.0,
			1.0, 1.0, 0.0,
			1.0, 2.0, 0.0,
			1.0, 3.0, 0.0;

//	zeta << 0.0, 0.0, 0.0,
//			0.0, 1.0, 0.0,
//			1.0, 0.0, 0.0,
//			1.0, 1.0, 0.0;

//	zeta << 0.0,    -1000.0, 0.0,
//			0.0,     1000.0, 0.0,
//			0.1,    -1000.0, 0.0,
//			0.1,     1000.0, 0.0;
//>>>>>>> rob/master

	const double* zetaPtr = zeta.data();

	VectorXd dZeta = VectorXd::Random(zeta.size())/100;
// <<<<<<< HEAD
// 	VectorXd zetaPdel = zeta + dZeta;
// 	const double* zetaPdelPtr = zetaPdel.data();

// 	genW(zetaPtr,M,N,W.data());

// 	genAstar(zeta,M,N,zeta,M,N,Astar);

// 	// test dAgam_dzeta subroutines
// 	// dxHat_dx
// 	Matrix3d dX;
// 	Vector3d x(1.0,1.0,1.0);
// 	Vector3d dx = Vector3d::Random()/(10.0*x.norm());
// 	Vector3d f = (x + dx).normalized();
// 	dX = dxHat_dx(x);
// 	Vector3d fApprox = x.normalized() + dX * dx;
// 	Vector3d diff = f-fApprox;
// 	cout << "dxHat_dx: ---------------" << endl << "dx:" << endl << dx
// 	     << endl << "norm (l-2) rel error:" << endl << diff.norm()/dx.norm()
// 	     << endl
// 	     << endl;

// 	// Init matrices for df_dGeom tests
// 	Vector3d f_r0, f_r1, f_r2, f_n;
// 	//df_dr0
// 	double r0[3] = {1.0,0.0,0.0};
// 	double r1[3] = {1.0,1.0,0.0};
// 	double r2[3] = {0.0,1.0,0.0};
// 	double n[3] = {0.0,0.0,1.0};
// 	// calc variations
// 	df_dgeom(r0,r1,r2,n,f_r0,f_r1,f_r2,f_n);
// 	double dr0[3] = {fRand(-1,1)/(10*NormTriad(r0)),
// 					  fRand(-1,1)/(10*NormTriad(r0)),
// 					  fRand(-1,1)/(10*NormTriad(r0))};
// 	double dr1[3] = {fRand(-1,1)/(10*NormTriad(r1)),
// 					  fRand(-1,1)/(10*NormTriad(r1)),
// 					  fRand(-1,1)/(10*NormTriad(r1))};
// 	double dr2[3] = {fRand(-1,1)/(10*NormTriad(r2)),
// 					  fRand(-1,1)/(10*NormTriad(r2)),
// 					  fRand(-1,1)/(10*NormTriad(r2))};
// 	double dn[3] = {fRand(-1,1)/(10*NormTriad(n)),
// 					  fRand(-1,1)/(10*NormTriad(n)),
// 					  fRand(-1,1)/(10*NormTriad(n))};

// 	// test vars in r0
// 	cout << "df_dr0: ---------------" << endl << "dr0:"
// 		 << endl << Vector3d(dr0[0],dr0[1],dr0[2])
// 		 << endl;
// 	double r0Test[3];
// 	AddTriad(r0,dr0,r0Test);
// 	double fGeomExactR0 = fGeom(r0Test,r1,r2,n);
// 	double fGeomAppDr0 = fGeom(r0,r1,r2,n);
// 	double fLinDeltaR0 = f_r0.transpose()*Vector3d(dr0[0],dr0[1],dr0[2]);
// 	cout << "error:" << endl << fGeomExactR0 - (fGeomAppDr0 + fLinDeltaR0)
// 	     << endl;

// 	// test vars in r1
// 	cout << "df_dr1: ---------------" << endl << "dr1:"
// 		 << endl << Vector3d(dr1[0],dr1[1],dr1[2])
// 		 << endl;
// 	double r1Test[3];
// 	AddTriad(r1,dr1,r1Test);
// 	double fGeomExactR1 = fGeom(r0,r1Test,r2,n);
// 	double fGeomAppDr1 = fGeom(r0,r1,r2,n);
// 	double fLinDeltaR1 = f_r1.transpose()*Vector3d(dr1[0],dr1[1],dr1[2]);
// 	cout << "error:" << endl << fGeomExactR1 - (fGeomAppDr1 + fLinDeltaR1)
// 		 << endl;

// 	// test vars in r2
// 	cout << "df_dr2: ---------------" << endl << "dr2:"
// 		 << endl << Vector3d(dr2[0],dr2[1],dr2[2])
// 		 << endl;
// 	double r2Test[3];
// 	AddTriad(r2,dr2,r2Test);
// 	double fGeomExactR2 = fGeom(r0,r1,r2Test,n);
// 	double fGeomAppDr2 = fGeom(r0,r1,r2,n);
// 	double fLinDeltaR2 = f_r2.transpose()*Vector3d(dr2[0],dr2[1],dr2[2]);
// 	cout << "error:" << endl << fGeomExactR2 - (fGeomAppDr2 + fLinDeltaR2)
// 		 << endl;

// 	// test vars in n
// 	cout << "df_dn: ----------------" << endl << "dn:"
// 		 << endl << Vector3d(dn[0],dn[1],dn[2])
// 		 << endl;
// 	double nTest[3];
// 	AddTriad(n,dn,nTest);
// 	double fGeomExactN = fGeom(r0,r1,r2,nTest);
// 	double fGeomAppDn = fGeom(r0,r1,r2,n);
// 	double fLinDeltaN = f_n.transpose()*Vector3d(dn[0],dn[1],dn[2]);
// 	cout << "error:" << endl << fGeomExactN - (fGeomAppDn + fLinDeltaN)
// 		 << endl;

// 	// test variations of the normal vector
// 	VectorXd normals(3*K);
// 	VectorXd approxNormals(3*K);
// 	VectorXd dNormals(3*K);
// 	getNormals(zeta,M,N,normals);
// 	getNormals(zetaPdel,M,N,approxNormals);
// 	for (int k = 0; k < K; k++) {
// 		// 1st diagonal
// 		Vector3d d = Vector3d(zeta[3*q_k(k,N,3)],
// 							  zeta[3*q_k(k,N,3)+1],
// 							  zeta[3*q_k(k,N,3)+2])
// 					-Vector3d(zeta[3*q_k(k,N,1)],
// 							  zeta[3*q_k(k,N,1)+1],
// 							  zeta[3*q_k(k,N,1)+2]);
// 		// 2nd diagonal (corner 4 to corner 2)
// 		Vector3d e = Vector3d(zeta[3*q_k(k,N,2)],
// 							  zeta[3*q_k(k,N,2)+1],
// 							  zeta[3*q_k(k,N,2)+2])
// 					-Vector3d(zeta[3*q_k(k,N,4)],
// 							  zeta[3*q_k(k,N,4)+1],
// 							  zeta[3*q_k(k,N,4)+2]);
// 		Matrix3d dN_dd = dn_dd(d,e);
// 		Matrix3d dN_de = dn_de(d,e);
// 		// calculate dNormals
// 		dNormals.block<3,1>(3*k,0) = dN_dd *
// 									 ( Vector3d(dZeta[3*q_k(k,N,3)],
// 									   		    dZeta[3*q_k(k,N,3)+1],
// 											    dZeta[3*q_k(k,N,3)+2])
// 									  -Vector3d(dZeta[3*q_k(k,N,1)],
// 											    dZeta[3*q_k(k,N,1)+1],
// 											    dZeta[3*q_k(k,N,1)+2]) ) +
// 									 dN_de *
// 									 ( Vector3d(dZeta[3*q_k(k,N,2)],
// 												dZeta[3*q_k(k,N,2)+1],
// 												dZeta[3*q_k(k,N,2)+2])
// 									  -Vector3d(dZeta[3*q_k(k,N,4)],
// 												dZeta[3*q_k(k,N,4)+1],
// 												dZeta[3*q_k(k,N,4)+2]) );
// 	}
// 	// add contribution of variations to approx result
// 	cout << "dn_dzeta: ----------------" << endl << "dZeta:"
// 			 << endl << dZeta
// 			 << endl;
// 	cout << "error in normals:" << endl << (dNormals - (approxNormals-normals))
// 	     << endl;

// 	// test A
// 	cout << endl << "A(zeta): ----------------" << endl;
// 	Eigen::MatrixXd A(K,K);
// 	AIC(zetaPtr,M,N,zetaPtr,M,N,false,A.data());
// =======
//	VectorXd dZeta(zeta.size()); dZeta(1) = 0.0123;
	VectorXd zetaPdel = zeta + dZeta;
	const double* zetaPdelPtr = zetaPdel.data();

//	genW(zetaPtr,M,N,W.data());
//
//	genAstar(zeta,M,N,zeta,M,N,Astar);
//
//	// test dAgam_dzeta subroutines
//	// dxHat_dx
//	Matrix3d dX;
//	Vector3d x(1.0,1.0,1.0);
//	Vector3d dx = Vector3d::Random()/(10.0*x.norm());
//	Vector3d f = (x + dx).normalized();
//	dX = dxHat_dx(x);
//	Vector3d fApprox = x.normalized() + dX * dx;
//	Vector3d diff = f-fApprox;
//	cout << "dxHat_dx: ---------------" << endl << "dx:" << endl << dx
//	     << endl << "norm (l-2) rel error:" << endl << diff.norm()/dx.norm()
//	     << endl
//	     << endl;
//
//	// Init matrices for df_dGeom tests
//	Vector3d f_r0, f_r1, f_r2, f_n;
//	//df_dr0
//	double r0[3] = {1.0,0.0,0.0};
//	double r1[3] = {1.0,1.0,0.0};
//	double r2[3] = {0.0,1.0,0.0};
//	double n[3] = {0.0,0.0,1.0};
//	// calc variations
//	df_dgeom(r0,r1,r2,n,f_r0,f_r1,f_r2,f_n);
//	double dr0[3] = {fRand(-1,1)/(10*NormTriad(r0)),
//					  fRand(-1,1)/(10*NormTriad(r0)),
//					  fRand(-1,1)/(10*NormTriad(r0))};
//	double dr1[3] = {fRand(-1,1)/(10*NormTriad(r1)),
//					  fRand(-1,1)/(10*NormTriad(r1)),
//					  fRand(-1,1)/(10*NormTriad(r1))};
//	double dr2[3] = {fRand(-1,1)/(10*NormTriad(r2)),
//					  fRand(-1,1)/(10*NormTriad(r2)),
//					  fRand(-1,1)/(10*NormTriad(r2))};
//	double dn[3] = {fRand(-1,1)/(10*NormTriad(n)),
//					  fRand(-1,1)/(10*NormTriad(n)),
//					  fRand(-1,1)/(10*NormTriad(n))};
//
//	// test vars in r0
//	cout << "df_dr0: ---------------" << endl << "dr0:"
//		 << endl << Vector3d(dr0[0],dr0[1],dr0[2])
//		 << endl;
//	double r0Test[3];
//	AddTriad(r0,dr0,r0Test);
//	double fGeomExactR0 = fGeom(r0Test,r1,r2,n);
//	double fGeomAppDr0 = fGeom(r0,r1,r2,n);
//	double fLinDeltaR0 = f_r0.transpose()*Vector3d(dr0[0],dr0[1],dr0[2]);
//	cout << "error:" << endl << fGeomExactR0 - (fGeomAppDr0 + fLinDeltaR0)
//	     << endl;
//
//	// test vars in r1
//	cout << "df_dr1: ---------------" << endl << "dr1:"
//		 << endl << Vector3d(dr1[0],dr1[1],dr1[2])
//		 << endl;
//	double r1Test[3];
//	AddTriad(r1,dr1,r1Test);
//	double fGeomExactR1 = fGeom(r0,r1Test,r2,n);
//	double fGeomAppDr1 = fGeom(r0,r1,r2,n);
//	double fLinDeltaR1 = f_r1.transpose()*Vector3d(dr1[0],dr1[1],dr1[2]);
//	cout << "error:" << endl << fGeomExactR1 - (fGeomAppDr1 + fLinDeltaR1)
//		 << endl;
//
//	// test vars in r2
//	cout << "df_dr2: ---------------" << endl << "dr2:"
//		 << endl << Vector3d(dr2[0],dr2[1],dr2[2])
//		 << endl;
//	double r2Test[3];
//	AddTriad(r2,dr2,r2Test);
//	double fGeomExactR2 = fGeom(r0,r1,r2Test,n);
//	double fGeomAppDr2 = fGeom(r0,r1,r2,n);
//	double fLinDeltaR2 = f_r2.transpose()*Vector3d(dr2[0],dr2[1],dr2[2]);
//	cout << "error:" << endl << fGeomExactR2 - (fGeomAppDr2 + fLinDeltaR2)
//		 << endl;
//
//	// test vars in n
//	cout << "df_dn: ----------------" << endl << "dn:"
//		 << endl << Vector3d(dn[0],dn[1],dn[2])
//		 << endl;
//	double nTest[3];
//	AddTriad(n,dn,nTest);
//	double fGeomPdelN = fGeom(r0,r1,r2,nTest);
//	double fGeomN = fGeom(r0,r1,r2,n);
//	double fLinDeltaN = f_n.transpose()*Vector3d(dn[0],dn[1],dn[2]);
//	cout << "fGeomPdelN:" << endl << fGeomPdelN << endl << "fGeomN:" << endl << fGeomN << endl;
//	cout << "ref + lin:" << endl << fGeomN + fLinDeltaN
//		 << endl;
//
//	// test variations of the normal vector
//	VectorXd normals(3*K);
//	VectorXd approxNormals(3*K);
//	VectorXd dNormals(3*K);
//	getNormals(zeta,M,N,normals);
//	getNormals(zetaPdel,M,N,approxNormals);
//	for (int k = 0; k < K; k++) {
//		// 1st diagonal
//		Vector3d d = Vector3d(zeta[3*q_k(k,N,3)],
//							  zeta[3*q_k(k,N,3)+1],
//							  zeta[3*q_k(k,N,3)+2])
//					-Vector3d(zeta[3*q_k(k,N,1)],
//							  zeta[3*q_k(k,N,1)+1],
//							  zeta[3*q_k(k,N,1)+2]);
//		// 2nd diagonal (corner 4 to corner 2)
//		Vector3d e = Vector3d(zeta[3*q_k(k,N,2)],
//							  zeta[3*q_k(k,N,2)+1],
//							  zeta[3*q_k(k,N,2)+2])
//					-Vector3d(zeta[3*q_k(k,N,4)],
//							  zeta[3*q_k(k,N,4)+1],
//							  zeta[3*q_k(k,N,4)+2]);
//		Matrix3d dN_dd = dn_dd(d,e);
//		Matrix3d dN_de = dn_de(d,e);
//		// calculate dNormals
//		dNormals.block<3,1>(3*k,0) = dN_dd *
//									 ( Vector3d(dZeta[3*q_k(k,N,3)],
//									   		    dZeta[3*q_k(k,N,3)+1],
//											    dZeta[3*q_k(k,N,3)+2])
//									  -Vector3d(dZeta[3*q_k(k,N,1)],
//											    dZeta[3*q_k(k,N,1)+1],
//											    dZeta[3*q_k(k,N,1)+2]) ) +
//									 dN_de *
//									 ( Vector3d(dZeta[3*q_k(k,N,2)],
//												dZeta[3*q_k(k,N,2)+1],
//												dZeta[3*q_k(k,N,2)+2])
//									  -Vector3d(dZeta[3*q_k(k,N,4)],
//												dZeta[3*q_k(k,N,4)+1],
//												dZeta[3*q_k(k,N,4)+2]) );
//	}
//	// add contribution of variations to approx result
//	cout << "dn_dzeta: ----------------" << endl << "dZeta:"
//			 << endl << dZeta
//			 << endl;
//	cout << "error in normals:" << endl << (dNormals - (approxNormals-normals))
//	     << endl;

	// test A
	cout << endl << "A(zeta): ----------------" << endl;
	EigDynMatrixRM A(K,K);
	AIC(zetaPtr,M,N,zetaPtr,M,N,true,A.data());
//>>>>>>> rob/master
	cout << A << endl;

	// test dAgamma_dZeta matrix
	cout << endl << "dAgamma_dZeta: ----------------" << endl;
//<<<<<<< HEAD
//	dAgamma0_dZeta(zeta.data(),M,N,gam.data(),zeta.data(),M,N,dAgam0_dzeta.data());
//	VectorXd wPdel(K);
//	MatrixXd Adel(K,K);
//	AIC(zetaPdelPtr,M,N,zetaPdelPtr,M,N,false,Adel.data());
//	wPdel = Adel*gam;
//	VectorXd w(K);
//	w = A*gam;
//	VectorXd wApprox(K);
//=======
	cout<<endl<<zeta+dZeta<<endl;
	dAgamma0_dZeta(zeta.data(),M,N,gam.data(),zeta.data(),M,N,true,dAgam0_dzeta.data());
	VectorXd wPdel(K); wPdel.setZero();
	EigDynMatrixRM Adel(K,K); Adel.setZero();
	AIC(zetaPdelPtr,M,N,zetaPdelPtr,M,N,true,Adel.data());
	cout << endl << "A(pDelZeta): --------" << endl << Adel << endl;
	wPdel = Adel*gam;
	VectorXd w(K); w.setZero();
	w = A*gam;
	VectorXd wApprox(K); wApprox.setZero();
//>>>>>>> rob/master
	wApprox = dAgam0_dzeta*dZeta;
	cout << "dw (approx):" << endl << wApprox << endl;
	cout << "dw (exact):" << endl << wPdel-w << endl;

//<<<<<<< HEAD
	// // test dWzetaPri0_dzeta
	// cout << endl << "dWzetaPri0_dzeta: ----------------" << endl;
	// // obtain linearization
	// VectorXd zetaPri0 = VectorXd::Random(zeta.size())*10;
	// const double* zetaPriPtr = zetaPri0.data();
	// EigDynMatrixRM dWzetaPri0_dzeta(K,3*K_zeta);
	// dWzetaPri0_dZeta(zetaPtr,M,N,zetaPriPtr,dWzetaPri0_dzeta.data());
	// VectorXd dwApprox = dWzetaPri0_dzeta*dZeta;
	// // check numerically
	// EigDynMatrixRM W1(K,3*K_zeta), W2(K,3*K_zeta);
	// genW(zetaPtr,M,N,W1.data());
	// genW(zetaPdelPtr,M,N,W2.data());
	// VectorXd dwExact = W2*zetaPri0 - W1*zetaPri0;
	// cout << "dw (approx):" << endl << dwApprox << endl;
	// cout << "dw (exact):" << endl << dwExact << endl;
	// cout << "dw (normed diff):" << endl << (dwApprox - dwExact).array()/(W1*zetaPri0).array() << endl;

	// // check elements of H matrix
	// const unsigned int M2=2;
	// const unsigned int N2=1;
	// const unsigned int K2=M2*N2;
	// const unsigned int K_zeta2=(M2+1)*(N2+1);
	// EigDynMatrixRM H = MatrixXd::Zero(3*K_zeta2,12*K2);
	// genH(M2,N2,H.data());
	// cout << endl << "H matrix: -------------------" << endl;
	// cout << endl << H << endl;

	// // print out 3 component AIC
	// EigDynMatrixRM aic3 = MatrixXd::Zero(3*K,K);
	// AIC3(zeta.data(),M,N,zeta.data(),M,N,aic3.data());
	// cout << endl << "AIC3 --------------------" << endl << aic3 << endl;

	// // print out dA3gamma0_dZeta matrix
	// EigDynMatrixRM dA3 = MatrixXd::Zero(3*K,3*K_zeta);
	// dA3gamma0_dZeta(zeta.data(),M,N,gam.data(),zeta.data(),M,N,dA3.data());
	// cout << endl << "dA3 -------------------------" << endl << endl << dA3 << endl;

	// // test numerically
	// EigDynMatrixRM aic3Pdel = MatrixXd::Zero(3*K,K);
	// AIC3(zetaPdelPtr,M,N,zetaPdelPtr,M,N,aic3Pdel.data());
	// dwExact = aic3Pdel*gam - aic3*gam;
	// dwApprox = dA3*dZeta;
	// cout << "dw (approx):" << endl << dwApprox << endl;
	// cout << "dw (exact):" << endl << dwExact << endl;
	// cout << "dw (normed diff):" << endl << (dwApprox - dwExact).array()/(aic3*gam).norm() << endl;
//=======
//	// test dAgamma_dZeta_num matrix
//	cout << endl << "dAgamma_dZeta_num: ----------------" << endl;
//	dAgamma0_dZeta_num(zeta.data(),M,N,gam.data(),zeta.data(),M,N,true,dAgam0_dzeta.data());
//	VectorXd wPdel(K); wPdel.setZero();
//	MatrixXd Adel(K,K); Adel.setZero();
//	AIC(zetaPdel.data(),M,N,zetaPdel.data(),M,N,true,Adel.data());
//	cout << endl << "A(pDelZeta): --------" << endl << Adel << endl;
//	wPdel = Adel*gam;
//	VectorXd w = A*gam;
//	VectorXd wApprox = dAgam0_dzeta*dZeta;
//	cout << "dAgam0_dzeta:" << endl << dAgam0_dzeta << endl;
//	cout << "dw (approx):" << endl << wApprox << endl;
//	cout << "dw (exact):" << endl << wPdel-w << endl;

//	// test dWzetaPri0_dzeta
//	cout << endl << "dWzetaPri0_dzeta: ----------------" << endl;
//	// obtain linearization
//	VectorXd zetaPri0 = VectorXd::Random(zeta.size())*10;
//	const double* zetaPriPtr = zetaPri0.data();
//	EigDynMatrixRM dWzetaPri0_dzeta(K,3*K_zeta);
//	dWzetaPri0_dZeta(zetaPtr,M,N,zetaPriPtr,dWzetaPri0_dzeta.data());
//	VectorXd dwApprox = dWzetaPri0_dzeta*dZeta;
//	// check numerically
//	EigDynMatrixRM W1(K,3*K_zeta), W2(K,3*K_zeta);
//	genW(zetaPtr,M,N,W1.data());
//	genW(zetaPdelPtr,M,N,W2.data());
//	VectorXd dwExact = W2*zetaPri0 - W1*zetaPri0;
//	cout << "dw (approx):" << endl << dwApprox << endl;
//	cout << "dw (exact):" << endl << dwExact << endl;
//	cout << "dw (normed diff):" << endl << (dwApprox - dwExact).array()/(W1*zetaPri0).array() << endl;

//	// check elements of H matrix
//	const unsigned int M2=2;
//	const unsigned int N2=1;
//	const unsigned int K2=M2*N2;
//	const unsigned int K_zeta2=(M2+1)*(N2+1);
//	EigDynMatrixRM H = MatrixXd::Zero(3*K_zeta2,12*K2);
//	genH(M2,N2,H.data());
//	cout << endl << "H matrix: -------------------" << endl;
//	cout << endl << H << endl;

	// print out 3 component AIC
//	EigDynMatrixRM aic3 = MatrixXd::Zero(3*K,K);
//	AIC3(zeta.data(),M,N,zeta.data(),M,N,aic3.data());
//	cout << endl << "AIC3 --------------------" << endl << aic3 << endl;
//
//	// print out dA3gamma0_dZeta matrix
//	EigDynMatrixRM dA3 = MatrixXd::Zero(3*K,3*K_zeta);
//	dA3gamma0_dZeta(zeta.data(),M,N,gam.data(),zeta.data(),M,N,dA3.data());
//	cout << endl << "dA3 -------------------------" << endl << endl << dA3 << endl;
//
//	// test numerically
//	EigDynMatrixRM aic3Pdel = MatrixXd::Zero(3*K,K);
//	AIC3(zetaPdelPtr,M,N,zetaPdelPtr,M,N,aic3Pdel.data());
//	dwExact = aic3Pdel*gam - aic3*gam;
//	dwApprox = dA3*dZeta;
//	cout << "dw (approx):" << endl << dwApprox << endl;
//	cout << "dw (exact):" << endl << dwExact << endl;
//	cout << "dw (normed diff):" << endl << (dwApprox - dwExact).array()/(aic3*gam).norm() << endl;

	// test dAs3gam0_dZeta_numerical
	EigDynMatrixRM aic3s(12*K,K);
	EigDynMatrixRM aic3sPdel(12*K,K);
	EigDynMatrixRM dA3sGam0_dZeta(12*K,3*K_zeta); dA3sGam0_dZeta.setZero();
	VectorXd v(12*K);
	VectorXd vPdel(12*K);
	VectorXd delVexact(12*K);
	VectorXd delVapprox(12*K);
	AIC3s(zetaPtr,M,N,zetaPtr,M,N,true,aic3s.data());
	AIC3s(zetaPdelPtr,M,N,zetaPdelPtr,M,N,true,aic3sPdel.data());
	v=aic3s*gam;
	vPdel=aic3sPdel*gam;
	delVexact=vPdel-v;
	dAs3gam0_dZeta_numerical(zetaPtr,M,N,gam.data(),zetaPtr,M,N,true,dA3sGam0_dZeta.data());
	delVapprox=dA3sGam0_dZeta*dZeta;

	cout << endl << "AIC3s --------------------" << endl << aic3s << endl;
	cout << endl << "AIC3sPdel --------------------" << endl << aic3sPdel << endl;
	cout << endl << "delVexact --------------------" << endl << delVexact << endl;
	cout << endl << "delVapprox --------------------" << endl << delVapprox << endl;
	cout << endl << "diff: --------------------" << endl << (delVapprox - delVexact).array()/delVexact.array() << endl;
//>>>>>>> rob/master
}
