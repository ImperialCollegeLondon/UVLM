/**@brief      AIC and downwash matrices from new UVLM description.
 * @author     Rob Simpson
 * @contact    r.simpson11@imperial.ac.uk
 * @version    0.0
 * @date       13/05/2015
 * @pre        None
 * @warning    None
 */

#include <aicMats.hpp>
#include <assert.h>
#include <indices.hpp>
#include <PanelTools.hpp>
#include <iostream>
#include <triads.hpp>
#include <stdexcept>
#include <vorticity.hpp>
using namespace Eigen;
using namespace std;

void genW(const double* zeta_,
		  const unsigned int M,
		  const unsigned int N,
		  double* W_) {
	/**@brief Populate W, downwash interpolation, matrix.
	 * @param zeta Vector of lattice vertices.
	 * @param M Chordwise panels.
	 * @param N Spanwise panels.
	 * @return W downwash interpolation matrix.
	 */

	// const cast Eigen outputs
	ConstMapVectXd zeta(zeta_,3*(M+1)*(N+1));
	EigenMapMatrixXd W(W_,M*N,3*(M+1)*(N+1));

	// temporary vars
	VectorXd normals = VectorXd::Zero(3*M*N);
	MatrixXd norMat = MatrixXd::Zero(M*N,3*M*N);
	MatrixXd Xi = MatrixXd::Zero(3*M*N,3*(M+1)*(N+1));
	Matrix3d XiKern = Matrix3d::Zero();

	// calculate normal vectors
	getNormals(zeta,M,N,normals);

	// get matrix of normal vectors
	unsigned int K=M*N;
	for (unsigned int k = 0; k < K; k++){
		norMat.block<1,3>(k,3*k)=normals.block<3,1>(3*k,0).transpose();
	}

	// get interpolation matrix, \Xi
	for (unsigned int k = 0; k < K; k++){
		for (unsigned int q = 0; q < W.cols()/3; q++){
			XiKernel(k,q,N,0.5,0.5,XiKern);
			Xi.block<3,3>(3*k,3*q)=XiKern;
		}
	}

	// calculate product of normal and interpolation matrices
	W = norMat*Xi;
	return;
}

void getNormals(const VectorXd& zeta,
		        const unsigned int M,
		        const unsigned int N,
		        const VectorXd& normals_) {
	/**@brief Get vector containing panel normal vectors.
	 * @param zeta Vector containing lattice vertices.
	 * @param M chordwise panels.
	 * @param N spanwise panels.
	 * @return normals Vector containing normal vectors.
	 */

	// const cast Eigen outputs
	VectorXd& normals = const_cast<VectorXd&> (normals_);

	// check inputs
	assert(zeta.size()==3*(M+1)*(N+1));
	assert(normals.size()==3*M*N);

	// declare corner points
	double x1[3] = {0.0,0.0,0.0};
	double x2[3] = {0.0,0.0,0.0};
	double x3[3] = {0.0,0.0,0.0};
	double x4[3] = {0.0,0.0,0.0};
	// declare normal vectors
	double n[3] = {0.0,0.0,0.0};

	unsigned int K = M*N;
	for (unsigned int k = 0; k < K; k++){

		// get corner points
		AssignTriad(x1,zeta(3*q_k(k,N,1)),zeta(3*q_k(k,N,1)+1),zeta(3*q_k(k,N,1)+2));
		AssignTriad(x2,zeta(3*q_k(k,N,2)),zeta(3*q_k(k,N,2)+1),zeta(3*q_k(k,N,2)+2));
		AssignTriad(x3,zeta(3*q_k(k,N,3)),zeta(3*q_k(k,N,3)+1),zeta(3*q_k(k,N,3)+2));
		AssignTriad(x4,zeta(3*q_k(k,N,4)),zeta(3*q_k(k,N,4)+1),zeta(3*q_k(k,N,4)+2));

		// calculate normal
		PanelNormal(x1,x2,x3,x4,n);

		//save to normals
		normals.block<3,1>(3*k,0)=Vector3d(n[0],n[1],n[2]);
	}
	return;
}

void genAstar(const VectorXd& zetaSrc,
				const unsigned int mSrc,
				const unsigned int nSrc,
		        const VectorXd& zetaTgt,
		        const unsigned int mTgt,
		        const unsigned int nTgt,
		        const EigDynMatrixRM& Astar_) {
	/**@brief Get vector containing panel normal vectors.
	 * @param zetaSrc Vector containing lattice vertices of source.
	 * @param mSrc Spanwise panels in the source lattice.
	 * @param nSrc Spanwise panels in the source lattice.
	 * @param zetaTgt Vector containing lattice vertices of target.
	 * @param mTgt Spanwise panels in the source lattice.
	 * @param nTgt Spanwise panels in the target lattice.
	 * @param N spanwise panels.
	 * @return Astar AIC matrix from source to target.
	 */

	// const cast Eigen outputs
	EigDynMatrixRM& Astar = const_cast<EigDynMatrixRM&> (Astar_);

	// reused panel numbers
	const unsigned int kSrc = mSrc*nSrc;
	const unsigned int kTgt_zeta = (mTgt+1)*(nTgt+1);

	// check sizes
	assert(zetaSrc.size()==3*(mSrc+1)*(nSrc+1));
	assert(zetaTgt.size()==3*kTgt_zeta);
	assert(Astar.rows()==3*kTgt_zeta);
	assert(Astar.cols()==kSrc);

	// init temps
	const double gamma = 1.0;
	unsigned int segStart = 0;
	unsigned int segEnd = 0;
	double xP[3] = {0.0, 0.0, 0.0};
	double x1[3] = {0.0, 0.0, 0.0};
	double x2[3] = {0.0, 0.0, 0.0};
	double vel[3] = {0.0, 0.0, 0.0};
	double sumVel[3] = {0.0, 0.0, 0.0};

	// loop through grid points, panels, segments
	for (unsigned int q = 0; q < kTgt_zeta; q++) {
		for (unsigned int k = 0; k < kSrc; k++) {
			// reset sum
			AssignTriad(sumVel,0.0,0.0,0.0);
			for (unsigned int l = 1; l < 5; l++) {
				// get segment corner points
				switch(l){
				case 1 :
					segStart=1;
					segEnd=2;
					break;
				case 2 :
					segStart=2;
					segEnd=3;
					break;
				case 3 :
					segStart=3;
					segEnd=4;
					break;
				case 4 :
					segStart=4;
					segEnd=1;
					break;
				default:
					throw std::invalid_argument("Segment index invalid.");
				}
				// get target point
				AssignTriad(xP,zetaTgt(3*q),zetaTgt(3*q+1),zetaTgt(3*q+2));
				// start of segment
				AssignTriad(x1,
						    zetaSrc(3*q_k(k,nSrc,segStart)),
						    zetaSrc(3*q_k(k,nSrc,segStart)+1),
						    zetaSrc(3*q_k(k,nSrc,segStart)+2));
				// end of segment
				AssignTriad(x2,
							zetaSrc(3*q_k(k,nSrc,segEnd)),
							zetaSrc(3*q_k(k,nSrc,segEnd)+1),
							zetaSrc(3*q_k(k,nSrc,segEnd)+2));
				// get velocity
				C_BiotSegment(xP,x1,x2,gamma,vel);
				// add to total
				AddTriad(sumVel,vel,sumVel);
			}
			// save to aStar
			Astar.block<3,1>(3*q,k)=Vector3d(sumVel[0],
													sumVel[1],
													sumVel[2]);
		}
	}
	return;
}

double fGeom(const double* r0,
          const double* r1,
	 	  const double* r2,
	 	  const double* n) {
	/**@brief Calculate geometry influence function f.
	 * @param r0 Vector r0.
	 * @param r1 Vector r1.
	 * @param r2 Vector r2.
	 * @param n Normal vector of target panel.
	 * @return Geometrical part of biot-Savart kernel.
	 */

	// temps
	double u[3] = {0.0,0.0,0.0};
	double r1hat[3] = {0.0,0.0,0.0}; // normalized r1
	double r2hat[3] = {0.0,0.0,0.0}; // normalized r2
	double normSub[3] = {0.0,0.0,0.0};
	double x[3] = {0.0,0.0,0.0};

	// get u
	CrossTriad(r1,r2,u);
//<<<<<<< HEAD
//	if (NormTriad(u) <= 1.0e-5) {
//=======
	if (NormTriad(u) <= 1.0e-8) {
//>>>>>>> rob/master
		return 0.0;
	} else {
		// get x
		MulTriad(u,1.0/pow(NormTriad(u),2.0),x);

		// get normalised triads
		NormaliseTriad(r1,r1hat);
		NormaliseTriad(r2,r2hat);

		// subtract them
		SubTriad(r1hat,r2hat,normSub);

		return DotTriad(r0,normSub)*DotTriad(x,n);
	}
}

void fGeom3(const double* r0,
          const double* r1,
	 	  const double* r2,
	 	  double* out) {
	/**@brief Calculate geometry influence function f.
	 * @param r0 Vector r0.
	 * @param r1 Vector r1.
	 * @param r2 Vector r2.
	 * @return out Output 'velocity' vector.
	 */

	// zero output
	out[0]=0.0;
	out[1]=0.0;
	out[2]=0.0;

	// temps
	double u[3] = {0.0,0.0,0.0};
	double r1hat[3] = {0.0,0.0,0.0}; // normalized r1
	double r2hat[3] = {0.0,0.0,0.0}; // normalized r2
	double normSub[3] = {0.0,0.0,0.0};
	double x[3] = {0.0,0.0,0.0};

	// get u
	CrossTriad(r1,r2,u);
// <<<<<<< HEAD
// 	if (NormTriad(u) <= 1.0e-5) {
// =======
	if (NormTriad(u) <= 1.0e-8) {
//>>>>>>> rob/master
		return;
	} else {
		// get x
		MulTriad(u,1.0/pow(NormTriad(u),2.0),x);

		// get normalised triads
		NormaliseTriad(r1,r1hat);
		NormaliseTriad(r2,r2hat);

		// subtract them
		SubTriad(r1hat,r2hat,normSub);

		// calculate output r0^T normSub * x.
		MulTriad(x,DotTriad(r0,normSub),out);
	}
	return;
}

void df_dgeom(const double* r0,
		        const double* r1,
			 	const double* r2,
			 	const double* n,
			 	const Vector3d& f_r0_,
			 	const Vector3d& f_r1_,
			 	const Vector3d& f_r2_,
			 	const Vector3d& f_n_) {
	/**@brief Calculate df/dr and df/dn, row vectors.
	 * @param r0 Vector r0.
     * @param r1 Vector r1.
     * @param r2 Vector r2.
     * @param n Normal vector of target panel.
     * @return f_r0 Row vector, gradient of f w.r.t r0.
     * @return f_r1 Row vector, gradient of f w.r.t r1.
     * @return f_r2 Row vector, gradient of f w.r.t r2.
     * @return f_n Row vector, gradient of f w.r.t n.
     */

	// const cast Eigen outputs
	Vector3d& f_r0 = const_cast<Vector3d&> (f_r0_);
	Vector3d& f_r1 = const_cast<Vector3d&> (f_r1_);
	Vector3d& f_r2 = const_cast<Vector3d&> (f_r2_);
	Vector3d& f_n = const_cast<Vector3d&> (f_n_);

	// temps
	double u[3] = {0.0,0.0,0.0};
	double r1hat[3] = {0.0,0.0,0.0}; // normalized r1
	double r2hat[3] = {0.0,0.0,0.0}; // normalized r2
	double normSub[3] = {0.0,0.0,0.0};
	double x[3] = {0.0,0.0,0.0};

	// get u
	CrossTriad(r1,r2,u);

	// get x
	MulTriad(u,1.0/pow(NormTriad(u),2.0),x);
	Vector3d xV(x[0],x[1],x[2]);

	// get normalised triads
	NormaliseTriad(r1,r1hat);
	NormaliseTriad(r2,r2hat);

	// sum them
	SubTriad(r1hat,r2hat,normSub);

	// if within cut-off radius, r-vector derivatives are set to zero
	if (NormTriad(u) <= 1.0e-5) {

		f_r0.setZero();
		f_r1.setZero();
		f_r2.setZero();

	} else {

		// temps
		double f_r0_triad[3] = {0.0,0.0,0.0};

		// get f_r0
		MulTriad(normSub,DotTriad(x,n),f_r0_triad);
		f_r0(0) = f_r0_triad[0];
		f_r0(1) = f_r0_triad[1];
		f_r0(2) = f_r0_triad[2];

		// get f_r1
		Vector3d r0v(r0[0],r0[1],r0[2]);
		Vector3d r1v(r1[0],r1[1],r1[2]);
		Vector3d r2v(r2[0],r2[1],r2[2]);
		Vector3d nV(n[0],n[1],n[2]);
		// calc dr1Hat_dr1 terms
		Vector3d preMul;
		preMul = DotTriad(x,n)*r0v;
		// calc d(u/|u|^2)/dr1 terms
		Vector3d preMul2;
		preMul2 = DotTriad(r0,normSub)*nV;
		f_r1 = preMul.transpose()*dxHat_dx(r1v)
			 + preMul2.transpose()*duHat2_dr1(r1v,r2v);

		// get f_r2
		f_r2 = -preMul.transpose()*dxHat_dx(r2v)
			 + preMul2.transpose()*duHat2_dr2(r1v,r2v);
	}

	// get f_n
	f_n = DotTriad(r0,normSub)*xV; // f_n needs to be transposed

	return;
}

void df3_dgeom(const double* r0,
		        const double* r1,
			 	const double* r2,
			 	const Matrix3d& f3_r0_,
			 	const Matrix3d& f3_r1_,
			 	const Matrix3d& f3_r2_) {
	/**@brief Calculate df/dr and df/dn, row vectors.
	 * @param r0 Vector r0.
     * @param r1 Vector r1.
     * @param r2 Vector r2.
     * @return f3_r0 Matrix, gradient of f3 w.r.t r0.
     * @return f3_r1 Matrix, gradient of f3 w.r.t r1.
     * @return f3_r2 Matrix, gradient of f3 w.r.t r2.
     */

	// const cast Eigen outputs
	Matrix3d& f3_r0 = const_cast<Matrix3d&> (f3_r0_);
	Matrix3d& f3_r1 = const_cast<Matrix3d&> (f3_r1_);
	Matrix3d& f3_r2 = const_cast<Matrix3d&> (f3_r2_);

	// temps
	double u[3] = {0.0,0.0,0.0};
	double r1hat[3] = {0.0,0.0,0.0}; // normalized r1
	double r2hat[3] = {0.0,0.0,0.0}; // normalized r2
	double normSub[3] = {0.0,0.0,0.0};
	double x[3] = {0.0,0.0,0.0};

	// get u
	CrossTriad(r1,r2,u);

	// get x
	MulTriad(u,1.0/pow(NormTriad(u),2.0),x);

	// get normalised triads
	NormaliseTriad(r1,r1hat);
	NormaliseTriad(r2,r2hat);

	// sum them
	SubTriad(r1hat,r2hat,normSub);

	// create Eigen maps
	ConstMapVect3d xV(x);
	ConstMapVect3d rHatV(normSub);
	ConstMapVect3d r0v(r0);
	ConstMapVect3d r1v(r1);
	ConstMapVect3d r2v(r2);

	// if within cut-off radius, r-vector derivatives are set to zero
	if (NormTriad(u) <= 1.0e-5) {

		f3_r0.setZero();
		f3_r1.setZero();
		f3_r2.setZero();

	} else {
		// get f3_r0
		f3_r0 = xV*rHatV.transpose();

		// get f_r1
		f3_r1 = r0v.dot(rHatV) * duHat2_dr1(r1v,r2v)
				+ xV*r0v.transpose() * dxHat_dx(r1v);

		// get f_r2
		f3_r2 = r0v.dot(rHatV) * duHat2_dr2(r1v,r2v)
				- xV*r0v.transpose() * dxHat_dx(r2v);
	}
	return;
}

Matrix3d dxHat_dx(const Vector3d& x) {
	/**@brief Calculate derivative of x/|x| w.r.t x.
	 * @param x Vector x.
	 * @return dX Matrix output.
	 */

	// init
	Matrix3d dX;

	if (x.norm() < 1.0e-5) {
		dX.setZero();
	} else {
		Matrix3d Eye = Matrix3d::Identity();
		dX = (1.0/x.norm())*(Eye - (1.0/pow(x.norm(),2.0)) * x*x.transpose());
	}
	return dX;
}

Matrix3d duHat2_dr1(const Vector3d& r1,
				     const Vector3d& r2) {
	/**@brief Calculate derivative of (r1xr2)/|r1xr2|^2 w.r.t r1.
	 * @param r1 Vector r1.
	 * @param r2 Vector r2.
	 * @return dX Matrix output.
	 */

	// init
	Matrix3d dX;

	// temps
	Vector3d u = r1.cross(r2);

	if (u.norm() < 1.0e-5) {
		dX.setZero();
	} else {
		Matrix3d Eye = Matrix3d::Identity();
		dX = - (1.0/pow(u.norm(),2.0))
			 * (Eye - (2.0/pow(u.norm(),2.0)) * u*u.transpose())
			 * skew(r2);
	}
	return dX;
}

Matrix3d duHat2_dr2(const Vector3d& r1,
				     const Vector3d& r2) {
	/**@brief Calculate derivative of (r1xr2)/|r1xr2|^2 w.r.t r2.
	 * @param r1 Vector r1.
	 * @param r2 Vector r2.
	 * @return dX Matrix output.
	 */

	// init
	Matrix3d dX;

	// temps
	Vector3d u = r1.cross(r2);

	if (u.norm() < 1.0e-5) {
		dX.setZero();
	} else {
		Matrix3d Eye = Matrix3d::Identity();
		dX =   (1.0/pow(u.norm(),2.0))
			 * (Eye - (2.0/pow(u.norm(),2.0)) * u*u.transpose())
			 * skew(r1);
	}
	return dX;
}

Matrix3d skew(const Vector3d& x) {
	/**@brief Calculate skew-symmetric matrix of vector x.
	 */
	Matrix3d X;
	X  << 0.0,  -x(2), x(1),
		  x(2),  0.0, -x(0),
		 -x(1),  x(0), 0.0;
	return X;
}

Matrix3d dn_dd(const Vector3d& d, const Vector3d& e) {
	/**@brief Calculate derivative of normal vector w.r.t d (diagonal 1).
	 * @param d Diagonal vector 1 (corner 1 to corner 3).
	 * @param e Diagonal vector 2 (corner 4 to corner 2).
	 * @return dX Matrix output.
	 */
	return -dxHat_dx(d.cross(e))*skew(e);
}

Matrix3d dn_de(const Vector3d& d, const Vector3d& e) {
	/**@brief Calculate derivative of normal vector w.r.t e (diagonal 2).
	 * @param d Diagonal vector 1 (corner 1 to corner 3).
	 * @param e Diagonal vector 2 (corner 4 to corner 2).
	 * @return dX Matrix output.
	 */
	return dxHat_dx(d.cross(e))*skew(d);
}

void dAgamma0_dZeta(const double* zetaSrc_,
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
					 double* dX_) {
	/**@brief Calculate tensor-free derivative of (A gamma_0) w.r.t zeta.
	 * @param zetaSrc Grid points of source lattice.
	 * @param mSrc chordwise panels on source lattice.
	 * @param nSrc spanwise panels on source lattice.
	 * @param gamma0 Reference circulation distribution on source lattice.
	 * @param zetaTgt Grid points of target lattice.
	 * @param mTgt chordwise panels on target lattice.
	 * @param nTgt spanwise panels on target lattice.
//<<<<<<< HEAD
//=======
	 * @param imageMethod Include influence of vorticity across x-plane.
//>>>>>>> rob/master
	 * @return dX K x 3K_{\zeta_{tgt}} matrix output.
	 * @warning If the zeta arguments are the same they must be the same object,
	 * therefore in the function call the arguments must be previously
	 * instantiated objects, e.g (zeta+delZeta, ..., zeta+delZeta, ...) is
	 * invalid because a two temps are instantiated.
	 */

	// eigen map output matrix
	ConstMapVectXd zetaSrc(zetaSrc_,3*(mSrc+1)*(nSrc+1));
	ConstMapVectXd gamma0(gamma0_,mSrc*nSrc);
	ConstMapVectXd zetaTgt(zetaTgt_,3*(mTgt+1)*(nTgt+1));
	EigenMapMatrixXd dX(dX_,mTgt*nTgt,3*(mTgt+1)*(nTgt+1));

	if (gamma0.isZero() == true) {
		return;
	}

	// Set dX to zero
	dX.setZero();

	// temps
	unsigned int kSrc = mSrc*nSrc;
	unsigned int kTgt = mTgt*nTgt;
	unsigned int qTgt = (mTgt+1)*(nTgt+1);
	unsigned int ll = 0; //segment counter
	unsigned int llp1 = 0; //segment counter
//<<<<<<< HEAD
//	Vector3d c1, c2, c3, c4, cp; // panel corner points, collocation point
//	Vector3d d, e, n; // panel diagonal and normal vectors
//	Matrix3d n_d, n_e;
//	Vector3d r0, r1, r2; // Biot-Savart kernel vectors
//	Vector3d f_r0, f_r1, f_r2, f_n; //dvtvs required for target/source variation
//	Matrix3d Xi; // interpolating matrix
//	double a; // prefactor (\gamma_0(k2)/(4 \pi)).
//=======
	Vector3d c1 = Vector3d::Zero();
	Vector3d c2 = Vector3d::Zero();
	Vector3d c3 = Vector3d::Zero();
	Vector3d c4 = Vector3d::Zero();
	Vector3d cp = Vector3d::Zero(); // panel corner points, collocation point
	Vector3d d  = Vector3d::Zero();
	Vector3d e = Vector3d::Zero();
	Vector3d n = Vector3d::Zero(); // panel diagonal and normal vectors
	Matrix3d n_d = Matrix3d::Zero();
	Matrix3d n_e = Matrix3d::Zero();
	Vector3d r0 = Vector3d::Zero();
	Vector3d r1 = Vector3d::Zero();
	Vector3d r2 = Vector3d::Zero(); // Biot-Savart kernel vectors
	Vector3d f_r0 = Vector3d::Zero();
	Vector3d f_r1 = Vector3d::Zero();
	Vector3d f_r2 = Vector3d::Zero();
	Vector3d f_n  = Vector3d::Zero(); //dvtvs required for target/source variation
	Matrix3d Xi = Matrix3d::Zero(); // interpolating matrix
	double a; // prefactor (\gamma_0(k2)/(4 \pi)).
	Vector3d x1 = Vector3d::Zero(); // Segment start/end
	Vector3d x2 = Vector3d::Zero();
//>>>>>>> rob/master

	// loop through DoFs to make (1x3) submatrices
	for (unsigned int k1 = 0; k1 < kTgt; k1++) {
		// calc n, dn_dd, dn_de, colloc point only once for each target panel
		c1 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,1),0);
		c2 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,2),0);
		c3 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,3),0);
		c4 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,4),0);
		// diagonals
		d = c3 - c1;
		e = c2 - c4;
		// calc n
		n = (d.cross(e)).normalized();
		// calc dn_dd, dn_de
		n_d = dn_dd(d,e);
		n_e = dn_de(d,e);
		// collocation point
		BilinearInterpTriad(c1.data(),
							c2.data(),
							c3.data(),
							c4.data(),
							cp.data(),
							0.5,0.5,false);
		for (unsigned int k2 = 0; k2 < kSrc; k2++) {
			if (gamma0(k2) == 0.0) {
				continue;
			}
			a = gamma0(k2)/(4*M_PI);
			for (unsigned int q = 0; q < qTgt; q++) {
				for (unsigned int l = 1; l < 5; l++) {
					// roll around segment index
					if (l < 4) {
						ll = l;
						llp1 = l+1;
					} else if (l == 4) {
						ll = l;
						llp1= 1;
					}
					// check if either k1 or k2 are active based on q
					if (q_k(k1,nTgt,1) == q || q_k(k1,nTgt,2) == q ||
						q_k(k1,nTgt,3) == q || q_k(k1,nTgt,4) == q ||
						q_k(k2,nSrc,ll) == q || q_k(k2,nSrc,llp1) == q) {

						// contributions at targets
//<<<<<<< HEAD
//						// calc r0, r1, r2
//						r0 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0)
//							-zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
//						// r1
//						r1 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
//						// r2
//						r2 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0);
//=======
						// segment endpoints
						x1 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
						x2 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0);
						// calc r0, r1, r2
						r0 = x2 - x1;
						// r1
						r1 = cp - x1;
						// r2
						r2 = cp - x2;
//>>>>>>> rob/master

						// calc f_r0, f_r1, f_r2, f_n
						df_dgeom(r0.data(),
								 r1.data(),
								 r2.data(),
								 n.data(),
								 f_r0,
								 f_r1,
								 f_r2,
								 f_n);

						// add effect to output matrix
						XiKernel(k1,q,nTgt,0.5,0.5,Xi);

						if (q_k(k1,nTgt,1) == q) {
							// add Xi, and Kronecker delta terms, c1
							dX.block<1,3>(k1,3*q) += a*(
								(f_r1 + f_r2).transpose()*Xi
							   - f_n.transpose()*n_d);
						} else if (q_k(k1,nTgt,2) == q) {
							// add Xi, and Kronecker delta terms, c2
							dX.block<1,3>(k1,3*q) += a*(
								(f_r1 + f_r2).transpose()*Xi
							   + f_n.transpose()*n_e);
						} else if (q_k(k1,nTgt,3) == q) {
							// add Xi, and Kronecker delta terms, c3
							dX.block<1,3>(k1,3*q) += a*(
								(f_r1 + f_r2).transpose()*Xi
							   + f_n.transpose()*n_d);
						} else if (q_k(k1,nTgt,4) == q) {
							// add Xi, and Kronecker delta terms, c4
							dX.block<1,3>(k1,3*q) += a*(
								(f_r1 + f_r2).transpose()*Xi
							   - f_n.transpose()*n_e);
						} // end if k1,q (targets)

						// contributions from sources
						if (zetaSrc.data() == zetaTgt.data()) {
							if (q_k(k2,nSrc,ll) == q) {
								// add Kronecker delta term, segment start
								dX.block<1,3>(k1,3*q) += -a*(f_r0+f_r1).transpose();
							} else if (q_k(k2,nSrc,llp1) == q) {
								// add Kronecker delta term, segment end
								dX.block<1,3>(k1,3*q) += a*(f_r0-f_r2).transpose();
							} // end if k2,q (sources)
						} // end if Src == Tgt

//<<<<<<< HEAD
//=======
						if (imageMeth == true) {
							r0(1)=-r0(1);
							x1(1)=-x1(1);
							x2(1)=-x2(1);
							r1=cp-x1;
							r2=cp-x2;
							df_dgeom(r0.data(),
									 r1.data(),
									 r2.data(),
									 n.data(),
									 f_r0,
									 f_r1,
									 f_r2,
									 f_n);
							if (q_k(k1,nTgt,1) == q) {
								// add Xi, and Kronecker delta terms, c1
								dX.block<1,3>(k1,3*q) += a*(
									(f_r1 + f_r2).transpose()*Xi
								   - f_n.transpose()*n_d);
							} else if (q_k(k1,nTgt,2) == q) {
								// add Xi, and Kronecker delta terms, c2
								dX.block<1,3>(k1,3*q) += a*(
									(f_r1 + f_r2).transpose()*Xi
								   + f_n.transpose()*n_e);
							} else if (q_k(k1,nTgt,3) == q) {
								// add Xi, and Kronecker delta terms, c3
								dX.block<1,3>(k1,3*q) += a*(
									(f_r1 + f_r2).transpose()*Xi
								   + f_n.transpose()*n_d);
							} else if (q_k(k1,nTgt,4) == q) {
								// add Xi, and Kronecker delta terms, c4
								dX.block<1,3>(k1,3*q) += a*(
									(f_r1 + f_r2).transpose()*Xi
								   - f_n.transpose()*n_e);
							} // end if k1,q (targets)
							// contributions from sources
							if (zetaSrc.data() == zetaTgt.data()) {
								if (q_k(k2,nSrc,ll) == q) {
									// add Kronecker delta term, segment start
									dX.block<1,3>(k1,3*q) += -a*(f_r0+f_r1).transpose();
								} else if (q_k(k2,nSrc,llp1) == q) {
									// add Kronecker delta term, segment end
									dX.block<1,3>(k1,3*q) += a*(f_r0-f_r2).transpose();
								} // end if k2,q (sources)
							} // end if Src == Tgt
						} // end if imageMeth

//>>>>>>> rob/master
					} else {
						continue;
					} // end if k,q
				} // for l
			} // for q
		} // for k2
	} // for k1
	return;
}

//<<<<<<< HEAD
//=======
void dAgamma0_dZeta_num(const double* zetaSrc_,
					 	   const unsigned int mSrc,
					 	   const unsigned int nSrc,
					 	   const double* gamma0_,
					 	   const double* zetaTgt_,
					 	   const unsigned int mTgt,
					 	   const unsigned int nTgt,
					 	   const bool imageMeth,
					 	   double* dX_) {
	/**@brief Calculate tensor-free derivative of (A gamma_0) w.r.t zeta numerically.
	 * @param zetaSrc Grid points of source lattice.
	 * @param mSrc chordwise panels on source lattice.
	 * @param nSrc spanwise panels on source lattice.
	 * @param gamma0 Reference circulation distribution on source lattice.
	 * @param zetaTgt Grid points of target lattice.
	 * @param mTgt chordwise panels on target lattice.
	 * @param nTgt spanwise panels on target lattice.
	 * @param imageMethod Include influence of vorticity across x-plane.
	 * @return dX K x 3K_{\zeta_{tgt}} matrix output.
	 * @warning If the zeta arguments are the same they must be the same object,
	 * therefore in the function call the arguments must be previously
	 * instantiated objects, e.g (zeta+delZeta, ..., zeta+delZeta, ...) is
	 * invalid because a two temps are instantiated.
	 */

	// eigen map output matrix
	ConstMapVectXd zetaSrc(zetaSrc_,3*(mSrc+1)*(nSrc+1));
	ConstMapVectXd gamma0(gamma0_,mSrc*nSrc);
	ConstMapVectXd zetaTgt(zetaTgt_,3*(mTgt+1)*(nTgt+1));
	EigenMapMatrixXd dX(dX_,mTgt*nTgt,3*(mTgt+1)*(nTgt+1));

	if (gamma0.isZero() == true) {
		return;
	}

	// Set dX to zero
	dX.setZero();

	// temps
	const unsigned int qTgt = 3*(mTgt+1)*(nTgt+1);
	const double del = 0.00001;
	VectorXd delZeta(3*(mTgt+1)*(nTgt+1)); delZeta.setZero();
	VectorXd zetaPdel(3*(mTgt+1)*(nTgt+1)); zetaPdel.setZero();

	// unperturbed downwash
	MatrixXd AIC0(mTgt*nTgt,mSrc*nSrc); AIC0.setZero();
	AIC(zetaSrc_,mSrc,nSrc,zetaTgt_,mTgt,nTgt,imageMeth,AIC0.data());
	VectorXd wRef(mTgt*nTgt); wRef.setZero();
	wRef=AIC0*gamma0;

	// pertubed entities
	VectorXd wDel(mTgt*nTgt); wDel.setZero();
	MatrixXd AIC_del(mTgt*nTgt,mSrc*nSrc); AIC_del.setZero();

	for (unsigned int q = 0; q < qTgt; q++) {
		delZeta.setZero();
		delZeta(q)=del;
		zetaPdel=zetaTgt+delZeta;
		if (zetaTgt_ == zetaSrc_) {
			AIC(zetaPdel.data(),mSrc,nSrc,zetaPdel.data(),mTgt,nTgt,imageMeth,AIC_del.data());
		} else {
			AIC(zetaSrc_,mSrc,nSrc,zetaPdel.data(),mTgt,nTgt,imageMeth,AIC_del.data());
		}
		wDel = (AIC_del*gamma0)-wRef;
		dX.col(q)=wDel/del;
	}
	return;
}

//>>>>>>> rob/master
void AIC(const double* zetaSrc_,
		  const unsigned int mSrc,
		  const unsigned int nSrc,
		  const double* zetaTgt_,
		  const unsigned int mTgt,
		  const unsigned int nTgt,
		  const bool imageMeth,
		  double* dX_) {
	/**@brief Calculate AIC matrix.
	 * @param zetaSrc Grid points of source lattice.
	 * @param mSrc chordwise panels on source lattice.
	 * @param nSrc spanwise panels on source lattice.
	 * @param zetaTgt Grid points of target lattice.
	 * @param mTgt chordwise panels on target lattice.
	 * @param nTgt spanwise panels on target lattice.
//<<<<<<< HEAD
//	 * @param imageMethod Include influence of vorticity across y-plane.
//=======
	 * @param imageMethod Include influence of vorticity across x-plane.
//>>>>>>> rob/master
	 * @return dX K_tgt x K_src matrix output.
	 */

	// Create Eigen maps to memory
	ConstMapVectXd zetaSrc(zetaSrc_,3*(mSrc+1)*(nSrc+1));
	ConstMapVectXd zetaTgt(zetaTgt_,3*(mTgt+1)*(nTgt+1));
	EigenMapMatrixXd dX(dX_,mTgt*nTgt,mSrc*nSrc);

	// Set dX to zero
	dX.setZero();

	// temps
	unsigned int kSrc = mSrc*nSrc;
	unsigned int kTgt = mTgt*nTgt;
	unsigned int ll = 0; //segment counter
	unsigned int llp1 = 0; //segment counter
	Vector3d c1 = Vector3d::Zero();// panel corner points, collocation point
	Vector3d c2 = Vector3d::Zero();
	Vector3d c3 = Vector3d::Zero();
	Vector3d c4 = Vector3d::Zero();
	Vector3d cp = Vector3d::Zero();
	Vector3d d = Vector3d::Zero(); // panel diagonal and normal vectors
	Vector3d e = Vector3d::Zero();
	Vector3d n = Vector3d::Zero();
	Vector3d r0 = Vector3d::Zero(); // Biot-Savart kernel vectors
	Vector3d r1 = Vector3d::Zero();
	Vector3d r2 = Vector3d::Zero();
	Vector3d x1 = Vector3d::Zero(); // Segment start/end
	Vector3d x2 = Vector3d::Zero();

	for (unsigned int k1 = 0; k1 < kTgt; k1++) {
		// calc n, colloc point only once for each target panel
		c1 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,1),0);
		c2 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,2),0);
		c3 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,3),0);
		c4 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,4),0);
		// diagonals
		d = c3 - c1;
		e = c2 - c4;
		// calc n
		n = (d.cross(e)).normalized();
		// collocation point
		BilinearInterpTriad(c1.data(),
							c2.data(),
							c3.data(),
							c4.data(),
							cp.data(),
							0.5,0.5,false);
		for (unsigned int k2 = 0; k2 < kSrc; k2++) {
			for (unsigned int l = 1; l < 5; l++) {
				// roll around segment index
				if (l < 4) {
					ll = l;
					llp1 = l+1;
				} else if (l == 4) {
					ll = l;
					llp1= 1;
				}

				// segment points
				x1 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
				x2 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0);

				// calc r0
				r0 = x2 - x1;
				// r1
				r1 = cp - x1;
				// r2
				r2 = cp - x2;
				// AIC entry
				dX(k1,k2)+=1.0/(4.0*M_PI)*fGeom(r0.data(),
										    r1.data(),
										    r2.data(),
										    n.data());

				if (imageMeth == true) {
					r0(1)=-r0(1);
					x1(1)=-x1(1);
					x2(1)=-x2(1);
					r1=cp-x1;
					r2=cp-x2;
					dX(k1,k2)+=-1.0/(4.0*M_PI)*fGeom(r0.data(),
													r1.data(),
													r2.data(),
													n.data());
				}
			}
		}
	}
	return;
}

//<<<<<<< HEAD
//=======
//TODO: Image method in functions below, plus testing.

//>>>>>>> rob/master
void dA3gamma0_dZeta(const double* zetaSrc_,
					 const unsigned int mSrc,
					 const unsigned int nSrc,
					 const double* gamma0_,
					 const double* zetaTgt_,
					 const unsigned int mTgt,
					 const unsigned int nTgt,
					 double* dX_) {
	/**@brief Calculate tensor-free derivative of (A^c gamma_0) w.r.t zeta.
	 * @param zetaSrc Grid points of source lattice.
	 * @param mSrc chordwise panels on source lattice.
	 * @param nSrc spanwise panels on source lattice.
	 * @param gamma0 Reference circulation distribution on source lattice.
	 * @param zetaTgt Grid points of target lattice.
	 * @param mTgt chordwise panels on target lattice.
	 * @param nTgt spanwise panels on target lattice.
	 * @return dX 3K x 3K_{\zeta_{tgt}} matrix output.
	 * @warning If the zeta arguments are the same they must be the same object,
	 * therefore in the function call the arguments must be previously
	 * instantiated objects, e.g (zeta+delZeta, ..., zeta+delZeta, ...) is
	 * invalid because a two temps are instantiated.
	 */

	// eigen map output matrix
	ConstMapVectXd zetaSrc(zetaSrc_,3*(mSrc+1)*(nSrc+1));
	ConstMapVectXd gamma0(gamma0_,mSrc*nSrc);
	ConstMapVectXd zetaTgt(zetaTgt_,3*(mTgt+1)*(nTgt+1));
	EigenMapMatrixXd dX(dX_,3*mTgt*nTgt,3*(mTgt+1)*(nTgt+1));

//<<<<<<< HEAD
//=======
	if (gamma0.isZero() == true) {
		return;
	}

//>>>>>>> rob/master
	// Set dX to zero
	dX.setZero();

	// temps
	unsigned int kSrc = mSrc*nSrc;
	unsigned int kTgt = mTgt*nTgt;
	unsigned int qTgt = (mTgt+1)*(nTgt+1);
	unsigned int ll = 0; //segment counter
	unsigned int llp1 = 0; //segment counter
	Vector3d c1 = Vector3d::Zero(); // panel corner points, collocation point
	Vector3d c2 = Vector3d::Zero();
	Vector3d c3 = Vector3d::Zero();
	Vector3d c4 = Vector3d::Zero();
	Vector3d cp = Vector3d::Zero();
	Vector3d r0 = Vector3d::Zero(); // Biot-Savart kernel vectors
	Vector3d r1 = Vector3d::Zero();
	Vector3d r2 = Vector3d::Zero();
	Matrix3d f3_r0 = Matrix3d::Zero(); //drvtvs. required for target/source variation
	Matrix3d f3_r1 = Matrix3d::Zero();
	Matrix3d f3_r2 = Matrix3d::Zero();
	Matrix3d Xi = Matrix3d::Zero(); // interpolating matrix
	double a = 0.0; // prefactor (\gamma_0(k2)/(4 \pi)).

//<<<<<<< HEAD
//	// loop through DoFs to make (1x3) submatrices
//=======
	// loop through DoFs to make (3x3) submatrices
//>>>>>>> rob/master
	for (unsigned int k1 = 0; k1 < kTgt; k1++) {
		// calc n, dn_dd, dn_de, colloc point only once for each target panel
		c1 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,1),0);
		c2 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,2),0);
		c3 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,3),0);
		c4 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,4),0);
		// collocation point
		BilinearInterpTriad(c1.data(),
							c2.data(),
							c3.data(),
							c4.data(),
							cp.data(),
							0.5,0.5,false);
		for (unsigned int k2 = 0; k2 < kSrc; k2++) {
			a = gamma0(k2)/(4*M_PI);
			for (unsigned int q = 0; q < qTgt; q++) {
				for (unsigned int l = 1; l < 5; l++) {
					// roll around segment index
					if (l < 4) {
						ll = l;
						llp1 = l+1;
					} else if (l == 4) {
						ll = l;
						llp1= 1;
					}
					// check if either k1 or k2 are active based on q
					if (q_k(k1,nTgt,1) == q || q_k(k1,nTgt,2) == q ||
						q_k(k1,nTgt,3) == q || q_k(k1,nTgt,4) == q ||
						q_k(k2,nSrc,ll) == q || q_k(k2,nSrc,llp1) == q) {

						// contributions at targets
						// calc r0, r1, r2
						r0 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0)
							-zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
						// r1
						r1 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
						// r2
						r2 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0);

						// calc f_r0, f_r1, f_r2, f_n
						df3_dgeom(r0.data(),
								 r1.data(),
								 r2.data(),
								 f3_r0,
								 f3_r1,
								 f3_r2);

						// add effect to output matrix
						XiKernel(k1,q,nTgt,0.5,0.5,Xi);

						if (q_k(k1,nTgt,1) == q || q_k(k1,nTgt,2) == q ||
							q_k(k1,nTgt,3) == q || q_k(k1,nTgt,4) == q) {
							// add Xi terms at corners c1-c4
							dX.block<3,3>(3*k1,3*q) += a*(f3_r1 + f3_r2)*Xi;
						} // end if k1,q (targets)

						// contributions from sources
						if (zetaSrc.data() == zetaTgt.data()) {
							if (q_k(k2,nSrc,ll) == q) {
								// add Kronecker delta term, segment start
								dX.block<3,3>(3*k1,3*q) += -a*(f3_r0+f3_r1);
							} else if (q_k(k2,nSrc,llp1) == q) {
								// add Kronecker delta term, segment end
								dX.block<3,3>(3*k1,3*q) += a*(f3_r0-f3_r2);
							} // end if k2,q (sources)
						} // end if Src == Tgt

					} else {
						continue;
					} // end if k,q
				} // for l
			} // for q
		} // for k2
	} // for k1
	return;
}

void AIC3(const double* zetaSrc_,
		  const unsigned int mSrc,
		  const unsigned int nSrc,
		  const double* zetaTgt_,
		  const unsigned int mTgt,
		  const unsigned int nTgt,
		  double* dX_) {
	/**@brief Calculate AIC matrix (3 components of velocity).
	 * @param zetaSrc Grid points of source lattice.
	 * @param mSrc chordwise panels on source lattice.
	 * @param nSrc spanwise panels on source lattice.
	 * @param zetaTgt Grid points of target lattice.
	 * @param mTgt chordwise panels on target lattice.
	 * @param nTgt spanwise panels on target lattice.
	 * @return dX 3*K_tgt x K_src matrix output.
	 */

	// Create Eigen maps to memory
	ConstMapVectXd zetaSrc(zetaSrc_,3*(mSrc+1)*(nSrc+1));
	ConstMapVectXd zetaTgt(zetaTgt_,3*(mTgt+1)*(nTgt+1));
	EigenMapMatrixXd dX(dX_,3*mTgt*nTgt,mSrc*nSrc);
	// Set dX to zero
	dX.setZero();

	// temps
	unsigned int kSrc = mSrc*nSrc;
	unsigned int kTgt = mTgt*nTgt;
	unsigned int ll = 0; //segment counter
	unsigned int llp1 = 0; //segment counter
	Vector3d c1, c2, c3, c4, cp; // panel corner points, collocation point
	Vector3d r0, r1, r2, v; // Biot-Savart kernel vectors, velocity.

	for (unsigned int k1 = 0; k1 < kTgt; k1++) {
		// calc n, colloc point only once for each target panel
		c1 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,1),0);
		c2 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,2),0);
		c3 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,3),0);
		c4 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,4),0);
		// collocation point
		BilinearInterpTriad(c1.data(),
							c2.data(),
							c3.data(),
							c4.data(),
							cp.data(),
							0.5,0.5,false);
		for (unsigned int k2 = 0; k2 < kSrc; k2++) {
			for (unsigned int l = 1; l < 5; l++) {
				// roll around segment index
				if (l < 4) {
					ll = l;
					llp1 = l+1;
				} else if (l == 4) {
					ll = l;
					llp1= 1;
				}

				// calc r0
				r0 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0)
					-zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
				// r1
				r1 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
				// r2
				r2 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0);
				// AIC entry (3 components)
				fGeom3(r0.data(),r1.data(),r2.data(),v.data());
				dX.block<3,1>(3*k1,k2) += 1.0/(4.0*M_PI)*v;
			}
		}
	}
	return;
}

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
		  	double* dX_) {
	/**@brief Calculate AIC matrix (3 components of velocity) at the segment
	 * midpoints.
	 * @param zetaSrc Grid points of source lattice.
	 * @param mSrc chordwise panels on source lattice.
	 * @param nSrc spanwise panels on source lattice.
	 * @param zetaTgt Grid points of target lattice.
	 * @param mTgt chordwise panels on target lattice.
	 * @param nTgt spanwise panels on target lattice.
//<<<<<<< HEAD
//=======
//	 * @param imageMeth use image method accross x-z plane.
//>>>>>>> rob/master
	 * @return dX 12*K_tgt x K_src matrix output.
	 */

	// Create Eigen maps to memory
	ConstMapVectXd zetaSrc(zetaSrc_,3*(mSrc+1)*(nSrc+1));
	ConstMapVectXd zetaTgt(zetaTgt_,3*(mTgt+1)*(nTgt+1));
	EigenMapMatrixXd dX(dX_,12*mTgt*nTgt,mSrc*nSrc);
	// Set dX to zero
	dX.setZero();

	// temps
	unsigned int kSrc = mSrc*nSrc;
	unsigned int kTgt = mTgt*nTgt;
	unsigned int ll_s = 0; //segment counter
	unsigned int llp1_s = 0; //segment counter
	unsigned int ll_t = 0; //segment counter
	unsigned int llp1_t = 0; //segment counter
	unsigned int s = 0; // total segment index (at target midpoints)
//<<<<<<< HEAD
//=======
	Vector3d x1 = Vector3d::Zero();
	Vector3d x2 = Vector3d::Zero();
//>>>>>>> rob/master
	Vector3d r0  = Vector3d::Zero(); //Biot-Savart kernel vectors
	Vector3d r1 = Vector3d::Zero();
	Vector3d r2 = Vector3d::Zero();
	Vector3d v = Vector3d::Zero(); // velocity.

	for (unsigned int k1 = 0; k1 < kTgt; k1++) {
		for (unsigned int lt = 1; lt < 5; lt++) {
			if (lt < 4) {
				ll_t = lt;
				llp1_t = lt+1;
			} else if (lt == 4) {
				ll_t = lt;
				llp1_t = 1;
			}
			for (unsigned int k2 = 0; k2 < kSrc; k2++) {
				for (unsigned int ls = 1; ls < 5; ls++) {
					if (ls < 4) {
						ll_s = ls;
						llp1_s = ls+1;
					} else if (ls == 4) {
						ll_s = ls;
						llp1_s = 1;
					}
//<<<<<<< HEAD
//					// calc r0
//					r0 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1_s),0)
//						-zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll_s),0);
//					// r1
//					r1 = 0.5*(  zetaTgt.block<3,1>(3*q_k(k1,nTgt,llp1_t),0)
//							  + zetaTgt.block<3,1>(3*q_k(k1,nTgt,ll_t),0)  )
//						 -zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll_s),0);
//					// r1
//					r2 = 0.5*(  zetaTgt.block<3,1>(3*q_k(k1,nTgt,llp1_t),0)
//							  + zetaTgt.block<3,1>(3*q_k(k1,nTgt,ll_t),0)  )
//						 -zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1_s),0);
//					// AIC entry (3 components)
//					fGeom3(r0.data(),r1.data(),r2.data(),v.data());
//					dX.block<3,1>(3*s,k2) += 1.0/(4.0*M_PI)*v;
//=======
					// segment endpoints
					x1 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll_s),0);
					x2 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1_s),0);
					// calc r0
					r0 = x2 - x1;
					// r1
					r1 = 0.5*(  zetaTgt.block<3,1>(3*q_k(k1,nTgt,llp1_t),0)
							  + zetaTgt.block<3,1>(3*q_k(k1,nTgt,ll_t),0)  )
						 -x1;
					// r1
					r2 = 0.5*(  zetaTgt.block<3,1>(3*q_k(k1,nTgt,llp1_t),0)
							  + zetaTgt.block<3,1>(3*q_k(k1,nTgt,ll_t),0)  )
						 -x2;
					// AIC entry (3 components)
					fGeom3(r0.data(),r1.data(),r2.data(),v.data());
					dX.block<3,1>(3*s,k2) += 1.0/(4.0*M_PI)*v;

					if (imageMeth == true) {
						r0(1)=-r0(1);
						x1(1)=-x1(1);
						x2(1)=-x2(1);
						r1=0.5*(  zetaTgt.block<3,1>(3*q_k(k1,nTgt,llp1_t),0)
								  + zetaTgt.block<3,1>(3*q_k(k1,nTgt,ll_t),0)  )
							 -x1;
						r2=0.5*(  zetaTgt.block<3,1>(3*q_k(k1,nTgt,llp1_t),0)
								  + zetaTgt.block<3,1>(3*q_k(k1,nTgt,ll_t),0)  )
							 -x2;
						fGeom3(r0.data(),
							   r1.data(),
							   r2.data(),
							   v.data());
						dX.block<3,1>(3*s,k2) += -1.0/(4.0*M_PI)*v;
					}
//>>>>>>> rob/master
				}
			}
			s++;
		}
	}
	return;
}

void AIC3s_noTE(const double* zetaSrc_,
		  	const unsigned int mSrc,
		  	const unsigned int nSrc,
		  	const double* zetaTgt_,
		  	const unsigned int mTgt,
		  	const unsigned int nTgt,
		  	const bool wakeSrc,
		  	double* dX_) {
	/**@brief Calculate AIC matrix (3 components of velocity) at the segment
	 * midpoints.
	 * @param zetaSrc Grid points of source lattice.
	 * @param mSrc chordwise panels on source lattice.
	 * @param nSrc spanwise panels on source lattice.
	 * @param zetaTgt Grid points of target lattice.
	 * @param mTgt chordwise panels on target lattice.
	 * @param nTgt spanwise panels on target lattice.
	 * @return dX 12*K_tgt x K_src matrix output.
	 */

	// Create Eigen maps to memory
	ConstMapVectXd zetaSrc(zetaSrc_,3*(mSrc+1)*(nSrc+1));
	ConstMapVectXd zetaTgt(zetaTgt_,3*(mTgt+1)*(nTgt+1));
	EigenMapMatrixXd dX(dX_,12*mTgt*nTgt,mSrc*nSrc);
	// Set dX to zero
	dX.setZero();

	// temps
	unsigned int kSrc = mSrc*nSrc;
	unsigned int kTgt = mTgt*nTgt;
	unsigned int ll_s = 0; //segment counter
	unsigned int llp1_s = 0; //segment counter
	unsigned int ll_t = 0; //segment counter
	unsigned int llp1_t = 0; //segment counter
	unsigned int s = 0; // total segment index (at target midpoints)
	unsigned int mmSrc = 0;
	Vector3d r0  = Vector3d::Zero(); //Biot-Savart kernel vectors
	Vector3d r1 = Vector3d::Zero();
	Vector3d r2 = Vector3d::Zero();
	Vector3d v = Vector3d::Zero(); // velocity.

	for (unsigned int k1 = 0; k1 < kTgt; k1++) {
		for (unsigned int lt = 1; lt < 5; lt++) {
			if (lt < 4) {
				ll_t = lt;
				llp1_t = lt+1;
			} else if (lt == 4) {
				ll_t = lt;
				llp1_t = 1;
			}
			for (unsigned int k2 = 0; k2 < kSrc; k2++) {
				mmSrc = k2/nSrc;
				for (unsigned int ls = 1; ls < 5; ls++) {
					if (   (mmSrc == mSrc-1 && ls == 3 && wakeSrc == false)
						|| (mmSrc == 0      && ls == 1 && wakeSrc == true)  ){
						continue; // if segment is at TE
					}
					if (ls < 4) {
						ll_s = ls;
						llp1_s = ls+1;
					} else if (ls == 4) {
						ll_s = ls;
						llp1_s = 1;
					}
					// calc r0
					r0 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1_s),0)
						-zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll_s),0);
					// r1
					r1 = 0.5*(  zetaTgt.block<3,1>(3*q_k(k1,nTgt,llp1_t),0)
							  + zetaTgt.block<3,1>(3*q_k(k1,nTgt,ll_t),0)  )
						 -zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll_s),0);
					// r1
					r2 = 0.5*(  zetaTgt.block<3,1>(3*q_k(k1,nTgt,llp1_t),0)
							  + zetaTgt.block<3,1>(3*q_k(k1,nTgt,ll_t),0)  )
						 -zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1_s),0);
					// AIC entry (3 components)
					fGeom3(r0.data(),r1.data(),r2.data(),v.data());
					dX.block<3,1>(3*s,k2) += 1.0/(4.0*M_PI)*v;
				}
			}
			s++;
		}
	}
	return;
}

void AIC3noTE(const double* zetaSrc_,
		  	  const unsigned int mSrc,
		  	  const unsigned int nSrc,
		  	  const double* zetaTgt_,
		  	  const unsigned int mTgt,
		  	  const unsigned int nTgt,
		  	  const bool wakeSrc,
		  	  double* dX_) {
	/**@brief Calculate AIC matrix (3 components of velocity).
	 * @param zetaSrc Grid points of source lattice.
	 * @param mSrc chordwise panels on source lattice.
	 * @param nSrc spanwise panels on source lattice.
	 * @param zetaTgt Grid points of target lattice.
	 * @param mTgt chordwise panels on target lattice.
	 * @param nTgt spanwise panels on target lattice.
	 * @param wakeSrc The source is a wake, TE induced vels omitted.
	 * @return dX 3*K_tgt x K_src matrix output.
	 * @note the velocity induced by vorticity at the TE is omitted.
	 */

	// Create Eigen maps to memory
	ConstMapVectXd zetaSrc(zetaSrc_,3*(mSrc+1)*(nSrc+1));
	ConstMapVectXd zetaTgt(zetaTgt_,3*(mTgt+1)*(nTgt+1));
	EigenMapMatrixXd dX(dX_,3*mTgt*nTgt,mSrc*nSrc);
	// Set dX to zero
	dX.setZero();

	// temps
	unsigned int kSrc = mSrc*nSrc;
	unsigned int kTgt = mTgt*nTgt;
	unsigned int ll = 0; //segment counter
	unsigned int llp1 = 0; //segment counter
	Vector3d c1, c2, c3, c4, cp; // panel corner points, collocation point
	Vector3d r0, r1, r2, v; // Biot-Savart kernel vectors, velocity.
	unsigned int mmSrc = 0;

	for (unsigned int k1 = 0; k1 < kTgt; k1++) {
		// calc n, colloc point only once for each target panel
		c1 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,1),0);
		c2 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,2),0);
		c3 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,3),0);
		c4 = zetaTgt.block<3,1>(3*q_k(k1,nTgt,4),0);
		// collocation point
		BilinearInterpTriad(c1.data(),
							c2.data(),
							c3.data(),
							c4.data(),
							cp.data(),
							0.5,0.5,false);
		for (unsigned int k2 = 0; k2 < kSrc; k2++) {
			mmSrc = k2/nSrc;
			for (unsigned int l = 1; l < 5; l++) {
				if (   (mmSrc == mSrc-1 && l == 3 && wakeSrc == false)
					|| (mmSrc == 0      && l == 1 && wakeSrc == true)  ){
					continue; // if segment is at TE
				}
				// roll around segment index
				if (l < 4) {
					ll = l;
					llp1 = l+1;
				} else if (l == 4) {
					ll = l;
					llp1= 1;
				}

				// calc r0
				r0 = zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0)
					-zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
				// r1
				r1 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,ll),0);
				// r2
				r2 = cp - zetaSrc.block<3,1>(3*q_k(k2,nSrc,llp1),0);
				// AIC entry (3 components)
				fGeom3(r0.data(),r1.data(),r2.data(),v.data());
				dX.block<3,1>(3*k1,k2) += 1.0/(4.0*M_PI)*v;
			}
		}
	}
	return;
}

void dWzetaPri0_dZeta(const double* zeta_,
						const unsigned int m,
						const unsigned int n,
						const double* zetaPri_,
						double* dX_) {
	/**@brief Calculate dWzetaPri0_dZeta matrix.
	 * @param zeta Lattice coordinates.
	 * @param m Chordwise panels.
	 * @param n Spanwise panels.
	 * @param zetaPri Lattice velocities.
	 * @return dX output.
	 */

	ConstMapVectXd zeta(zeta_,3*(m+1)*(n+1));
	ConstMapVectXd zetaPri(zetaPri_,3*(m+1)*(n+1));
	EigenMapMatrixXd dX(dX_,m*n,3*(m+1)*(n+1));

	// Set Eigen memory to zero just in case
	dX.setZero();

	// temps
	unsigned int K=m*n;
	unsigned int Q1=(m+1)*(n+1);
	unsigned int Q2=(m+1)*(n+1);
	Vector3d zetaPriQ1, b, d, e;
	Matrix3d Xi, dbHat_db;

	for (unsigned int k = 0; k < K; k++) {
		d = zeta.block<3,1>(3*q_k(k,n,3),0)
		   -zeta.block<3,1>(3*q_k(k,n,1),0);
		e = zeta.block<3,1>(3*q_k(k,n,2),0)
		   -zeta.block<3,1>(3*q_k(k,n,4),0);
		b = d.cross(e);
		dbHat_db = dxHat_dx(b);
		for (unsigned int q1 = 0; q1 < Q1; q1++) {
			if (q_k(k,n,1) == q1 || q_k(k,n,2) == q1 ||
				q_k(k,n,3) == q1 || q_k(k,n,4) == q1) {
				// get vars
				zetaPriQ1 = zetaPri.block<3,1>(3*q1,0);
				XiKernel(k,q1,n,0.5,0.5,Xi);
				for (unsigned int q2 = 0; q2 < Q2; q2++) {
					//check if q2 is active based on k and corner point
					if (q_k(k,n,1) == q2) {
						// movement of corner 1
						dX.block<1,3>(k,3*q2) += -zetaPriQ1.transpose()
								             *Xi.transpose()
								             *dbHat_db
								             *skew(-e);
					} else if (q_k(k,n,2) == q2) {
						// movement of corner 2
						dX.block<1,3>(k,3*q2) += zetaPriQ1.transpose()
								             *Xi.transpose()
								             *dbHat_db
								             *skew(d);
					} else if (q_k(k,n,3) == q2) {
						// movement of corner 3
						dX.block<1,3>(k,3*q2) += zetaPriQ1.transpose()
											 *Xi.transpose()
											 *dbHat_db
											 *skew(-e);
					} else if (q_k(k,n,4) == q2) {
						// movement of corner 4
						dX.block<1,3>(k,3*q2) += -zetaPriQ1.transpose()
											 *Xi.transpose()
											 *dbHat_db
											 *skew(d);
					}
				}
			}
		}
	}
	return;
}

void genH(const unsigned int m,
		   const unsigned int n,
		   double* H_){
	/**@brief Calculate the segment to lattice vertex distribution matrix, H.
	 * @param m Chordwise panels.
	 * @param n Spanwise panels.
	 * @return H Segment to lattice vertex distribution matrix.
	 */

	// temps
	unsigned int k = m*n;
	unsigned int kZeta = (m+1)*(n+1);
	unsigned int kk, ll;
	Matrix3d hKern = Matrix3d::Zero();

	// create Eigen map to memory
	EigenMapMatrixXd H(H_,3*kZeta,12*k);

	for (unsigned int q = 0; q < kZeta; q++) {
		for (unsigned int s = 0; s < 4*k; s++) {
			// infer kk, ll
			kk = s/4; // panel index
			ll = s%4 + 1; // sub-panel segments are numbered 1,2,3,4
			hKernel(q,kk,ll,n,hKern);
			H.block<3,3>(3*q,3*s) = hKern;
		}
	}
	return;
}

void genXi(const unsigned int m,
		    const unsigned int n,
		    const double eta1,
		    const double eta2,
		    double* Xi_) {
	/**@brief Calculate the interpolation matrix, Xi.
	 * @param m Chordwise panels.
	 * @param n Spanwise panels.
	 * @param eta1 Computational coord in chordwise sense.
	 * @param eta2 computational coord in spanwise sense.
	 * @return Xi Lattice to sub-panel interpolation/distribution matrix.
	 */

	// map Eigen types
	EigenMapMatrixXd Xi(Xi_,3*m*n,3*(m+1)*(n+1));

	//temp
	Matrix3d xiKern = Matrix3d::Zero();

	for (unsigned int k = 0; k < m*n; k++) {
		for (unsigned int q = 0; q < (m+1)*(n+1); q++) {
			XiKernel(k,q,n,eta1,eta2,xiKern);
			Xi.block<3,3>(3*k,3*q) = xiKern;
		}
	}
}

void Y1(const double* vM_,
		const double* zeta_,
		const unsigned int m,
		const unsigned int n,
		double* Y1_) {
	/**@brief Calculate velocity-segment cross-product matrix for dGamma, Y1.
	 * @param vM 3-component fluid-grid relative vel at the collocation points.
	 * @param zeta Bound lattice vertices.
	 * @param m Chordwise panels.
	 * @param n Spanwise panels.
	 * @return Y1 output.
	 * @note If at the trailing edge the direct effect of dGamma is set zero.
	 */

	// map Eigen types
	ConstMapVectXd vM(vM_,12*m*n);
	ConstMapVectXd zeta(zeta_,3*(m+1)*(n+1));
	EigenMapMatrixXd Y1(Y1_,12*m*n,m*n);

	//temps
	unsigned int mm = 0;
	Vector3d c1 = Vector3d::Zero();
	Vector3d c2 = Vector3d::Zero();
	Vector3d c3 = Vector3d::Zero();
	Vector3d c4 = Vector3d::Zero();
	VectorXd yKern = VectorXd::Zero(12);

	// loop through panels
	for (unsigned int kk = 0; kk < m*n; kk++) {
		// infer mm
		mm = kk/n;
		// corner points
		c1 = zeta.block<3,1>(3*q_k(kk,n,1),0);
		c2 = zeta.block<3,1>(3*q_k(kk,n,2),0);
		c3 = zeta.block<3,1>(3*q_k(kk,n,3),0);
		c4 = zeta.block<3,1>(3*q_k(kk,n,4),0);
		// kernel
		yKern.block<3,1>(0,0) = vM.block<3,1>(12*kk,0).cross(c2-c1);
		yKern.block<3,1>(3,0) = vM.block<3,1>(12*kk+3,0).cross(c3-c2);
		if (mm < m-1) { // at TE
			yKern.block<3,1>(6,0) = vM.block<3,1>(12*kk+6,0).cross(c4-c3);
		} else {
			yKern.block<3,1>(6,0) = Vector3d::Zero();
		}
		yKern.block<3,1>(9,0) = vM.block<3,1>(12*kk+9,0).cross(c1-c4);
		// add to full matrix
		Y1.block<12,1>(12*kk,kk)=yKern;
	}
}

void Y2(const double* gamma_,
		const double* vM_,
		const unsigned int m,
		const unsigned int n,
		double* Y2_) {
	/**@brief Calculate force due to segment length vector variation (dZeta), Y2.
	 * @param gamma Panel circulation strengths.
	 * @param vM 3-component fluid-grid relative vel at the collocation points.
	 * @param m Chordwise panels.
	 * @param n Spanwise panels.
	 * @return Y2 output.
	 */

	//map Eigen types
	ConstMapVectXd gamma(gamma_,m*n);
	ConstMapVectXd vM(vM_,12*m*n);
	EigenMapMatrixXd Y2(Y2_,12*m*n,3*(m+1)*(n+1));

	// set to zero
	Y2.setZero();

	//temps
	Matrix3d yKern = Matrix3d::Zero();
	unsigned int ss = 0;
	unsigned int mm = 0;

	// loop through panels
	for (unsigned int kk = 0; kk < m*n; kk++) {
		// loop through segments
		mm = kk/n;
		for (unsigned int ll = 1; ll < 5; ll++) {
			yKern = gamma(kk)*skew(vM.block<3,1>(3*ss,0));
			for (unsigned int qq = 0; qq < (m+1)*(n+1); qq++) {
				if (ll < 4) {
					if (mm == m-1 && ll == 3) {
						continue; // at TE
					} else {
						if (qq == q_k(kk,n,ll)) {
							Y2.block<3,3>(3*ss,3*qq) = -yKern;
						} else if (qq == q_k(kk,n,ll+1)) {
							Y2.block<3,3>(3*ss,3*qq) = yKern;
						}
					}
				} else {
					if (qq == q_k(kk,n,4)) {
						Y2.block<3,3>(3*ss,3*qq) = -yKern;
					} else if (qq == q_k(kk,n,1)) {
						Y2.block<3,3>(3*ss,3*qq) = yKern;
					}
				}
			}
			ss++;
		}
	}
}

void Y3(const double* gamma_,
		const double* zeta_,
		const unsigned int m,
		const unsigned int n,
		double* Y3_) {
	/**@brief Calculate gamma *skew(segment) matrix for dvM, Y3.
	 * @param gamma Panel circulation strengths.
	 * @param zeta Bound lattice vertices.
	 * @param m Chordwise panels.
	 * @param n Spanwise panels.
	 * @return Y3 Output.
<<<<<<< HEAD
	 * @note If at the trailing edge the direct effect of dGamma is set zero.
=======
	 * @note If at the trailing edge the direct effect is set zero.
>>>>>>> rob/master
	 */

	// map Eigen types
	ConstMapVectXd gamma(gamma_,m*n);
	ConstMapVectXd zeta(zeta_,3*(m+1)*(n+1));
	EigenMapMatrixXd Y3(Y3_,12*m*n,12*m*n);

	// temps
	unsigned int mm = 0;
	Vector3d c1 = Vector3d::Zero();
	Vector3d c2 = Vector3d::Zero();
	Vector3d c3 = Vector3d::Zero();
	Vector3d c4 = Vector3d::Zero();
	MatrixXd yKern = MatrixXd::Zero(12,12);

	// loop through panels
	for (unsigned int kk = 0; kk < m*n; kk++) {
		// infer mm
		mm = kk/n;
		// corner points
		c1 = zeta.block<3,1>(3*q_k(kk,n,1),0);
		c2 = zeta.block<3,1>(3*q_k(kk,n,2),0);
		c3 = zeta.block<3,1>(3*q_k(kk,n,3),0);
		c4 = zeta.block<3,1>(3*q_k(kk,n,4),0);
		// kernel
		yKern.block<3,3>(0,0) = gamma(kk)*skew(c2-c1);
		yKern.block<3,3>(3,3) = gamma(kk)*skew(c3-c2);
		if (mm < m-1) {
			yKern.block<3,3>(6,6) = gamma(kk)*skew(c4-c3);
		} else { // at TE
			yKern.block<3,3>(6,6) = Matrix3d::Zero();
		}
		yKern.block<3,3>(9,9) = gamma(kk)*skew(c1-c4);
		// add to full matrix
		Y3.block<12,12>(12*kk,12*kk)=yKern;
	}
}

void Y4(const double* zeta_,
		const unsigned int m,
		const unsigned int n,
		double* Y4_) {
	/**@brief Calculate A_k n_k matrix for dGammaPrime, Y4.
	 * @param zeta Bound lattice vertices.
	 * @param m Chordwise panels.
	 * @param n Spanwise panels.
	 * @return Y4 Output.
	 */

	// map Eigen types
	ConstMapVectXd zeta(zeta_,3*(m+1)*(n+1));
	EigenMapMatrixXd Y4(Y4_,3*m*n,m*n);

	//temps
	Vector3d d = Vector3d::Zero();
	Vector3d e = Vector3d::Zero();

	// loop through panels
	for (unsigned int kk = 0; kk < m*n; kk++) {
		d = zeta.block<3,1>(3*q_k(kk,n,3),0)
				-zeta.block<3,1>(3*q_k(kk,n,1),0);
		e = zeta.block<3,1>(3*q_k(kk,n,2),0)
			    -zeta.block<3,1>(3*q_k(kk,n,4),0);
		Y4.block<3,1>(3*kk,kk) = 0.5*d.cross(e);
	}
}

void Y5(const double* gammaPri_,
		const double* zeta_,
		const unsigned int m,
		const unsigned int n,
		double* Y5_) {
	/**@brief Calculate force due to A_k n_k variation (dZeta), Y5.
	 * @param gammaPri Rate of change of panel circulation strengths.
	 * @param zeta Lattice vertices.
	 * @param m Chordwise panels.
	 * @param n Spanwise panels.
	 * @return Y5 output.
	 */

	// map Eigen types
	ConstMapVectXd gammaPri(gammaPri_,m*n);
	ConstMapVectXd zeta(zeta_,3*(m+1)*(n+1));
	EigenMapMatrixXd Y5(Y5_,3*m*n,3*(m+1)*(n+1));

	//temps
	Vector3d d = Vector3d::Zero();
	Vector3d e = Vector3d::Zero();
	Matrix3d dFac = Matrix3d::Zero();
	Matrix3d eFac = Matrix3d::Zero();

	// loop through panels
	for (unsigned int kk = 0; kk < m*n; kk++) {
		d = zeta.block<3,1>(3*q_k(kk,n,3),0)
				-zeta.block<3,1>(3*q_k(kk,n,1),0);
		e = zeta.block<3,1>(3*q_k(kk,n,2),0)
				-zeta.block<3,1>(3*q_k(kk,n,4),0);
		dFac=0.5*gammaPri(kk)*skew(d);
		eFac=0.5*gammaPri(kk)*skew(e);
		for (unsigned int qq = 0; qq < (m+1)*(n+1); qq++) {
			if (qq == q_k(kk,n,1)) {
				Y5.block<3,3>(3*kk,3*qq) = eFac;
			} else if (qq == q_k(kk,n,2)) {
				Y5.block<3,3>(3*kk,3*qq) = dFac;
			} else if (qq == q_k(kk,n,3)) {
				Y5.block<3,3>(3*kk,3*qq) = -eFac;
			} else if (qq == q_k(kk,n,4)) {
				Y5.block<3,3>(3*kk,3*qq) = -dFac;
			}
		}
	}
}

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
								double* dX_) {
	/**@brief Numerically calculate tensor-free derivative of (A^s gamma_0) w.r.t zeta.
	* @param zetaSrc Grid points of source lattice.
	* @param mSrc chordwise panels on source lattice.
	* @param nSrc spanwise panels on source lattice.
	* @param gamma0 Reference circulation distribution on source lattice.
	* @param zetaTgt Grid points of target lattice.
	* @param mTgt chordwise panels on target lattice.
	* @param nTgt spanwise panels on target lattice.
<<<<<<< HEAD
=======
	* @param imageMeth image method across x-z plane.
>>>>>>> rob/master
	* @return dX 3K x 3K_{\zeta_{tgt}} matrix output.
	* @warning If the zeta arguments are the same they must be the same object,
	* therefore in the function call the arguments must be previously
	* instantiated objects, e.g (zeta+delZeta, ..., zeta+delZeta, ...) is
	* invalid because a two temps are instantiated.
	*/

	// eigen maps
	ConstMapVectXd zetaSrc(zetaSrc_,3*(mSrc+1)*(nSrc+1));
	ConstMapVectXd gamma0(gamma0_,mSrc*nSrc);
	ConstMapVectXd zetaTgt(zetaTgt_,3*(mTgt+1)*(nTgt+1));
	EigenMapMatrixXd dX(dX_,12*mTgt*nTgt,3*(mTgt+1)*(nTgt+1));

	// check if gamma0 is all zeros
	if (gamma0.isZero() == true) {
		return;
	}


	// set to zero
	dX.setZero();

	// temps
	const unsigned int S = 12*mTgt*nTgt;
//<<<<<<< HEAD
//	const double del = 0.00001; // small perturbation
//=======
	const double del = 0.001; // small perturbation
//>>>>>>> rob/master
	VectorXd delZeta = VectorXd::Zero(3*(mTgt+1)*(nTgt+1));
	VectorXd zetaPdel = VectorXd::Zero(3*(mTgt+1)*(nTgt+1));
	VectorXd u = VectorXd::Zero(S);
	VectorXd dU = VectorXd::Zero(S);
	EigDynMatrixRM aic3s = EigDynMatrixRM::Zero(S,mSrc*nSrc);
	EigDynMatrixRM aic3sDel = EigDynMatrixRM::Zero(S,mSrc*nSrc);

	// calculate reference velocities
//	if (zetaSrc.data() == zetaTgt.data()) {
//		AIC3s_noTE(zetaSrc_,mSrc,nSrc,zetaTgt_,mTgt,nTgt,false,aic3s.data());
//	} else {
//		AIC3s_noTE(zetaSrc_,mSrc,nSrc,zetaTgt_,mTgt,nTgt,true,aic3s.data());
//	}

//<<<<<<< HEAD
//	AIC3s(zetaSrc_,mSrc,nSrc,zetaTgt_,mTgt,nTgt,aic3s.data());
//=======
	AIC3s(zetaSrc_,mSrc,nSrc,zetaTgt_,mTgt,nTgt,imageMeth,aic3s.data());
//>>>>>>> rob/master
	u = aic3s*gamma0;

	for (unsigned int qPri = 0; qPri < 3*(mTgt+1)*(nTgt+1); qPri++) {

		// set zero and permutate delZeta
		delZeta.setZero();
		delZeta[qPri] = del;

		// new aic
		if (zetaSrc.data() == zetaTgt.data()) {
			zetaPdel = zetaSrc + delZeta; // src and target lattices are same body
//<<<<<<< HEAD
//			AIC3s(zetaPdel.data(),mSrc,nSrc,zetaPdel.data(),mTgt,nTgt,aic3sDel.data());
//   //			AIC3s_noTE(zetaPdel.data(),mSrc,nSrc,zetaPdel.data(),mTgt,nTgt,false,aic3sDel.data());
//		} else {
//			zetaPdel = zetaTgt + delZeta; // src lattice is wake
//			AIC3s(zetaSrc_,mSrc,nSrc,zetaPdel.data(),mTgt,nTgt,aic3sDel.data());
//=======
			AIC3s(zetaPdel.data(),mSrc,nSrc,zetaPdel.data(),mTgt,nTgt,imageMeth,aic3sDel.data());
//			AIC3s_noTE(zetaPdel.data(),mSrc,nSrc,zetaPdel.data(),mTgt,nTgt,false,aic3sDel.data());
		} else {
			zetaPdel = zetaTgt + delZeta; // src lattice is wake
			AIC3s(zetaSrc_,mSrc,nSrc,zetaPdel.data(),mTgt,nTgt,imageMeth,aic3sDel.data());
//>>>>>>> rob/master
//			AIC3s_noTE(zetaSrc_,mSrc,nSrc,zetaPdel.data(),mTgt,nTgt,true,aic3sDel.data());
		}

		// diff downwash
		dU = aic3sDel*gamma0 - u;

		// add to matrix
		dX.block(0,qPri,S,1) = (1/del)*dU;
	}

	return;
}
