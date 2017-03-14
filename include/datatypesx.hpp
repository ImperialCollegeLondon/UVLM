/*
 * datatypesx.h
 *
 *  Created on: 30 Jan 2013
 *      Author: rjs10
 */

#ifndef DATATYPESX_H_
#define DATATYPESX_H_

#include<boost/multi_array.hpp>
#include<Eigen/Dense>

class VMopts {
public:
	unsigned int M;
	unsigned int N;
	bool ImageMethod;
	unsigned int Mstar;
	bool Steady;
	bool KJMeth;
	bool NewAIC;
	double DelTime;
	bool Rollup;
	unsigned int NumCores;
};

typedef boost::multi_array<double, 2> BoostArray2D;
typedef boost::multi_array<double, 3> BoostArray3D;
typedef Eigen::Map<Eigen::VectorXd> EigenMapVectXd;
typedef Eigen::Map<const Eigen::VectorXd> ConstMapVectXd;
typedef Eigen::Map<Eigen::Vector3d> MapVect3d;
typedef Eigen::Map<const Eigen::Vector3d> ConstMapVect3d;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>  EigDynMatrixRM;
typedef Eigen::Map<EigDynMatrixRM> EigenMapMatrixXd;

class VortexSegment {
	/**@brief Segment class for use with KJ method and fast rollup.
	 * @details SegmentNo refers to the position of segment within panel
	 *(Paneli,Panelj) ordered 1-2-3-4 which is i,j -> i,j+1 -> i+1,j+1 -> i+1,j.
	 */
public:
	VortexSegment(double* p1_, double* p2_,\
				  double* GammaMain_, double* GammaAdj_,\
				  bool atTE_, bool atImage_,\
				  unsigned int Paneli_, unsigned int Panelj_,\
				  unsigned int SegmentNo_);
	double* p1;
	double* p2;
	double* GammaMain;
	double* GammaAdj;
	bool atTE;
	bool atImage;
	unsigned int Paneli, Panelj;
	unsigned int SegmentNo;

	double Gamma(void);
	void PanelEtas(double &eta1, double &eta2);
	void BiotSavart(double* pX, double* Uind, bool ImageMethod);
	void Force(double* Uinc, double* F);
};

#endif /* DATATYPESX_H_ */
