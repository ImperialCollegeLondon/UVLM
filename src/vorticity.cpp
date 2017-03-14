/**
 * vorticity.cpp
 *
 *  Created on: 28 Jan 2013
 *      Author: rjs10
 */
#include <vorticity.hpp>
#include <triads.hpp>
#include <cmath>
#include <stdio.h>

Eigen::Vector3d BiotSegment(const Eigen::Vector3d& xp, \
							const Eigen::Vector3d& x1, \
							const Eigen::Vector3d& x2, \
							const double& gam) {
	/** @brief Calculate the velocity induced at point p by a vortex line
	  * with start and end points 1 and 2.
	  * @brief Velocity induced by segment x1->x2 at point xP.
      *	@param xp Point at which we require velocity.
      *	@param x1 Vortex segment start point.
      *	@param x2 Vortex segment end point.
      *	@param gam Circulation strength of vortex segment.
      *	@returns Eigen::Vector3d Velocity induced at point xP.
	  * @details Solved using Eigen 3rd party library.
	  * @details If r is less than the vortex core radius then u = v = w = 0.
	  */


    if (gam == 0) {return Eigen::Vector3d(0,0,0);}


    //define position vectors
    Eigen::Vector3d r0 = (x2-x1);
    Eigen::Vector3d r1 = (xp-x1);
    Eigen::Vector3d r2 = (xp-x2);
    Eigen::Vector3d rx = r1.cross(r2);


    // if the point P is within the core radius of p1 p2 or p1Xp2 return zeros
    if((r1.norm() <= 1.0e-5) || (r2.norm() <= 1.0e-5) || (rx.norm() <= 1.0e-5)){
        return Eigen::Vector3d(0,0,0);
    }


    //vector dot products
    double r0d1 = r0.dot(r1);
    double r0d2 = r0.dot(r2);

    //K defined as in Katz & Plotkin
    double K = (gam/(4*M_PI*pow(rx.norm(),2))) *
                    ((r0d1/r1.norm())-(r0d2/r2.norm()));

    return Eigen::Vector3d( K*rx(0), K*rx(1), K*rx(2) );
}

void BiotSegment_map(const Eigen::Map<Eigen::Vector3d>& xp, \
								const Eigen::Map<Eigen::Vector3d>& x1, \
								const Eigen::Map<Eigen::Vector3d>& x2, \
								const double& gam,\
								const Eigen::Map<Eigen::Vector3d>& Uind) {
    /** @brief * @brief Calculate the velocity induced at point p by a vortex line
	* with start and end points 1 and 2.
	* @brief Velocity induced by segment x1->x2 at point xP.
    * @param xp Point at which we require velocity.
    * @param x1 Vortex segment start point.
    * @param x2 Vortex segment end point.
    * @param gam Circulation strength of vortex segment.
    * @param Uind Induced velocity (Eigen pointer to existing memory).
    */

    if (gam == 0.0) {
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[0] = 0.0;
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[1] = 0.0;
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[2] = 0.0;
    	return;
    }


    //define position vectors
    Eigen::Vector3d r0 = (x2-x1);
    Eigen::Vector3d r1 = (xp-x1);
    Eigen::Vector3d r2 = (xp-x2);
    Eigen::Vector3d rx = r1.cross(r2);


    // if the point P is within the core radius of p1 p2 or p1Xp2 return zeros
    if((r1.norm() <= 1.0e-5) || (r2.norm() <= 1.0e-5) || (rx.norm() <= 1.0e-5)){
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[0] = 0.0;
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[1] = 0.0;
    	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[2] = 0.0;
    	return;
    }


    //vector dot products
    double r0d1 = r0.dot(r1);
    double r0d2 = r0.dot(r2);

    //K defined as in Katz & Plotkin
    double K = (gam/(4*M_PI*pow(rx.norm(),2))) *
                    ((r0d1/r1.norm())-(r0d2/r2.norm()));


	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[0] = K*rx[0];
	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[1] = K*rx[1];
	const_cast<Eigen::Map<Eigen::Vector3d>&>(Uind)[2] = K*rx[2];
    return;
}

void C_BiotSegment(const double* xP,\
					  const double* x1,\
					  const double* x2,\
					  const double& Gamma,\
					  double* Uind) {
	/** @brief Calculate the velocity induced at point p by a vortex line
	  * with start and end points 1 and 2.
	  * @brief Velocity induced by segment x1->x2 at point xP.
      *	@param xP Point at which we require velocity.
      *	@param x1 Vortex segment start point.
      *	@param x2 Vortex segment end point.
      *	@param Gamma Circulation strength of vortex segment.
      *	@param Uind Velocity induced at point xP.
	  * @details Solved using pointers and dereferencing only.
	  * @details If r is less than the vortex core radius then u = v = w = 0.
	  */

	if (Gamma == 0.0) {
		Uind[0] = 0.0;
		Uind[1] = 0.0;
		Uind[2] = 0.0;
		return;
	}

	double r0[] = {0.0,0.0,0.0};
	double r1[] = {0.0,0.0,0.0};
	double r2[] = {0.0,0.0,0.0};
	double rx[] = {0.0,0.0,0.0};

	SubTriad(x2,x1,r0);
	SubTriad(xP,x1,r1);
	SubTriad(xP,x2,r2);
	CrossTriad(r1,r2,rx);

	if((NormTriad(r1) <= 1.0e-5) || (NormTriad(r2) <= 1.0e-5) || \
			(NormTriad(rx) <= 1.0e-5)) {
		Uind[0] = 0.0;
		Uind[1] = 0.0;
		Uind[2] = 0.0;
	    return;
	}


	double r0d1 = DotTriad(r0,r1);
	double r0d2 = DotTriad(r0,r2);

	double K = (Gamma/(4*M_PI*pow(NormTriad(rx),2))) *
	                    ((r0d1/NormTriad(r1)-(r0d2/NormTriad(r2))));

	//overwrite elements of Uind
	Uind[0] = K*rx[0];
	Uind[1] = K*rx[1];
	Uind[2] = K*rx[2];
}

void C_BiotSegment_ImageYZ(const double* xP,\
		  const double* x1,\
		  const double* x2,\
		  double& Gamma,\
		  double* Uind) {
	/** @brief calculate effect due to mirror of segment in YZ-plane.
	 *  @details Uind is overwritten.
	 */
	// create target image location (reverse x ordinate)
	double xImage[3] = {0.0,0.0,0.0};
	xImage[0] = -xP[0];
	xImage[1] = xP[1];
	xImage[2] = xP[2];

	//create result
	double VelImage[3] = {0.0,0.0,0.0};

	// segment 2 (point 2 -> point 3)
	C_BiotSegment(xImage, \
				  x1,x2, \
				  Gamma, VelImage);

	//reverse x velocity
	VelImage[0] = -VelImage[0];

	//copy to output
	CopyTriad(Uind,VelImage);
}


void BiotSavartSurf(const double* Zeta_Vec, const double* Gamma_Vec, \
					const double* TargetTriad,\
					const unsigned int Mstart, const unsigned int Nstart,
					const unsigned int Mend, const unsigned int Nend,
					const unsigned int Mfull, const unsigned int Nfull, \
					const bool ImageMethod, \
					double* Uout) {
	/**@brief Velocity induced by grid.
	 * @param Zeta_Vec Aero grid - size (M+1)*(N+1)*3.
	 * @param Gamma_Vec Vortex ring circulation strengths - size M*N.
	 * @param TargetTriad Point at which velocity is required - size 3.
	 */

	//create some useful pointers
	double (*Zeta)[Nfull+1][3] = (double (*)[Nfull+1][3]) Zeta_Vec;
	double (*Gamma)[Nfull] = (double (*)[Nfull]) Gamma_Vec;

	//slow method (every segment of every ring)
	//TODO: differencing for individual segment strengths, required for KJ methods.

	//set Uout to zero
	Uout[0] = 0.0;
	Uout[1] = 0.0;
	Uout[2] = 0.0;

	//define temp variables
	double Temp1[3] = {0.0,0.0,0.0};
	double Temp2[3] = {0.0,0.0,0.0};

	//TODO: to make this parallel use locally defined variables NOT pointers
	//to global data
	for (unsigned int i = Mstart; i < Mend; i++) {
		for (unsigned int j = Nstart; j < Nend; j++) {
			// loop through each panel
			// reset running total
			Temp2[0] = 0.0;
			Temp2[1] = 0.0;
			Temp2[2] = 0.0;

			// get effect of all segments
			// segment 1 (point 1 -> point 2)
			C_BiotSegment(TargetTriad,Zeta[i][j],Zeta[i][j+1],
						  Gamma[i][j], Temp1);
			AddTriad(Temp2,Temp1,Temp2);

			// segment 2 (point 2 -> point 3)
			C_BiotSegment(TargetTriad,Zeta[i][j+1],Zeta[i+1][j+1],
						  Gamma[i][j], Temp1);
			AddTriad(Temp2,Temp1,Temp2);

			// segment 3 (point 3 -> point 4)
			C_BiotSegment(TargetTriad,Zeta[i+1][j+1],Zeta[i+1][j],
						  Gamma[i][j], Temp1);
			AddTriad(Temp2,Temp1,Temp2);

			// segment 3 (point 4 -> point 1)
			C_BiotSegment(TargetTriad,Zeta[i+1][j],Zeta[i][j],
						  Gamma[i][j], Temp1);
			AddTriad(Temp2,Temp1,Temp2);

			if (ImageMethod == 1) {
				// get effect of all segments
				// segment 1 (point 1 -> point 2) image
				C_BiotSegment_ImageYZ(TargetTriad,Zeta[i][j],Zeta[i][j+1],
							  Gamma[i][j], Temp1);
				AddTriad(Temp2,Temp1,Temp2);

				// segment 2 (point 2 -> point 3) image
				C_BiotSegment_ImageYZ(TargetTriad,Zeta[i][j+1],Zeta[i+1][j+1],
							  Gamma[i][j], Temp1);
				AddTriad(Temp2,Temp1,Temp2);

				// segment 3 (point 3 -> point 4) image
				C_BiotSegment_ImageYZ(TargetTriad,Zeta[i+1][j+1],Zeta[i+1][j],
							  Gamma[i][j], Temp1);
				AddTriad(Temp2,Temp1,Temp2);

				// segment 3 (point 4 -> point 1) image
				C_BiotSegment_ImageYZ(TargetTriad,Zeta[i+1][j],Zeta[i][j],
							  Gamma[i][j], Temp1);
				AddTriad(Temp2,Temp1,Temp2);
			}

			//Add to Uout
			AddTriad(Uout,Temp2,Uout);
		} // END for j
	} //END for i
}



void BiotSurfaceSegments(std::vector<VortexSegment>& Seg,\
						 unsigned int SurfM,\
						 unsigned int SurfN,\
						 bool ImageMethod,\
						 double* Target,\
						 double* Uind) {
	/**@brief Velocity induced by segments initialised on M*N grid.
	 * @param Target Point at which velocity is required - size 3.
	 * @param Uind Induced velocity.
	 */
	for (unsigned int k = 0; k < (SurfM+1)*SurfN + SurfM*(SurfN+1); k++) {
		double Temp[3];
		Seg[k].BiotSavart(Target,Temp,ImageMethod);
		AddTriad(Uind,Temp,Uind);
	}
}
