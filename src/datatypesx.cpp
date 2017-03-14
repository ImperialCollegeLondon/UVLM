/*
 * datatypesx.cpp
 *
 *  Created on: 30 Jan 2013
 *      Author: rjs10
 */

#include<datatypesx.hpp>
#include<stdio.h>
#include<vorticity.hpp>
#include<triads.hpp>

VortexSegment::VortexSegment(double* p1_, double* p2_,\
							 double* GammaMain_, double* GammaAdj_,\
							 bool atTE_, bool atImage_,\
							 unsigned int Paneli_, unsigned int Panelj_,\
							 unsigned int SegmentNo_)
							: p1(p1_), p2(p2_), \
							  GammaMain(GammaMain_), GammaAdj(GammaAdj_), \
							  atTE(atTE_), atImage(atImage_),
							  Paneli(Paneli_), Panelj(Panelj_),
							  SegmentNo(SegmentNo_) {
	//nothing else to do
}

double VortexSegment::Gamma() {
	if (GammaMain == NULL) {
		fprintf(stderr,"Un-init GammaMain in segment");
		exit(1);
	} else if (GammaAdj == NULL) {
		return *GammaMain;
	} else {
		return *GammaMain - *GammaAdj;
	} //END if, else
}

void VortexSegment::PanelEtas(double& eta1, double& eta2) {
	switch(SegmentNo) {
	case 1:
		eta1 = 0.0;
		eta2 = 0.5;
		break;
	case 2:
		eta1 = 0.5;
		eta2 = 1.0;
		break;
	case 3:
		eta1 = 1.0;
		eta2 = 0.5;
		break;
	case 4:
		eta1 = 0.5;
		eta2 = 0.0;
		break;
	default:
		fprintf(stderr,"PanelEtas: not a valid SegmentNo");
		exit(1);
	}
}

void VortexSegment::BiotSavart(double* pX, double* Uind, bool ImageMethod) {
	/**@brief Calculate self induced velocity at pX.
	 *
	 */

	//check if contribution is required
	if (this->Gamma() == 0.0 || atImage == true) {
		Uind[0] = 0.0;
		Uind[1] = 0.0;
		Uind[2] = 0.0;
		return;
	}

	//get object Gamma
	double Gamma = this->Gamma();

	if (ImageMethod == false) {
		C_BiotSegment(pX,p1,p2,Gamma,Uind);
	} else if (ImageMethod == true) {
		C_BiotSegment(pX,p1,p2,Gamma,Uind);

		//add image effect
		double TempTriad[3] = {0.0,0.0,0.0};
		C_BiotSegment_ImageYZ(pX,p1,p2,Gamma,TempTriad);
		AddTriad(Uind,TempTriad,Uind);
	}
}

void VortexSegment::Force(double* Uinc, double* F){
	/**@brief Calculate force on segment.
	 *
	 */

	//check if at TE or image plane
	if (atTE == true || atImage == true) {
		F[0] = 0.0;
		F[1] = 0.0;
		F[2] = 0.0;
		return;
	}

	double r0[3];

	//get orientation and length of segment
	SubTriad(p2,p1,r0);

	//get cross product
	CrossTriad(Uinc,r0,F);

	// multiply by gamma
	MulTriad(F,this->Gamma(),F);
}
