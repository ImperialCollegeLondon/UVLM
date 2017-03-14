/**@brief      quasi-steady VLM solution with fixed-wake.
 * @author     Rob Simpson
 * @contact    r.simpson11@imperial.ac.uk
 * @version    0.0
 * @date       30/01/2013
 * @pre        None
 * @warning    None
 */

#include <stdio.h>
#include <fstream>
#include <datatypesx.hpp>
#include <triads.hpp>
#include <vorticity.hpp>
#include <PanelTools.hpp>
#include <vector>
#include <assert.h>
#include "omp.h"
#include <aicMats.hpp>
#include <iostream>


void KJMethodForces(const double* Zeta_Vec, const double* Gamma_Vec,\
		const double* ZetaStar_Vec, const double* GammaStar_Vec,\
		const double* ZetaDot_Vec, \
		const double* Uext_Vec, \
		VMopts VMOPTS,\
		const double* Gamma_tm1_Vec,\
		double* Forces_Vec) {
	/**@brief calculate unsteady panel forces then segment forces
	 * @param Zeta_Vec Grid points.
	 * @param Gamma_Vec Bound vortex ring gammas.
	 * @param ZetaStar_Vec Wake grid points.
	 * @param GammaStar_Vec Wake vortex ring gammas.
	 * @param VMOPTS Simulation options.
	 * @param Gamma_tm1_Vec Bound vortex ring gammas from previous timestep.
	 * @param Forces_Vec Unsteady Forces at grid points.
	 */

	// cast vectors into useful pointers
	double (*Zeta)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Zeta_Vec;
	double (*ZetaDot)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) ZetaDot_Vec;
	double (*Uext)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Uext_Vec;
	double (*Gamma)[VMOPTS.N] = (double (*)[VMOPTS.N]) Gamma_Vec;
	double (*Gamma_tm1)[VMOPTS.N] = (double (*)[VMOPTS.N]) Gamma_tm1_Vec;
	double (*Forces)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Forces_Vec;

	// declare local, automatically-managed dynamic memory using boost ...
	BoostArray2D Area_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray2D dGamma_dt_(boost::extents[VMOPTS.M][VMOPTS.N]);

	// ... and get useful pointers to that data
	double (*Area)[VMOPTS.N] = (double (*)[VMOPTS.N]) Area_.data();
	double (*dGamma_dt)[VMOPTS.N] = (double (*)[VMOPTS.N]) dGamma_dt_.data();

	//Temporary variables
	double Normal[3] = {0.0,0.0,0.0};
	double DeltaP = 0.0;
	double NormalForce = 0.0;
	double UnstForceTemp[3] = {0.0,0.0,0.0};
	double* NullDouble = NULL;

	//get unsteady component of pressure forces and map to nodes
	//loop through collocation points i,j
	for (unsigned int i = 0; i < VMOPTS.M; i++) {
		for (unsigned int j = 0; j < VMOPTS.N; j++) {

			//calculate panel normals
			PanelNormal(Zeta[i][j], Zeta[i][j+1], \
						Zeta[i+1][j+1], Zeta[i+1][j],
									Normal);


			// calculate panel areas
			Area[i][j] = PanelArea(Zeta[i][j], Zeta[i][j+1],\
								   Zeta[i+1][j+1], Zeta[i+1][j]);


			// calculate panel d(Gamma)/dt
			if (VMOPTS.Steady == true) {
				dGamma_dt[i][j] = 0.0;
			} else if (VMOPTS.Steady == false) {
				dGamma_dt[i][j] = (Gamma[i][j] - Gamma_tm1[i][j]) / \
									VMOPTS.DelTime;
			}

			// pressure jump
			DeltaP = dGamma_dt[i][j];

			//calculate corresponding normal force
			NormalForce = DeltaP*Area[i][j];


			// apply along panel normal
			MulTriad(Normal,NormalForce,UnstForceTemp);

			// map to corner points
			BilinearInterpTriad(Forces[i][j],Forces[i][j+1],\
								Forces[i+1][j+1],Forces[i+1][j],\
								UnstForceTemp,\
								0.5,0.5,true);

		} //END for j
	}//END for i

	//get vector of segments
	std::vector<VortexSegment> Seg;

	//calculate expected number of segements
	unsigned int K = (VMOPTS.M+1)*VMOPTS.N + VMOPTS.M*(VMOPTS.N+1);

	//initialise segment list
	for (unsigned int i = 0; i < VMOPTS.M; i++) {
		for (unsigned int j = 0; j < VMOPTS.N; j++) {

			//segment 4 of panels (:,:)
			if (VMOPTS.ImageMethod == 1 && j==0) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j][0],&Zeta[i][j][0],\
								  	  	  	&Gamma[i][j], &Gamma[i][j],\
								  	  	  	false,true,i,j,4));

			} else if (VMOPTS.ImageMethod == 0 && j==0) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j][0],&Zeta[i][j][0],\
											&Gamma[i][j], NullDouble,\
											false,false,i,j,4));

			} else if (j>0) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j][0],&Zeta[i][j][0],\
											&Gamma[i][j], &Gamma[i][j-1],\
											false,false,i,j,4));

			} //end if, else-if, else-if

			//segment 1 of panels (:,:)
			if (i == 0) {
				Seg.push_back(VortexSegment(&Zeta[i][j][0],&Zeta[i][j+1][0],\
											&Gamma[i][j], NullDouble,\
											false,false,i,j,1));

			} else if (i > 0) {
				Seg.push_back(VortexSegment(&Zeta[i][j][0],&Zeta[i][j+1][0],\
											&Gamma[i][j], &Gamma[i-1][j],\
											false,false,i,j,1));

			} // end if, else-if

			//segment 2 of panels(:,N-1)
			if (j == VMOPTS.N - 1) {
				Seg.push_back(VortexSegment(&Zeta[i][j+1][0],&Zeta[i+1][j+1][0],\
											&Gamma[i][j], NullDouble,\
											false,false,i,j,2));

			} //end if

			//add segment 3 if at the TE
			if (i == VMOPTS.M - 1) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j+1][0],&Zeta[i+1][j][0],\
											&Gamma[i][j], NullDouble,\
											true,false,i,j,3));
			} //end if
		}
	}

	//check length
	assert(Seg.size() == K);

	//print segment info
//	for (unsigned int k = 0; k < (VMOPTS.M+1)*VMOPTS.N + VMOPTS.M*(VMOPTS.N+1); k++) {
//		printf("Seg %d : i=%d ,j=%d, No=%d, Gamma=%f p1,p2 = ",\
//				k, Seg[k].Paneli,Seg[k].Panelj,Seg[k].SegmentNo,Seg[k].Gamma());
//		PrintTriad(Seg[k].p1);
//		PrintTriad(Seg[k].p2);
//		printf("\n");
//	}


	//test segments v surface
//	double Target[3] = {0.5,-0.5,1.0};
//	double Surf[3] = {0.0,0.0,0.0};
//	double Segs[3] = {0.0,0.0,0.0};
//
//	BiotSavartSurf(Zeta_Vec, Gamma_Vec, Target,\
//				   0, 0, \
//				   VMOPTS.M, VMOPTS.N,\
//				   VMOPTS.M, VMOPTS.N,\
//				   VMOPTS.ImageMethod,\
//				   Surf);
//
//	PrintTriad(Surf);
//
//	BiotSurfaceSegments(Seg, VMOPTS.M, VMOPTS.N, VMOPTS.ImageMethod,\
//						Target,Segs);
//
//	PrintTriad(Segs);


	//loop through each segment and add force to Forces
	#pragma omp parallel for
	for (unsigned int k = 0; k < K; k++) {
		//get panel i,j
		unsigned int i = Seg[k].Paneli;
		unsigned int j = Seg[k].Panelj;

		//declare local vars
		double ZetaDotSeg[3] = {0.0,0.0,0.0};
		double UextSeg[3] = {0.0,0.0,0.0};
		double UwakeSeg[3] = {0.0,0.0,0.0};
		double UsurfSeg[3] = {0.0,0.0,0.0};
		double UincSeg[3] = {0.0,0.0,0.0};
		double ForceSeg[3] = {0.0,0.0,0.0};

		//get midpoint location in panel dimensionless coords
		double Eta1 = 0.0;
		double Eta2 = 0.0;
		Seg[k].PanelEtas(Eta1,Eta2);

		//calc incident veloc
		// -ZetaDot + Uext + Uwake + Ubound

		//ZetaDot at midpoint
		BilinearInterpTriad(ZetaDot[i][j], ZetaDot[i][j+1], \
						 	ZetaDot[i+1][j+1], ZetaDot[i+1][j], \
						 	ZetaDotSeg,Eta1,Eta2,false);

		//Uext at midpoint
		BilinearInterpTriad(Uext[i][j], Uext[i][j+1], \
							Uext[i+1][j+1], Uext[i+1][j], \
							UextSeg,Eta1,Eta2,false);

		//Uwake at midpoint using segment list
		//get midpoint
		double Midpoint[3] = {0.0,0.0,0.0};
		BilinearInterpTriad(Zeta[i][j], Zeta[i][j+1], \
							Zeta[i+1][j+1], Zeta[i+1][j], \
							Midpoint,Eta1,Eta2,false);

		//Uwake
		BiotSavartSurf(ZetaStar_Vec,GammaStar_Vec,Midpoint,\
					   0,0,\
					   VMOPTS.Mstar,VMOPTS.N,\
					   VMOPTS.Mstar,VMOPTS.N,\
					   VMOPTS.ImageMethod, UwakeSeg);

		//Usurf
		BiotSurfaceSegments(Seg,VMOPTS.M,VMOPTS.N,\
				   	   	   	VMOPTS.ImageMethod,Midpoint,\
				   	   	   	UsurfSeg);

		//Combine Uext - ZetaDot
		SubTriad(UextSeg,ZetaDotSeg,UincSeg);
		//Add on Uwake
		AddTriad(UincSeg,UwakeSeg,UincSeg);
		//Add on Usurf
		AddTriad(UincSeg,UsurfSeg,UincSeg);

//		std::cout << "downwash at seg " << Seg[k].SegmentNo << std::endl;
//		PrintTriad(UincSeg);

		//Get force
		Seg[k].Force(UincSeg,ForceSeg);


//		printf("Seg %d : i=%d ,j=%d, No=%d, Gamma=%f Uinc,Force,Uwake = ",\
//				k, Seg[k].Paneli,Seg[k].Panelj,Seg[k].SegmentNo,Seg[k].Gamma());
//		PrintTriad(UincSeg);
//		PrintTriad(ForceSeg);
//		PrintTriad(UwakeSeg);
////		printf("eta1=%f\teta2=%f",Eta1,Eta2);
//		printf("\n");

		//Map to corner points
		#pragma omp critical
		BilinearInterpTriad(Forces[i][j], Forces[i][j+1], \
							Forces[i+1][j+1], Forces[i+1][j], \
							ForceSeg,Eta1,Eta2,true);
	}

}

void KJMethodForces_vC(const double* Zeta_Vec, const double* Gamma_Vec,\
		const double* ZetaStar_Vec, const double* GammaStar_Vec,\
		const double* ZetaDot_Vec, \
		const double* Uext_Vec, \
		VMopts VMOPTS,\
		const double* Gamma_tm1_Vec,\
		double* Forces_Vec) {
	/**@brief calculate unsteady panel forces then segment forces with collocation velocity approximation.
	 * @param Zeta_Vec Grid points.
	 * @param Gamma_Vec Bound vortex ring gammas.
	 * @param ZetaStar_Vec Wake grid points.
	 * @param GammaStar_Vec Wake vortex ring gammas.
	 * @param VMOPTS Simulation options.
	 * @param Gamma_tm1_Vec Bound vortex ring gammas from previous timestep.
	 * @param Forces_Vec Unsteady Forces at grid points.
	 */

	// cast vectors into useful pointers
	double (*Zeta)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Zeta_Vec;
	double (*ZetaDot)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) ZetaDot_Vec;
	double (*Uext)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Uext_Vec;
	double (*Gamma)[VMOPTS.N] = (double (*)[VMOPTS.N]) Gamma_Vec;
	double (*Gamma_tm1)[VMOPTS.N] = (double (*)[VMOPTS.N]) Gamma_tm1_Vec;
	double (*Forces)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Forces_Vec;

	// declare local, automatically-managed dynamic memory using boost ...
	BoostArray2D Area_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray2D dGamma_dt_(boost::extents[VMOPTS.M][VMOPTS.N]);

	// ... and get useful pointers to that data
	double (*Area)[VMOPTS.N] = (double (*)[VMOPTS.N]) Area_.data();
	double (*dGamma_dt)[VMOPTS.N] = (double (*)[VMOPTS.N]) dGamma_dt_.data();

	//Temporary variables
	double Normal[3] = {0.0,0.0,0.0};
	double DeltaP = 0.0;
	double NormalForce = 0.0;
	double UnstForceTemp[3] = {0.0,0.0,0.0};
	double* NullDouble = NULL;

	//get unsteady component of pressure forces and map to nodes
	//loop through collocation points i,j
	for (unsigned int i = 0; i < VMOPTS.M; i++) {
		for (unsigned int j = 0; j < VMOPTS.N; j++) {

			//calculate panel normals
			PanelNormal(Zeta[i][j], Zeta[i][j+1], \
						Zeta[i+1][j+1], Zeta[i+1][j],
									Normal);


			// calculate panel areas
			Area[i][j] = PanelArea(Zeta[i][j], Zeta[i][j+1],\
								   Zeta[i+1][j+1], Zeta[i+1][j]);


			// calculate panel d(Gamma)/dt
			if (VMOPTS.Steady == true) {
				dGamma_dt[i][j] = 0.0;
			} else if (VMOPTS.Steady == false) {
				dGamma_dt[i][j] = (Gamma[i][j] - Gamma_tm1[i][j]) / \
									VMOPTS.DelTime;
			}

			// pressure jump
			DeltaP = dGamma_dt[i][j];

			//calculate corresponding normal force
			NormalForce = DeltaP*Area[i][j];


			// apply along panel normal
			MulTriad(Normal,NormalForce,UnstForceTemp);

			// map to corner points
			BilinearInterpTriad(Forces[i][j],Forces[i][j+1],\
								Forces[i+1][j+1],Forces[i+1][j],\
								UnstForceTemp,\
								0.5,0.5,true);

		} //END for j
	}//END for i

	//get vector of segments
	std::vector<VortexSegment> Seg;

	//calculate expected number of segements
	unsigned int K = (VMOPTS.M+1)*VMOPTS.N + VMOPTS.M*(VMOPTS.N+1);

	//initialise segment list
	for (unsigned int i = 0; i < VMOPTS.M; i++) {
		for (unsigned int j = 0; j < VMOPTS.N; j++) {

			//segment 4 of panels (:,:)
			if (VMOPTS.ImageMethod == 1 && j==0) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j][0],&Zeta[i][j][0],\
								  	  	  	&Gamma[i][j], &Gamma[i][j],\
								  	  	  	false,true,i,j,4));

			} else if (VMOPTS.ImageMethod == 0 && j==0) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j][0],&Zeta[i][j][0],\
											&Gamma[i][j], NullDouble,\
											false,false,i,j,4));

			} else if (j>0) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j][0],&Zeta[i][j][0],\
											&Gamma[i][j], &Gamma[i][j-1],\
											false,false,i,j,4));

			} //end if, else-if, else-if

			//segment 1 of panels (:,:)
			if (i == 0) {
				Seg.push_back(VortexSegment(&Zeta[i][j][0],&Zeta[i][j+1][0],\
											&Gamma[i][j], NullDouble,\
											false,false,i,j,1));

			} else if (i > 0) {
				Seg.push_back(VortexSegment(&Zeta[i][j][0],&Zeta[i][j+1][0],\
											&Gamma[i][j], &Gamma[i-1][j],\
											false,false,i,j,1));

			} // end if, else-if

			//segment 2 of panels(:,N-1)
			if (j == VMOPTS.N - 1) {
				Seg.push_back(VortexSegment(&Zeta[i][j+1][0],&Zeta[i+1][j+1][0],\
											&Gamma[i][j], NullDouble,\
											false,false,i,j,2));

			} //end if

			//add segment 3 if at the TE
			if (i == VMOPTS.M - 1) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j+1][0],&Zeta[i+1][j][0],\
											&Gamma[i][j], NullDouble,\
											true,false,i,j,3));
			} //end if
		}
	}

	//check length
	assert(Seg.size() == K);

	//loop through each segment and add force to Forces
	#pragma omp parallel for
	for (unsigned int k = 0; k < K; k++) {
		//get panel i,j
		unsigned int i = Seg[k].Paneli;
		unsigned int j = Seg[k].Panelj;

		//declare local vars
		double ZetaDotSeg[3] = {0.0,0.0,0.0};
		double UextSeg[3] = {0.0,0.0,0.0};
		double UwakeSeg[3] = {0.0,0.0,0.0};
		double UsurfSeg[3] = {0.0,0.0,0.0};
		double UincSeg[3] = {0.0,0.0,0.0};
		double ForceSeg[3] = {0.0,0.0,0.0};

		//set collocation location in panel dimensionless coords
		double Eta1 = 0.5;
		double Eta2 = 0.5;

		//calc incident veloc
		// -ZetaDot + Uext + Uwake + Ubound

		//ZetaDot at collocation
		BilinearInterpTriad(ZetaDot[i][j], ZetaDot[i][j+1], \
						 	ZetaDot[i+1][j+1], ZetaDot[i+1][j], \
						 	ZetaDotSeg,Eta1,Eta2,false);

		//Uext at collocation
		BilinearInterpTriad(Uext[i][j], Uext[i][j+1], \
							Uext[i+1][j+1], Uext[i+1][j], \
							UextSeg,Eta1,Eta2,false);

		//Uwake at collocation using segment list
		//get collocation
		double Col[3] = {0.0,0.0,0.0};
		BilinearInterpTriad(Zeta[i][j], Zeta[i][j+1], \
							Zeta[i+1][j+1], Zeta[i+1][j], \
							Col,Eta1,Eta2,false);

		//Uwake
		BiotSavartSurf(ZetaStar_Vec,GammaStar_Vec,Col,\
					   0,0,\
					   VMOPTS.Mstar,VMOPTS.N,\
					   VMOPTS.Mstar,VMOPTS.N,\
					   VMOPTS.ImageMethod, UwakeSeg);

		//Usurf
		BiotSurfaceSegments(Seg,VMOPTS.M,VMOPTS.N,\
				   	   	   	VMOPTS.ImageMethod,Col,\
				   	   	   	UsurfSeg);

		//Combine Uext - ZetaDot
		SubTriad(UextSeg,ZetaDotSeg,UincSeg);
		//Add on Uwake
		AddTriad(UincSeg,UwakeSeg,UincSeg);
		//Add on Usurf
		AddTriad(UincSeg,UsurfSeg,UincSeg);

//		std::cout << "downwash at seg " << Seg[k].SegmentNo << std::endl;
//		PrintTriad(UincSeg);

		//Get force
		Seg[k].Force(UincSeg,ForceSeg);

		//Map to corner points
		#pragma omp critical
		BilinearInterpTriad(Forces[i][j], Forces[i][j+1], \
							Forces[i+1][j+1], Forces[i+1][j], \
							ForceSeg,Eta1,Eta2,true);
	}

}

void KJMethodForces_vC_mod(const double* Zeta_Vec, const double* Gamma_Vec,\
		const double* ZetaStar_Vec, const double* GammaStar_Vec,\
		const double* ZetaDot_Vec, \
		const double* Uext_Vec, \
		VMopts VMOPTS,\
		const double* Gamma_tm1_Vec,\
		double* Forces_Vec) {
	/**@brief calculate unsteady panel forces then segment forces with
	 * collocation velocity approximation and modification omitting vortex
	 * segment self influence.
	 * @param Zeta_Vec Grid points.
	 * @param Gamma_Vec Bound vortex ring gammas.
	 * @param ZetaStar_Vec Wake grid points.
	 * @param GammaStar_Vec Wake vortex ring gammas.
	 * @param VMOPTS Simulation options.
	 * @param Gamma_tm1_Vec Bound vortex ring gammas from previous timestep.
	 * @param Forces_Vec Unsteady Forces at grid points.
	 */

	// cast vectors into useful pointers
	double (*Zeta)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Zeta_Vec;
	double (*ZetaDot)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) ZetaDot_Vec;
	double (*Uext)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Uext_Vec;
	double (*Gamma)[VMOPTS.N] = (double (*)[VMOPTS.N]) Gamma_Vec;
	double (*Gamma_tm1)[VMOPTS.N] = (double (*)[VMOPTS.N]) Gamma_tm1_Vec;
	double (*Forces)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Forces_Vec;

	// declare local, automatically-managed dynamic memory using boost ...
	BoostArray2D Area_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray2D dGamma_dt_(boost::extents[VMOPTS.M][VMOPTS.N]);

	// ... and get useful pointers to that data
	double (*Area)[VMOPTS.N] = (double (*)[VMOPTS.N]) Area_.data();
	double (*dGamma_dt)[VMOPTS.N] = (double (*)[VMOPTS.N]) dGamma_dt_.data();

	//Temporary variables
	double Normal[3] = {0.0,0.0,0.0};
	double DeltaP = 0.0;
	double NormalForce = 0.0;
	double UnstForceTemp[3] = {0.0,0.0,0.0};
	double* NullDouble = NULL;

	//get unsteady component of pressure forces and map to nodes
	//loop through collocation points i,j
	for (unsigned int i = 0; i < VMOPTS.M; i++) {
		for (unsigned int j = 0; j < VMOPTS.N; j++) {

			//calculate panel normals
			PanelNormal(Zeta[i][j], Zeta[i][j+1], \
						Zeta[i+1][j+1], Zeta[i+1][j],
									Normal);


			// calculate panel areas
			Area[i][j] = PanelArea(Zeta[i][j], Zeta[i][j+1],\
								   Zeta[i+1][j+1], Zeta[i+1][j]);


			// calculate panel d(Gamma)/dt
			if (VMOPTS.Steady == true) {
				dGamma_dt[i][j] = 0.0;
			} else if (VMOPTS.Steady == false) {
				dGamma_dt[i][j] = (Gamma[i][j] - Gamma_tm1[i][j]) / \
									VMOPTS.DelTime;
			}

			// pressure jump
			DeltaP = dGamma_dt[i][j];

			//calculate corresponding normal force
			NormalForce = DeltaP*Area[i][j];


			// apply along panel normal
			MulTriad(Normal,NormalForce,UnstForceTemp);

			// map to corner points
			BilinearInterpTriad(Forces[i][j],Forces[i][j+1],\
								Forces[i+1][j+1],Forces[i+1][j],\
								UnstForceTemp,\
								0.5,0.5,true);

		} //END for j
	}//END for i

	//get vector of segments
	std::vector<VortexSegment> Seg;

	//calculate expected number of segements
	unsigned int K = (VMOPTS.M+1)*VMOPTS.N + VMOPTS.M*(VMOPTS.N+1);

	//initialise segment list
	for (unsigned int i = 0; i < VMOPTS.M; i++) {
		for (unsigned int j = 0; j < VMOPTS.N; j++) {

			//segment 4 of panels (:,:)
			if (VMOPTS.ImageMethod == 1 && j==0) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j][0],&Zeta[i][j][0],\
								  	  	  	&Gamma[i][j], &Gamma[i][j],\
								  	  	  	false,true,i,j,4));

			} else if (VMOPTS.ImageMethod == 0 && j==0) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j][0],&Zeta[i][j][0],\
											&Gamma[i][j], NullDouble,\
											false,false,i,j,4));

			} else if (j>0) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j][0],&Zeta[i][j][0],\
											&Gamma[i][j], &Gamma[i][j-1],\
											false,false,i,j,4));

			} //end if, else-if, else-if

			//segment 1 of panels (:,:)
			if (i == 0) {
				Seg.push_back(VortexSegment(&Zeta[i][j][0],&Zeta[i][j+1][0],\
											&Gamma[i][j], NullDouble,\
											false,false,i,j,1));

			} else if (i > 0) {
				Seg.push_back(VortexSegment(&Zeta[i][j][0],&Zeta[i][j+1][0],\
											&Gamma[i][j], &Gamma[i-1][j],\
											false,false,i,j,1));

			} // end if, else-if

			//segment 2 of panels(:,N-1)
			if (j == VMOPTS.N - 1) {
				Seg.push_back(VortexSegment(&Zeta[i][j+1][0],&Zeta[i+1][j+1][0],\
											&Gamma[i][j], NullDouble,\
											false,false,i,j,2));

			} //end if

			//add segment 3 if at the TE
			if (i == VMOPTS.M - 1) {
				Seg.push_back(VortexSegment(&Zeta[i+1][j+1][0],&Zeta[i+1][j][0],\
											&Gamma[i][j], NullDouble,\
											true,false,i,j,3));
			} //end if
		}
	}

	//check length
	assert(Seg.size() == K);

	//loop through each segment and add force to Forces
	for (unsigned int k = 0; k < K; k++) {
		//get panel i,j
		unsigned int i = Seg[k].Paneli;
		unsigned int j = Seg[k].Panelj;

		//declare local vars
		double ZetaDotSeg[3] = {0.0,0.0,0.0};
		double UextSeg[3] = {0.0,0.0,0.0};
		double UwakeSeg[3] = {0.0,0.0,0.0};
		double UsurfSeg[3] = {0.0,0.0,0.0};
		double UincSeg[3] = {0.0,0.0,0.0};
		double ForceSeg[3] = {0.0,0.0,0.0};

		//set collocation location in panel dimensionless coords
		double Eta1 = 0.5;
		double Eta2 = 0.5;

		//calc incident veloc
		// -ZetaDot + Uext + Uwake + Ubound

		//ZetaDot at collocation
		BilinearInterpTriad(ZetaDot[i][j], ZetaDot[i][j+1], \
						 	ZetaDot[i+1][j+1], ZetaDot[i+1][j], \
						 	ZetaDotSeg,Eta1,Eta2,false);

		//Uext at collocation
		BilinearInterpTriad(Uext[i][j], Uext[i][j+1], \
							Uext[i+1][j+1], Uext[i+1][j], \
							UextSeg,Eta1,Eta2,false);

		//Uwake at collocation using segment list
		//get collocation
		double Col[3] = {0.0,0.0,0.0};
		BilinearInterpTriad(Zeta[i][j], Zeta[i][j+1], \
							Zeta[i+1][j+1], Zeta[i+1][j], \
							Col,Eta1,Eta2,false);

		//Uwake
		BiotSavartSurf(ZetaStar_Vec,GammaStar_Vec,Col,\
					   0,0,\
					   VMOPTS.Mstar,VMOPTS.N,\
					   VMOPTS.Mstar,VMOPTS.N,\
					   VMOPTS.ImageMethod, UwakeSeg);

		//Usurf
		// set self influence to zero
		bool ignore = Seg[k].atImage;
		Seg[k].atImage = true;
		BiotSurfaceSegments(Seg,VMOPTS.M,VMOPTS.N,\
				   	   	   	VMOPTS.ImageMethod,Col,\
				   	   	   	UsurfSeg);
		Seg[k].atImage = ignore;

		//Combine Uext - ZetaDot
		SubTriad(UextSeg,ZetaDotSeg,UincSeg);
		//Add on Uwake
		AddTriad(UincSeg,UwakeSeg,UincSeg);
		//Add on Usurf
		AddTriad(UincSeg,UsurfSeg,UincSeg);

		//Get force
		Seg[k].Force(UincSeg,ForceSeg);

		//Map to corner points
		BilinearInterpTriad(Forces[i][j], Forces[i][j+1], \
							Forces[i+1][j+1], Forces[i+1][j], \
							ForceSeg,Eta1,Eta2,true);
	}

}


void KatzForces(const double* Zeta_Vec, const double* Gamma_Vec,\
				const double* ZetaStar_Vec, const double* GammaStar_Vec,\
				const double* ZetaDot_Vec, \
				const double* Uext_Vec, \
				VMopts VMOPTS,\
				const double* Gamma_tm1_Vec,\
				const double* Downwash,\
				const double* LiftVector_Vec,\
				double* Forces_Vec) {
	/** @brief Calculate panel forces at collocation points.
	 * @param Zeta_Vec Grid points.
	 * @param Gamma_Vec Bound vortex ring gammas.
	 * @param ZetaStar_Vec Wake grid points.
	 * @param GammaStar_Vec Wake vortex ring gammas.
	 * @param VMOPTS Simulation options.
	 * @param Gamma_tm1_Vec Bound vortex ring gammas from previous timestep.
	 * @param Forces_Vec Unsteady Forces at grid points.
	 */

	// cast vectors into useful pointers
	double (*Zeta)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Zeta_Vec;
	double (*ZetaDot)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) ZetaDot_Vec;
	double (*Uext)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Uext_Vec;
	double (*Gamma)[VMOPTS.N] = (double (*)[VMOPTS.N]) Gamma_Vec;
	double (*Gamma_tm1)[VMOPTS.N] = (double (*)[VMOPTS.N]) Gamma_tm1_Vec;
	double (*LiftVector)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) LiftVector_Vec;
	double (*Forces)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Forces_Vec;

	// declare local, automatically-managed dynamic memory using boost ...
	BoostArray2D alpha_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray3D Tau_c_(boost::extents[VMOPTS.M][VMOPTS.N][3]);
	BoostArray3D Tau_s_(boost::extents[VMOPTS.M][VMOPTS.N][3]);
	BoostArray2D Del_c_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray2D Del_s_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray2D Area_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray2D dGamma_dt_(boost::extents[VMOPTS.M][VMOPTS.N]);
	BoostArray3D Uwake_(boost::extents[VMOPTS.M][VMOPTS.N][3]);

	// ... and get useful pointers to that data
	double (*alpha)[VMOPTS.N] = (double (*)[VMOPTS.N]) alpha_.data();
	double (*Tau_c)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) Tau_c_.data();
	double (*Tau_s)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) Tau_s_.data();
	double (*Del_c)[VMOPTS.N] = (double (*)[VMOPTS.N]) Del_c_.data();
	double (*Del_s)[VMOPTS.N] = (double (*)[VMOPTS.N]) Del_s_.data();
	double (*Area)[VMOPTS.N] = (double (*)[VMOPTS.N]) Area_.data();
	double (*dGamma_dt)[VMOPTS.N] = (double (*)[VMOPTS.N]) dGamma_dt_.data();
	double (*Uwake)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) Uwake_.data();

	//loop through collocation points i,j
	#pragma omp parallel for
	for (unsigned int i = 0; i < VMOPTS.M; i++) {
		for (unsigned int j = 0; j < VMOPTS.N; j++) {

			//Temporary variables
			double ZetaDotCol[3] = {0.0,0.0,0.0};
			double UextCol[3] = {0.0,0.0,0.0};
			double Uinc[3] = {0.0,0.0,0.0};
			double Normal[3] = {0.0,0.0,0.0};
			double Collocation[3] = {0.0,0.0,0.0};
			double LiftVel[3] = {0.0,0.0,0.0};
			double TempC = 0.0;
			double TempS = 0.0;
			double TempGamma_i = 0.0;
			double TempGamma_j = 0.0;
			double DeltaP = 0.0;
			double LiftLocal = 0.0;
			double DragLocal1 = 0.0;
			double DragLocalWakeDown = 0.0;
			double DragLocal2 = 0.0;
			double LiftTemp[3] = {0.0,0.0,0.0};
			double DragTemp[3] = {0.0,0.0,0.0};
			double NormUinc[3] = {0.0,0.0,0.0};

			//calculate velocity of collocation
			BilinearMapTriad(ZetaDot[i][j], ZetaDot[i][j+1], \
							 ZetaDot[i+1][j+1], ZetaDot[i+1][j], \
							 ZetaDotCol);

			//External fluid velocities at collocation
			BilinearMapTriad(Uext[i][j], Uext[i][j+1], \
							 Uext[i+1][j+1], Uext[i+1][j], \
							 UextCol);

			// set incident velocity
			Uinc[0] = -ZetaDotCol[0] + UextCol[0];
			Uinc[1] = -ZetaDotCol[1] + UextCol[1];
			Uinc[2] = -ZetaDotCol[2] + UextCol[2];


			//calculate panel collocation
			BilinearMapTriad(Zeta[i][j], Zeta[i][j+1], \
							 Zeta[i+1][j+1], Zeta[i+1][j], \
							 Collocation);


			//calculate panel normals
			PanelNormal(Zeta[i][j], Zeta[i][j+1], \
						Zeta[i+1][j+1], Zeta[i+1][j],
								    Normal);


			// calculate panel tangential vectors
			PanelTau_c(Zeta[i][j], Zeta[i][j+1],\
					   Zeta[i+1][j+1], Zeta[i+1][j],
					   Tau_c[i][j]);

			PanelTau_s(Zeta[i][j], Zeta[i][j+1],\
					   Zeta[i+1][j+1], Zeta[i+1][j],
					   Tau_s[i][j]);


			// calculate local AoA
			alpha[i][j] = atan2(DotTriad(Uinc,Normal),\
							    DotTriad(Uinc,Tau_c[i][j]));


			// calculate panel \Delta chord and \Delta span
			Del_c[i][j] = PanelDeltaC(Zeta[i][j], Zeta[i][j+1],\
					   Zeta[i+1][j+1], Zeta[i+1][j]);

			Del_s[i][j] = PanelDeltaS(Zeta[i][j], Zeta[i][j+1],\
								   Zeta[i+1][j+1], Zeta[i+1][j]);


			// calculate panel areas
			Area[i][j] = PanelArea(Zeta[i][j], Zeta[i][j+1],\
								   Zeta[i+1][j+1], Zeta[i+1][j]);


			// calculate panel d(Gamma)/dt
			if (VMOPTS.Steady == true) {
				dGamma_dt[i][j] = 0.0;
			} else if (VMOPTS.Steady == false) {
				dGamma_dt[i][j] = (Gamma[i][j] - Gamma_tm1[i][j]) / \
										VMOPTS.DelTime;
			}


			// calculate wake induced velocity
			if (VMOPTS.Steady == true) {
				Uwake[i][j][0] = 0.0;
				Uwake[i][j][1] = 0.0;
				Uwake[i][j][2] = 0.0;
			} else if (VMOPTS.Steady == false) {
				BiotSavartSurf(ZetaStar_Vec, GammaStar_Vec, Collocation, \
						0, 0, \
						VMOPTS.Mstar, VMOPTS.N, \
						VMOPTS.Mstar,VMOPTS.N, VMOPTS.ImageMethod,\
						Uwake[i][j]);
			}


			// calculate local lift from Simpson (2013) plus external velocities
			// (-ZetaDot+Uext+Uwake)
			AddTriad(Uinc,Uwake[i][j],LiftVel);

			// dot product with tangential vectors
			TempC = DotTriad(LiftVel,Tau_c[i][j]);
			TempS = DotTriad(LiftVel,Tau_s[i][j]);

			//calculate spatial Delta Gamma_i
			if ( i==0 ) {
				TempGamma_i = Gamma[i][j];
			} else {
				TempGamma_i = Gamma[i][j] - Gamma[i-1][j];
			}

			//calculate spatial Delta Gamma_j
			// firstly for if there is an image plane at wing root
			if (VMOPTS.ImageMethod == 1) {
				if (j == 0) {
					TempGamma_j = Gamma[i][j];
				} else if (j > 0) {
					TempGamma_j = Gamma[i][j] - Gamma[i][j-1];
					//if the whole wing is modelled then
				} else if (VMOPTS.ImageMethod == 0) {
				//TODO: this.
				}
			}

			// pressure jump
			DeltaP = (TempC * TempGamma_i / Del_c[i][j] + \
					TempS * TempGamma_j / Del_s[i][j] + \
					dGamma_dt[i][j]);
//			printf("\tU.Tau_C:%f\tTempGamma_1:%f\tDel_c:%f\n",TempC,TempGamma_i,Del_c[i][j]);


			// calculate local lift
			LiftLocal = DeltaP * Area[i][j] * cos(alpha[i][j]);
//			printf("\tDeltaP:%f \tArea:%f\talpha:%f\n",DeltaP,Area[i][j],alpha[i][j]);


			// calculate local drag
			if (VMOPTS.Steady == true) {
				DragLocal1 = -Downwash[i*VMOPTS.N + j]*TempGamma_i*Del_s[i][j];
			} else if (VMOPTS.Steady == false) {
				//wake velocity down local lift must be added
				DragLocalWakeDown = DotTriad(Uwake[i][j],LiftVector[i][j]);
				DragLocal1 = -(Downwash[i*VMOPTS.N + j] + DragLocalWakeDown) * \
							 TempGamma_i*Del_s[i][j];
			}

			DragLocal2 = dGamma_dt[i][j] * Area[i][j] * sin(alpha[i][j]);
//			printf("\tdGamma_dt%f\tdrag1:%f\tdrag2:%f\n",dGamma_dt[i][j],DragLocal1,DragLocal2);

			// total force
			//lift vector
			MulTriad(LiftVector[i][j],LiftLocal,LiftTemp);
//			printf("lift vector = ");
//			PrintTriad(LiftVector[i][j]);
//			PrintTriad(LiftTemp);
//			printf("\n");

			//drag vector
			// subtract wake velocity (this is set to zero above for steady)
			NormaliseTriad(Uinc,NormUinc);
			if (NormTriad(Uinc) == 0.0) {
				std::cerr << "VLM: Warning incident velocity Zero!" \
										  << std::endl;
			}
			double Foo = DragLocal1+DragLocal2;
			MulTriad(NormUinc,Foo,DragTemp);
//			printf("drag vector = ");
//			PrintTriad(NormUinc);
//			PrintTriad(DragTemp);
//			printf("\n");

			//combine and save
			AddTriad(LiftTemp,DragTemp,LiftTemp);

//			printf("\tlift:%f\tdrag:%f\tVector:",LiftLocal,DragLocal1+DragLocal2);
//			PrintTriad(LiftTemp);
//			printf("\n");

			// map panel force (on LE segment) to nodes
			#pragma omp critical
			BilinearInterpTriad(Forces[i][j],Forces[i][j+1],\
								Forces[i+1][j+1],Forces[i+1][j],\
								LiftTemp,\
								0.0,0.5,true);
		}
	}
}


void cpp_solver_vlm(const double* Zeta_Vec, const double* ZetaDot_Vec, \
				    const double* Uext_Vec, double* ZetaStar_Vec, \
				    VMopts VMOPTS, \
				    double* Forces_Vec, \
				    double* Gamma_Vec, double* GammaStar_Vec,\
				    double* AIC_Vec,\
				    double* BIC_Vec) {
	/**@warning If VMOPTS.NewAIC == false and dynamic solutions are sought where
	 * the local lift vector may be changing the Katz and Plotkin force method
	 * may return wrong drag as new BIC is required.
	 */

	/* Calculate and save panel normal vectors, location of collocation point,
	 * velocity at collocation point, and the External fluid velocity at
	 * collocation points (using bilinear mapping on grid) and relative
	 *  normal-wash at collocation points.
	 *  Use this to calculate BIC then AIC matrix if desired.
	 *  Extract forces.
	 *  Rollup if desired. */

	/** Here the contiguous data provided by Python - 'const double* Zeta_Vec' -
	 * needs to be accessed as if it were a 3D C-array -'double Zeta[M][N][3]'.
	 * To do this WITHOUT copying any data we create an appropriate pointer
	 * to the 3D array
	 * - 'double (*Zeta)[VMOPTS_N+1][3]' - that is initialised by 'casting'
	 * the contiguous data pointer 'Zeta_Vec' to the type required for our
	 * 3D array - '(*Zeta)[VMOPTS.N+1][3]'.
	 */
	double (*Zeta)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Zeta_Vec;
	double (*ZetaDot)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) ZetaDot_Vec;
	double (*ZetaStar)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) ZetaStar_Vec;
	double (*Uext)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Uext_Vec;
	double (*Forces)[VMOPTS.N+1][3] = (double (*)[VMOPTS.N+1][3]) Forces_Vec;

	//Declare dynamic memory for local variables using boost multi_array...
	/**@warning Once this function ends this data does not exist anymore */
	BoostArray3D ZetaCol_(boost::extents[VMOPTS.M][VMOPTS.N][3]);
	BoostArray3D ZetaDotCol_(boost::extents[VMOPTS.M][VMOPTS.N][3]);
	BoostArray3D UextCol_(boost::extents[VMOPTS.M][VMOPTS.N][3]);
	BoostArray3D Normal_(boost::extents[VMOPTS.M][VMOPTS.N][3]);
	BoostArray3D Uinc_(boost::extents[VMOPTS.M][VMOPTS.N][3]);
	BoostArray3D LocalLift_(boost::extents[VMOPTS.M][VMOPTS.N][3]);


	// ... and useful array pointers to that memory.
	double (*ZetaCol)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) ZetaCol_.data();
	double (*ZetaDotCol)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) ZetaDotCol_.data();
	double (*UextCol)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) UextCol_.data();
	double (*Normal)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) Normal_.data();
	double (*Uinc)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) Uinc_.data();
	double (*LocalLift)[VMOPTS.N][3] = (double (*)[VMOPTS.N][3]) LocalLift_.data();


	// Declare Eigen types for linear algebra
	Eigen::VectorXd RHS(VMOPTS.M*VMOPTS.N);
	EigenMapVectXd Gamma(Gamma_Vec,VMOPTS.M*VMOPTS.N);
	EigenMapVectXd GammaStar(GammaStar_Vec,VMOPTS.Mstar*VMOPTS.N);
	EigenMapMatrixXd Aic(AIC_Vec,VMOPTS.M*VMOPTS.N,VMOPTS.M*VMOPTS.N);
	EigenMapMatrixXd BIC(BIC_Vec,VMOPTS.M*VMOPTS.N,VMOPTS.M*VMOPTS.N);
	Eigen::VectorXd Downwash(VMOPTS.M*VMOPTS.N);


	//set OMP environment
	omp_set_num_threads(VMOPTS.NumCores);

	//use Eigen to set GammaStar_Vec to 1.0 for all j
	if (VMOPTS.Steady == true) {
		GammaStar.setOnes();
	}


	//loop through collocation points i,j
	#pragma omp parallel for
	for (unsigned int i = 0; i < VMOPTS.M; i++) {
		for (unsigned int j = 0; j < VMOPTS.N; j++) {

			//temporary variables
			double GammaOne = 1.0;
			double Temp1[3] = {0.0,0.0,0.0};
			double Temp2[3] = {0.0,0.0,0.0};
			double Temp3[3] = {0.0,0.0,0.0};
			Eigen::Matrix3d Proj;
			Eigen::Matrix3d Eye = Eigen::Matrix3d::Identity();
			Eigen::Vector3d Uincident(0.0,0.0,0.0);
			Eigen::Vector3d LocalLiftEig(0.0,0.0,0.0);
			Eigen::Vector3d NormalEig(0.0,0.0,0.0);

			//collocation points
			BilinearMapTriad(Zeta[i][j], Zeta[i][j+1], \
							 Zeta[i+1][j+1], Zeta[i+1][j], \
							 ZetaCol[i][j]);

			//relative inertial velocity at collocation points
			BilinearMapTriad(ZetaDot[i][j], ZetaDot[i][j+1], \
							 ZetaDot[i+1][j+1], ZetaDot[i+1][j], \
							 ZetaDotCol[i][j]);

			//External fluid velocities at collocation
			BilinearMapTriad(Uext[i][j], Uext[i][j+1], \
							 Uext[i+1][j+1], Uext[i+1][j], \
							 &UextCol[i][j][0]);

			//panel normals
			PanelNormal(Zeta[i][j], Zeta[i][j+1], \
					    Zeta[i+1][j+1], Zeta[i+1][j],
					    Normal[i][j]);

			//Normal wash at Collocations points
			SubTriad(UextCol[i][j],ZetaDotCol[i][j],Uinc[i][j]);
//			printf("Uext - ZetaDot = ");
//			PrintTriad(Uinc[i][j]);

			if (VMOPTS.Steady == false) {
				BiotSavartSurf(ZetaStar_Vec, GammaStar_Vec, ZetaCol[i][j],\
								0,0,\
								VMOPTS.Mstar,VMOPTS.N,\
								VMOPTS.Mstar,VMOPTS.N,\
								VMOPTS.ImageMethod, Temp1);
				AddTriad(Uinc[i][j],Temp1,Uinc[i][j]);
			}

			// Fill RHS Vector
			RHS(i*VMOPTS.N + j) = - DotTriad(Uinc[i][j],
									 	 	 Normal[i][j]);

			// calc local lift vector for BIC calculation (as in Simpson, 2013)
			// calc orthogonal project operator (Eigen types)
			Uincident(0) = -ZetaDotCol[i][j][0] + Uext[i][j][0];
			Uincident(1) = -ZetaDotCol[i][j][1] + Uext[i][j][1];
			Uincident(2) = -ZetaDotCol[i][j][2] + Uext[i][j][2];

			//check incident velocity and define projection operator
			if (Uincident.norm() == 0.0) {
				std::cerr << "VLM: Warning incident velocity Zero!" \
						  << std::endl;
				Proj = Eye;
			} else {
				Proj = Eye - ( Uincident.normalized() ) * \
							( Uincident.normalized().transpose() );
			} // END if, else

			// get local lift vector
			NormalEig(0) = Normal[i][j][0];
			NormalEig(1) = Normal[i][j][1];
			NormalEig(2) = Normal[i][j][2];
			LocalLiftEig = Proj*NormalEig;
			LocalLiftEig.normalize();
			// populate local lift LocalLift
			LocalLift[i][j][0] = LocalLiftEig[0];
			LocalLift[i][j][1] = LocalLiftEig[1];
			LocalLift[i][j][2] = LocalLiftEig[2];

			// Loop through each vortex ring (panel) ii, jj
			// TODO: if this is to be parallelised write to variables local to
			// the scope NOT the global pointers!

			if (VMOPTS.NewAIC == true) {
				for (unsigned int ii = 0; ii < VMOPTS.M; ii++) {
					for (unsigned int jj = 0; jj < VMOPTS.N; jj++) {
						//Induced vel of ring ii,jj to panel i,j

						// reset running total to zero
						Temp3[0] = 0.0;
						Temp3[1] = 0.0;
						Temp3[2] = 0.0;


						// chordwise orientated vorticity only (for BIC)

						// segment 2 (point 2 -> point 3)
						C_BiotSegment(ZetaCol[i][j],Zeta[ii][jj+1],Zeta[ii+1][jj+1],
									  GammaOne, Temp1);

						// segment 4 (point 4 -> point 1)
						C_BiotSegment(ZetaCol[i][j],Zeta[ii+1][jj],Zeta[ii][jj],
									  GammaOne, Temp2);

						// sum of segments 2 and 4
						AddTriad(Temp1,Temp2,Temp3); //Temp3 overwritten here

						//image method
						if (VMOPTS.ImageMethod == 1) {
							// segment 2 (point 2 -> point 3) image
							C_BiotSegment_ImageYZ(ZetaCol[i][j],\
												  Zeta[ii][jj+1],Zeta[ii+1][jj+1], \
												  GammaOne, Temp1);

							AddTriad(Temp3,Temp1,Temp3);


							// segment 4 (point 4 -> point 1) image
							C_BiotSegment_ImageYZ(ZetaCol[i][j],\
												  Zeta[ii+1][jj],Zeta[ii][jj],\
												  GammaOne, Temp1);

							AddTriad(Temp3,Temp1,Temp3);
						}


						// if we're at the trailing edge then segment3 must be added
						// also the wake effect is added here
						if (ii == VMOPTS.M-1) {
							// add TE segment
							C_BiotSegment(ZetaCol[i][j],\
										  Zeta[ii+1][jj+1],\
										  Zeta[ii+1][jj],\
										  GammaOne,\
										  Temp1);

							AddTriad(Temp3,Temp1,Temp3);

							//image method
							if (VMOPTS.ImageMethod == 1) {
								C_BiotSegment_ImageYZ(ZetaCol[i][j],\
													  Zeta[ii+1][jj+1],\
													  Zeta[ii+1][jj],\
													  GammaOne, Temp1);
								AddTriad(Temp3,Temp1,Temp3);
							} // END if Image Method

							// add wake effect
							if (VMOPTS.Steady == true) {
								BiotSavartSurf(ZetaStar_Vec, GammaStar_Vec, \
											   ZetaCol[i][j],\
											   0, jj, \
											   VMOPTS.Mstar, jj+1, \
											   VMOPTS.Mstar, VMOPTS.N, VMOPTS.ImageMethod,\
											   Temp1);
								AddTriad(Temp3,Temp1,Temp3);
							}
						} // END if VMOPTS.M-1


						// element of BIC matrix
						BIC(i*VMOPTS.N + j,ii*VMOPTS.N + jj) = \
								DotTriad(Temp3,LocalLift[i][j]);

	//					printf("BIC velocity = ");
	//					PrintTriad(Temp3);
	//					printf("\n");


						//calculate induced velocity for AIC calculation
						// both spanwise-orientated segments must be added UNLESS
						// we are at the trailing edge. At the TE we add the only
						// remaining bound segment.

						if (ii < VMOPTS.M-1) {
							//add segments 1 and 3
							C_BiotSegment(ZetaCol[i][j],Zeta[ii][jj],Zeta[ii][jj+1],
														GammaOne, Temp1);
							C_BiotSegment(ZetaCol[i][j],Zeta[ii+1][jj+1],Zeta[ii+1][jj],
														GammaOne, Temp2);
							AddTriad(Temp3,Temp1,Temp3);
							AddTriad(Temp3,Temp2,Temp3);

							if (VMOPTS.ImageMethod == 1) {
								//add segments 1 and 3 image
								C_BiotSegment_ImageYZ(ZetaCol[i][j],\
													  Zeta[ii][jj],\
													  Zeta[ii][jj+1],\
													  GammaOne,\
													  Temp1);

								C_BiotSegment_ImageYZ(ZetaCol[i][j],\
													  Zeta[ii+1][jj+1],\
													  Zeta[ii+1][jj],\
													  GammaOne,\
													  Temp2);
								AddTriad(Temp3,Temp1,Temp3);
								AddTriad(Temp3,Temp2,Temp3);

							} // END if Image Method

						} else if (ii == VMOPTS.M-1) {
							//add segments 1 only
							C_BiotSegment(ZetaCol[i][j],Zeta[ii][jj],Zeta[ii][jj+1],
														GammaOne, Temp1);

							AddTriad(Temp3,Temp1,Temp3);
							if (VMOPTS.ImageMethod == 1) {
								//add segment 1 image
								C_BiotSegment_ImageYZ(ZetaCol[i][j],\
													  Zeta[ii][jj],\
													  Zeta[ii][jj+1],
													  GammaOne, Temp1);

								AddTriad(Temp3,Temp1,Temp3);

							} // END if Image Method

						} //END if, else if VMOPTS.M-1

						// element of AIC matrix
						Aic(i*VMOPTS.N + j,ii*VMOPTS.N + jj) = \
													DotTriad(Temp3,Normal[i][j]);

					} //END for jj
				} //END for ii
			} //END if New AIC
		} //END for j
	} //END for i


	// declare memory for Gamma_tm1_
	BoostArray2D Gamma_tm1_(boost::extents[VMOPTS.M][VMOPTS.N]);

	// copy Gamma_Vec to Gamma_tm1
	unsigned int k = 0;
	for (unsigned int i = 0; i < VMOPTS.M; i++) {
		for (unsigned int j = 0; j < VMOPTS.N; j++) {
			k = i*VMOPTS.N + j;
			Gamma_tm1_[i][j] = Gamma(k);
		}
	}


	bool writeDebug = false;
	if (writeDebug == true) {
		std::ofstream myfile;
		myfile.open("/home/rjs10/git/SHARPy/output/temp/rhs.dat");
		myfile << RHS;
		myfile.close();
	}

	if (false) {// I'm testing the new AIC routine
		AIC(Zeta_Vec,VMOPTS.M,VMOPTS.N,Zeta_Vec,VMOPTS.M,VMOPTS.N,false,AIC_Vec);
		// must add wake effect in steady calcs!
	}

	//solve for gamma
//	std::cout << "AIC matrix:" << std::endl << Aic << std::endl;
	Gamma = Aic.colPivHouseholderQr().solve(RHS);

	//set GammaStar to TE velocities
	//check if steady set all to TE gamma
	if (VMOPTS.Steady == 1 && VMOPTS.Mstar == 1) {
		//wake gamma set very easily
		GammaStar = Gamma.tail(VMOPTS.N);
	} else if (VMOPTS.Steady == 1 && VMOPTS.Mstar > 1) {
		//each column must be set
		for (unsigned int i = 0; i < VMOPTS.Mstar; i++) {
			GammaStar.segment(i*VMOPTS.N,VMOPTS.N) = Gamma.tail(VMOPTS.N);
		}
	}


	//calculate downwash
	Downwash = BIC * (Gamma);

//	std::cout << Gamma << std::endl;
//
//
//	std::cout << AIC << std::endl << std::endl;
//	std::cout << BIC << std::endl;



	// Calculate forces
	double* LocalLift_Vec = LocalLift_.data();
	double* Downwash_ptr = Downwash.data();
	double* Gamma_tm1_vec = Gamma_tm1_.data();

	if (VMOPTS.KJMeth == 0) {
		KatzForces(Zeta_Vec, Gamma_Vec, ZetaStar_Vec, GammaStar_Vec,\
				   ZetaDot_Vec, Uext_Vec, VMOPTS, Gamma_tm1_vec, Downwash_ptr,\
				   LocalLift_Vec, Forces_Vec);

	} else if (VMOPTS.KJMeth == 1) {
		KJMethodForces(Zeta_Vec, Gamma_Vec, ZetaStar_Vec, GammaStar_Vec,\
				   	   ZetaDot_Vec, Uext_Vec, VMOPTS, Gamma_tm1_vec,\
				   	   Forces_Vec);
	}



	if (VMOPTS.Rollup == false) {
		/* do nothing (external velocities at the wake corner points must be
		 * accounted for outside cpp solver for now)
		 */
	} else if (VMOPTS.Rollup == true) {
		/* calculate rollup at wake grid corner points*/

		//loop through wake corner points
		#pragma omp parallel for
		for (unsigned int i = 0; i < VMOPTS.Mstar + 1; i++) {
			for (unsigned int j = 0; j < VMOPTS.N + 1; j++) {
				double Vel1[3] = {0.0,0.0,0.0};
				double Vel2[3] = {0.0,0.0,0.0};

				//surface contribution
				BiotSavartSurf(Zeta_Vec, Gamma_Vec, ZetaStar[i][j],\
							   0,0,\
							   VMOPTS.M,VMOPTS.N,\
							   VMOPTS.M,VMOPTS.N,\
							   VMOPTS.ImageMethod,\
							   Vel1);

				//wake contribution
				BiotSavartSurf(ZetaStar_Vec, GammaStar_Vec, ZetaStar[i][j],\
							   0,0,\
							   VMOPTS.Mstar,VMOPTS.N,\
							   VMOPTS.Mstar,VMOPTS.N,\
							   VMOPTS.ImageMethod,\
							   Vel2);

				//Add and multiply by timestep then add to existing ZetaStar
				AddTriad(Vel1,Vel2,Vel1);
				MulTriad(Vel1,VMOPTS.DelTime,Vel1);
				AddTriad(ZetaStar[i][j],Vel1,ZetaStar[i][j]);
			}
		}

	}
}

