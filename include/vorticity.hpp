/*
 * vorticity.hpp
 *
 *  Created on: 28 Jan 2013
 *      Author: rjs10
 */

#ifndef VORTICITY_HPP_
#define VORTICITY_HPP_

#include <Eigen/Dense>
#include <vector>
#include <datatypesx.hpp>

Eigen::Vector3d BiotSegment(const Eigen::Vector3d& xp, \
							const Eigen::Vector3d& x1, \
							const Eigen::Vector3d& x2, \
							const double& gam);

void BiotSegment_map(const Eigen::Map<Eigen::Vector3d>& xp, \
								const Eigen::Map<Eigen::Vector3d>& x1, \
								const Eigen::Map<Eigen::Vector3d>& x2, \
								const double& gam,
								const Eigen::Map<Eigen::Vector3d>& Uind_triad);

void C_BiotSegment(const double* xP,\
					  const double* x1,\
					  const double* x2,\
					  const double& Gamma,\
					  double* Uind);

void C_BiotSegment_ImageYZ(const double* xP,\
		  const double* x1,\
		  const double* x2,\
		  double& Gamma,\
		  double* Uind);


void BiotSavartSurf(const double* Zeta_Vec, const double* Gamma_Vec, \
					const double* TargetTriad,\
					const unsigned int Mstart, const unsigned int Nstart,
					const unsigned int Mend, const unsigned int Nend,
					const unsigned int Mfull, const unsigned int Nfull, \
					const bool ImageMethod, \
					double* Uout);

void BiotSurfaceSegments(std::vector<VortexSegment>& Seg,\
						 unsigned int SurfM,\
						 unsigned int SurfN,\
						 bool ImageMethod,\
						 double* Target,\
						 double* Uind);


#endif /* VORTICITY_HPP_ */
