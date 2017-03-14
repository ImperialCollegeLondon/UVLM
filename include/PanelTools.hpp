/*
 * PanelTools.hpp
 *
 *  Created on: 13 Feb 2013
 *      Author: rjs10
 */

#ifndef PANELTOOLS_HPP_
#define PANELTOOLS_HPP_

#include<Eigen/Dense>

void PanelNormal(const double* p1, const double* p2, \
				 const double* p3, const double* p4, \
				 double* normal);

void PanelTau_c(const double* p1, const double* p2, \
				 const double* p3, const double* p4, \
				 double* Tau_c);

void PanelTau_s(const double* p1, const double* p2, \
				 const double* p3, const double* p4, \
				 double* Tau_s);


double PanelDeltaC(const double* p1, const double* p2, \
				 const double* p3, const double* p4);

double PanelDeltaS(const double* p1, const double* p2, \
				 const double* p3, const double* p4);


double PanelArea(const double* p1, const double* p2, \
				 const double* p3, const double* p4);

void XiKernel(const unsigned int k,
		      const unsigned int q,
		      const unsigned int N,
		      const double eta1,
		      const double eta2,
		      const Eigen::Matrix3d& XiKern_);

void hKernel(const unsigned int q,
			  const unsigned int k,
			  const unsigned int l,
			  const unsigned int N,
			  const Eigen::Matrix3d& hKern_);

#endif /* PANELTOOLS_HPP_ */
