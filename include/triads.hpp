/*
 * triads.hpp
 *
 *  Created on: 31 Jan 2013
 *      Author: rjs10
 */

#ifndef TRIADS_HPP_
#define TRIADS_HPP_

void AddTriad(const double* x1, const double* x2, double* xOut);

void SubTriad(const double* x1, const double* x2, double* xOut);

void MulTriad(const double* x1, const double Factor, double* xOut);

void DivTriad(const double* x1, const double Factor, double* xOut);

double DotTriad(const double* x1, const double* x2);

double NormTriad(const double* x1);

void NormaliseTriad(const double* x1, double* xOut);

void CrossTriad(const double* x1, const double* x2, double* xOut);

void BilinearMapTriad(const double* p1, const double* p2, \
					  const double* p3, const double* p4, \
					  double* pOut);

void BilinearInterpTriad( double* x00, double* x01, \
					  	  double* x11, double* x10, \
					  	  double* xP,\
					  	  const double eta1, const double eta2,\
					  	  const bool reverse);

void CopyTriad(double* pTarget, const double* pSrc);

void PrintTriad(const double* x);

void AssignTriad(double* x,
		           const double x0, const double x1, const double x2);

#endif /* TRIADS_HPP_ */
