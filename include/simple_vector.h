/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#ifndef vector_h
#define vector_h

#include <vector>
#include <iostream>

namespace XMOF2D {

void   daxpy(double alpha, const std::vector<double>& x, std::vector<double>& y);
void   dscal(double alpha, std::vector<double>& x);
double ddot(const std::vector<double>& x, const std::vector<double>& y);
double dnrm2(const std::vector<double>& x);

void   daxpy(double alpha, const double* x, double* y, int n);
void   dscal(double alpha, double* x, int n);
double ddot(const double* x, const double* y, int n);

std::ostream& operator<<(std::ostream& os, const std::vector<double>& v);

}
#endif /* vector_h */
