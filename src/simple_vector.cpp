/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#include "simple_vector.h"
#include "tools.h"
#include "exception.h"

namespace XMOF2D {

void daxpy(double alpha, const double* x, double* y, int n) {
  for (int i = 0; i < n; i++)
    y[i] += alpha*x[i];
}

void dscal(double alpha, double* x, int n) {
  for (int i = 0; i < n; i++)
    x[i] *= alpha;
}

double ddot(const double* x, const double* y, int n) {
  double s = 0.;
  for (int i = 0; i < n; i++)
    s += x[i]*y[i];
  return s;
}

void daxpy(double alpha, const std::vector<double>& x, std::vector<double>& y) {
  int n = (int) x.size();
  XMOF2D_ASSERT(y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");
  XMOF2D_ASSERT(n, "Zero length vectors");
  
  daxpy(alpha, &x[0], &y[0], n);
}

void dscal(double alpha, std::vector<double>& x) {
  int n = (int) x.size();
  
  if (n && is_not_equal(alpha, 1.0))
    dscal(alpha, &x[0], n);
}

double ddot(const std::vector<double>& x, const std::vector<double>& y) {
  int n = (int) x.size();
  XMOF2D_ASSERT(y.size() == n, "Different sizes: x (" << x.size() << "), y (" << y.size() << ")");
  XMOF2D_ASSERT(n, "Zero length vectors");
  
  return ddot(&x[0], &y[0], n);
}

double dnrm2(const std::vector<double>& x) {
  return sqrt(ddot(x, x));
}

std::ostream& operator<<(std::ostream& os, const std::vector<double>& v) {
  const int n = (int) v.size();
  os << "Size = " << n << std::endl;
  for (int i = 0; i < n; i++)
    os << "   " << i << ": " << v[i] << std::endl;
  return os;
}

}
