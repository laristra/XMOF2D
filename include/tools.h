/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#ifndef tools_h
#define tools_h

#include <limits>
#include <vector>
#include <cmath>

namespace XMOF2D {

static const double epsilon = std::numeric_limits<double>::epsilon();
static const double PI = 3.14159265358979;

inline bool is_equal(double x, double y) {
  return (std::fabs(x - y) <= epsilon) ? true : false;
}

inline bool is_not_equal(double x, double y) {
  return (std::fabs(x - y) > epsilon) ? true : false;
}

template<typename T>
std::vector< std::pair<int,int> > isamei(const std::vector<T>& a, const std::vector<T>& b) {
  std::vector< std::pair<int,int> > isame_ind;
  for (int ia = 0; ia < a.size(); ia++)
    for (int ib = 0; ib < b.size(); ib++)
      if (a[ia] == b[ib])
        isame_ind.push_back(std::make_pair(ia, ib));
  
  return isame_ind;
}

template<typename T>
inline T pow2(const T& x) {
  return x*x;
}

}
#endif /* tools_h */
