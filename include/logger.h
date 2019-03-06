/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#ifndef logger_h
#define logger_h

#include <vector>
#include <iostream>

namespace XMOF2D {

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
  os << " size = " << v.size() << std::endl;
  for (typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); it++)
    os << " " << it - v.begin() << ": " << *it << std::endl;
  return os;
}

template<typename T1,typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1,T2>& p) {
  return os << "(" << p.first << "," << p.second << ")";
}
}
#endif /* logger_h */
