/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Los Alamos National Security, LLC.
 All rights reserved.
*/

#include "mesh.h"
#include <algorithm>

namespace XMOF2D {

Factory::Factory(const MeshBC& _mesh, std::vector<Face>& _faces) : mesh(_mesh), faces(_faces) {
  typedef std::pair< std::vector<int>, int > pair;
  for (int iface = 0; iface < faces.size(); iface++) {
    std::vector<int> r = faces[iface].get_nodes();
    std::sort(r.begin(), r.end());
    int rmin = r[0];
    if (rmin >= sorted.size()) {
      sorted.resize(rmin + 1);
      sorted[rmin].push_back(pair(r, iface));
      continue;
    }
    std::vector<pair>& rs = sorted[rmin];
    int j;
    for (j = 0; j < rs.size(); j++)
      if (rs[j].first == r)
        break;
    if (j == rs.size())
      rs.push_back(pair(r, iface));
  }
}

int Factory::get_face(int cell_ind, const std::vector<int>& v) {
  /* Construct r = sorted v */
  std::vector<int> r = v;
  std::sort(r.begin(), r.end());
  
  typedef std::pair< std::vector<int>, int > pair;
  
  /* Check if such smallest index was already in sorted */
  int rmin = r[0];
  if (rmin >= sorted.size()) {
    /* Create new face */
    int nface = (int) faces.size();
    faces.push_back(Face(v, nface, mesh));
    faces[nface].set_cell_index(false, cell_ind, true);
    
    sorted.resize(rmin + 1);
    sorted[rmin].push_back(pair(r, nface));
    
    return nface;
  }
  
  std::vector<pair>& rs = sorted[rmin];
  int j;
  for (j = 0; j < rs.size(); j++)
    if (rs[j].first == r)
      break;
  if (j == rs.size()) {
    /* Create new face */
    int nface = (int) faces.size();
    faces.push_back(Face(v, nface, mesh));
    faces[nface].set_cell_index(false, cell_ind, true);
    
    rs.push_back(pair(r, nface));
    
    return nface;
  }
  
  /* old face */
  faces[rs[j].second].set_cell_index(true, cell_ind, true);
  return rs[j].second;
}

}
