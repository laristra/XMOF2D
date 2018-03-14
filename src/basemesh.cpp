/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Los Alamos National Security, LLC.
 All rights reserved.
*/

#include "mesh.h"
#include "exception.h"

namespace XMOF2D {

BaseMesh::BaseMesh(const MeshConfig& config, const IRTolerances& tolerances) {
  dist_eps_ = tolerances.dist_eps;
  ddot_eps_ = tolerances.ddot_eps;
  div_eps_ = tolerances.div_eps;
  vfrac_eps_ = tolerances.vfrac_eps;
  ang_eps_ = tolerances.ang_eps;
  mof_max_iter_ = tolerances.mof_max_iter;
  
  cells_material = config.cells_material;
  
  int nnodes = (int) config.nodes_coords.size();
  for (int inode = 0; inode < nnodes; inode++)
    nodes.push_back(Node(config.nodes_coords[inode], inode, *this));
  
  int nfaces = (int) config.ifaces_nodes.size();
  for (int iface = 0; iface < nfaces; iface++)
    faces.push_back(Face({config.ifaces_nodes[iface][0], config.ifaces_nodes[iface][1]}, iface, *this));
    
  int ncells = (int) config.icells_faces.size();
  for (int icell = 0; icell < ncells; icell++) {
    int nsides = (int) config.icells_faces[icell].size();
    for (int iside = 0; iside < nsides; iside++) {
      int iface = config.icells_faces[icell][iside];
      if (get_face(iface).get_cell_index(false) == -1)
        faces[iface].set_cell_index(false, icell, true);
      else if (get_face(iface).get_cell_index(true) == -1)
        faces[iface].set_cell_index(true, icell, true);
      else {
        THROW_EXCEPTION("Face " << iface << " is shared by more than two cells!");
      }
    }
    std::pair< std::vector<int>, std::vector<bool> > ccw_faces_ordering =
      get_ccw_faces_ordering(config.icells_faces[icell]);
    const std::vector<int>& icellfaces = ccw_faces_ordering.first;
    for (int iside = 0; iside < nsides; iside++) {
      Face& cur_face = faces[icellfaces[iside]];
      if (cur_face.normal_is_out(icell) != ccw_faces_ordering.second[iside])
        cur_face.swap_nodes_order();
    }
    cells.push_back(Cell(icellfaces, icell, *this));
  }

  set_nodes_ifaces();
  set_nodes_icells();
}

void BaseMesh::ConstructMinimesh(int icell, const std::vector<int>& imats, const std::vector<double>& vfracs, const std::vector<Point2D>& centroids) {
  cells[icell].delete_minimesh();
  if (imats.size() > 1) {
    cells[icell].ConstructMinimesh(imats, vfracs, centroids);
    cells_material[icell] = -1;
  }
  else cells_material[icell] = imats[0];
}

}
