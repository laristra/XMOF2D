/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#include "xmof2D.h"

namespace XMOF2D {

XMOF_Reconstructor::XMOF_Reconstructor(const MeshConfig& mesh_config, const IRTolerances& ir_tolerances) : mesh(mesh_config, ir_tolerances) {}

void XMOF_Reconstructor::set_materials_data(const CellsMatData& mat_data) {
  XMOF2D_ASSERT(mat_data.cells_materials.size() == mesh.ncells(),
    "Provided material data is for a mesh with a different number of cells!");
  for (int icell = 0; icell < mesh.ncells(); icell++) {
    int ncmats = (int) mat_data.cells_vfracs[icell].size();
    for (int icmat = 0; icmat < ncmats; icmat++) {
      double cmat_vol = mat_data.cells_vfracs[icell][icmat]*mesh.get_cell(icell).size();
      XMOF2D_ASSERT(cmat_vol > mesh.area_eps(), "Given area for material " << 
        mat_data.cells_materials[icell][icmat] << " in cell " << icell <<
        " is " << cmat_vol << ", which should be above the area tolerance of " << mesh.area_eps());
    }
  }

  cells_materials = mat_data.cells_materials;
  cells_vfracs = mat_data.cells_vfracs;
  cells_centroids = mat_data.cells_centroids;
}

void XMOF_Reconstructor::set_global_ind(const GlobalIndData& global_ind_data) {
  nodes_global_ind = global_ind_data.nodes_global_ind;
  faces_global_ind = global_ind_data.faces_global_ind;
  cells_global_ind = global_ind_data.cells_global_ind;
}

const std::vector<int>& XMOF_Reconstructor::get_cell_materials(int icell) const {
  return cells_materials[icell];
}
const std::vector<double>& XMOF_Reconstructor::get_cell_vfracs(int icell) const {
  return cells_vfracs[icell];
}
const std::vector<Point2D>& XMOF_Reconstructor::get_cell_centroids(int icell) const {
  return cells_centroids[icell];
}

int XMOF_Reconstructor::construct_minimesh(int icell) {
  try {
    mesh.ConstructMinimesh(icell, cells_materials[icell], cells_vfracs[icell], cells_centroids[icell]);
  }
  catch (Exception e) {
    std::cerr << "Interface reconstruction for cell " << icell << 
      " has failed with the following message" << std::endl <<
      e.what() << std::endl;
    return 1;
  } 

  return 0;
}

int XMOF_Reconstructor::cell_global_ind(int icell) const { 
  XMOF2D_ASSERT_SIZE_LESS(icell, mesh.ncells());
  return cells_global_ind[icell]; 
}

int XMOF_Reconstructor::face_global_ind(int iface) const { 
  XMOF2D_ASSERT_SIZE_LESS(iface, mesh.nfaces());
  return faces_global_ind[iface];
}
int XMOF_Reconstructor::node_global_ind(int inode) const { 
  XMOF2D_ASSERT_SIZE_LESS(inode, mesh.nnodes());
  return nodes_global_ind[inode];
}

int XMOF_Reconstructor::cell_parent_global_ind(const Cell& cell) const {
  int parent_ind = cell.iparent();
  return (parent_ind == -1) || (cells_global_ind.empty()) ? parent_ind : cells_global_ind[parent_ind];
}

int XMOF_Reconstructor::face_parent_global_ind(const Face& face) const {
  int parent_ind = face.iparent();
  return (parent_ind == -1)  || (faces_global_ind.empty()) ? parent_ind : faces_global_ind[parent_ind];
}

int XMOF_Reconstructor::node_parent_global_ind(const Node& node) const {
  int parent_ind = node.iparent();
  return (parent_ind == -1)  || (nodes_global_ind.empty()) ? parent_ind : nodes_global_ind[parent_ind];
}

}
