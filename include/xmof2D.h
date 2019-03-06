/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#ifndef xmof_ir_h
#define xmof_ir_h

#include "mesh.h"

namespace XMOF2D {

struct CellsMatData {
  std::vector< std::vector<int> >    cells_materials; //indices of materials in each cell
  std::vector< std::vector<double> > cells_vfracs;    //volume fractions of materials in each cell
  std::vector< std::vector<Point2D> >  cells_centroids; //centroids of materials in each cells
};

struct GlobalIndData {
  std::vector<int> nodes_global_ind;  //indices of base mesh nodes in the global mesh
  std::vector<int> faces_global_ind;  //indices of base mesh faces in the global mesh
  std::vector<int> cells_global_ind;  //indices of base mesh cells in the global mesh
};

class XMOF_Reconstructor {
private:
  BaseMesh mesh;
  
  std::vector< std::vector<int> >     cells_materials;
  std::vector< std::vector<double> >  cells_vfracs;
  std::vector< std::vector<Point2D> > cells_centroids;
  
  std::vector<int> nodes_global_ind;
  std::vector<int> faces_global_ind;
  std::vector<int> cells_global_ind;
  
public:
  XMOF_Reconstructor(const MeshConfig& mesh_config, const IRTolerances& ir_tolerances);
  void set_materials_data(const CellsMatData& mat_data);
  void set_global_ind(const GlobalIndData& global_ind_data);
  
  void construct_minimesh(int icell);
  
  const std::vector<int>& get_cell_materials(int icell) const;
  const std::vector<double>& get_cell_vfracs(int icell) const;
  const std::vector<Point2D>& get_cell_centroids(int icell) const;
  
  int cell_global_ind(int icell) const;
  int face_global_ind(int iface) const;
  int node_global_ind(int inode) const;

  int cell_parent_global_ind(const Cell& cell) const;
  int face_parent_global_ind(const Face& face) const;
  int node_parent_global_ind(const Node& node) const;
  
  const BaseMesh& get_base_mesh() const { return mesh; }
};

}
#endif /* xmof_ir_h */
