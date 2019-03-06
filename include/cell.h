/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#ifndef cell_h
#define cell_h

#include <vector>
#include "geometry.h"
#include "mesh.h"

namespace XMOF2D {
class MiniMesh;

class Cell {
protected:
  const MeshBC& mesh;
  int cell_ind;
  int parent_ind;
  
  std::vector<int> nodes;
  std::vector<int> faces;
  
  MiniMesh* minimesh;
  
  double size_;
  Point2D center_;
  
  void nodes_from_faces_geo();
  void nodes_from_faces_top();
  void add_face(int inew_face, int ineigh_face, bool after = true);
  void InitializeMinimesh();
  
public:
  Cell(int cell_ind_, const MeshBC& mesh_);
  Cell(const std::vector<int>& faces_, int cell_ind_, const MeshBC& mesh_);
  Cell(const std::vector<int>& nodes_, const std::vector<int>& faces_, int cell_ind_, const MeshBC& mesh_);
  Cell(const Cell& c);
  ~Cell();
  
  const Cell& operator=(const Cell& c);
  
  void calculate_size();
  void calculate_center();
  
  double size() const { return size_; }
  Point2D center() const { return center_; }
  
  int index() const { return cell_ind; }
  int iparent() const { return parent_ind; }
  int nfaces() const { return (int) faces.size(); }
  const std::vector<int>& get_nodes() const { return nodes; }
  const std::vector<int>& get_faces() const { return faces; }
  const MeshBC& get_mesh() const { return mesh; }

  bool has_minimesh() const { return minimesh != NULL; }
  void delete_minimesh();
  
  int get_node_index(int i) const;
  const Point2D& get_node_crd(int i) const;
  const Node& get_node(int i) const;
  int get_face_index(int i) const;
  const Face& get_face(int i) const;
  int get_material_index() const;
  const MiniMesh& get_minimesh() const;

  SimpleConvex as_simpleconvex() const;
  
  void ConstructMinimesh(const std::vector<int>& imats, const std::vector<double>& vfracs,
                         const std::vector<Point2D>& centroids);
  
  friend std::ostream& operator<<(std::ostream& os, const XMOF2D::Cell& cell);
  friend std::ofstream& operator<<(std::ofstream& ofs, const XMOF2D::Cell& cell);
  friend class MiniMesh;
};

}
#endif /* cell_h */
