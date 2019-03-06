/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#ifndef face_h
#define face_h

#include <vector>
#include "geometry.h"

namespace XMOF2D {

class Cell;
class MeshBC;

class Face {
private:
  const MeshBC& mesh;
  int face_ind;
  std::vector<int> nodes;
  int cell0, cell1;
  int parent_ind;
  
  double size_;
  Point2D  center_;
  
public:
  Face(const std::vector<int>& _nodes, int _face_ind, const MeshBC& _mesh);
  
  void calculate_size();
  void calculate_center();
  
  double size() const { return size_; }
  Point2D center() const { return center_; }
  
  int index() const { return face_ind; }
  int iparent() const { return parent_ind; }
  int nnodes() const { return (int) nodes.size(); }
  const std::vector<int>& get_nodes() const { return nodes; }
  
  std::vector<double> normal() const;
  int get_neigh_cell(int icell) const;
  bool normal_is_out(int icell) const;
  bool is_boundary() const;
  
  bool check_ndir() const;
  void swap_nodes_order();
  
  int get_node_index(int i) const;
  const Point2D& get_node_crd(int i) const;
  int get_cell_index(bool in_normal) const;
  const Cell& get_cell(bool in_normal) const;
  
  int get_other_node(int inode) const;
  
  void set_parent(int parent);
  void set_node_index(int i4face, int i4mesh);
  void set_cell_index(bool in_normal, int icell, bool bypass_ndir_check = false);

  Segment as_segment() const;
  
  friend class Factory;
  friend std::ostream& operator<<(std::ostream& os, const Face& face);
  friend std::ofstream& operator<<(std::ofstream& ofs, const Face& face);
};

}
#endif /* face_h */
