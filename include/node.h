/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Los Alamos National Security, LLC.
 All rights reserved.
*/

#ifndef node_h
#define node_h

#include <vector>
#include "geometry.h"

namespace XMOF2D {

class Face;
class Cell;
class MeshBC;

class Node {
private:
  Point2D coord;
  
  const MeshBC& mesh;
  int node_ind;
  std::vector<int> faces;
  std::vector<int> cells;
  int parent_ind;
  
public:
  Node(const Point2D& _coord, int _node_ind, const MeshBC& _mesh);
  
  int index() const { return node_ind; }
  int iparent() const { return parent_ind; }
  const Point2D& get_crd() const { return coord; };
  bool has_valid_crd() const;
  
  bool is_boundary() const;
  std::vector<int> iparent_faces() const;
  
  int nfaces() const { return (int) faces.size(); }
  int ncells() const { return (int) cells.size(); }
  const std::vector<int>& get_faces() const { return faces; }
  const std::vector<int>& get_cells() const { return cells; }
  
  const Face& get_face(int i) const;
  const Cell& get_cell(int i) const;
  const MeshBC& get_mesh() const {return mesh; }

  void set_coord(const Point2D& new_coord);
  void set_parent(int parent);
  void set_faces(const std::vector<int>& ifaces) { faces = ifaces; }
  void set_cells(const std::vector<int>& icells) { cells = icells; }
  
  bool faces_in_ccw_order(int iface1, int iface2) const;
  bool cells_in_ccw_order(int icell1, int icell2) const;
  
  friend std::ostream& operator<<(std::ostream& os, const Node& node);
  friend std::ofstream& operator<<(std::ofstream& ofs, const Node& node);
};

}
#endif /* node_h */
