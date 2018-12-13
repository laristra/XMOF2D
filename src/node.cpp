/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Los Alamos National Security, LLC.
 All rights reserved.
*/

#include "node.h"
#include "mesh.h"
#include "exception.h"
#include "tools.h"
#include <fstream>
#include <cmath>

namespace XMOF2D {

Node::Node(const Point2D& _coord, int _node_ind, const MeshBC& _mesh) : mesh(_mesh), node_ind(_node_ind) {
  XMOF2D_ASSERT(_coord != BAD_POINT, "Invalid node coordinate!");
  coord = _coord;
  parent_ind = -1;
}

const Face& Node::get_face(int i) const {
  XMOF2D_ASSERT_SIZE_LESS(i, nfaces());
  return mesh.get_face(faces[i]);
}

const Cell& Node::get_cell(int i) const {
  XMOF2D_ASSERT_SIZE_LESS(i, ncells());
  return mesh.get_cell(cells[i]);
}

bool Node::has_valid_crd() const {
  return (coord != BAD_POINT) && isfinite(coord.x) && isfinite(coord.y);
}

void Node::set_coord(const Point2D& new_coord) {
  XMOF2D_ASSERT(new_coord != BAD_POINT, "Invalid node coordinate!");
  coord = new_coord;
}

void Node::set_parent(int parent) {
  XMOF2D_ASSERT(parent >= 0, "Negative parent node index: " << parent);
  parent_ind = parent;
}

bool Node::is_boundary() const {
  XMOF2D_ASSERT(nfaces() > 0, "Node connectivity information is not established!");
  for (int inf = 0; inf < nfaces(); inf++)
    if (get_face(inf).is_boundary())
      return true;
  
  return false;
}

std::vector<int> Node::iparent_faces() const {
  XMOF2D_ASSERT((ncells() > 0) && (nfaces() > 0),
         "Node connectivity information is not established!");
  std::vector<int> ipfaces;
  if (iparent() != -1) {
    const MiniMesh& minimesh = dynamic_cast<const MiniMesh&>(get_mesh());
    ipfaces = minimesh.get_parent_mesh().get_node(iparent()).get_faces();
  }
  else
    for (int inf = 0; inf < nfaces(); inf++) {
      const Face& cur_face = get_face(inf);
      if (cur_face.is_boundary()) {
        ipfaces.push_back(cur_face.iparent());
        break;
      }
    }
  return ipfaces;
}
  
bool Node::faces_in_ccw_order(int iface1, int iface2) const {
  double polar_angle1, polar_angle2;
  polar_angle1 = atan2(
    mesh.get_node(mesh.get_face(iface1).get_other_node(index())).get_crd().y - coord.y,
    mesh.get_node(mesh.get_face(iface1).get_other_node(index())).get_crd().x - coord.x);
  polar_angle2 = atan2(
    mesh.get_node(mesh.get_face(iface2).get_other_node(index())).get_crd().y - coord.y,
    mesh.get_node(mesh.get_face(iface2).get_other_node(index())).get_crd().x - coord.x);

  return polar_angle1 < polar_angle2;
}

bool Node::cells_in_ccw_order(int icell1, int icell2) const {
  double polar_angle1, polar_angle2;
  polar_angle1 = atan2(
    mesh.get_cell(icell1).center().y - coord.y,
    mesh.get_cell(icell1).center().x - coord.x);
  polar_angle2 = atan2(
    mesh.get_cell(icell2).center().y - coord.y,
    mesh.get_cell(icell2).center().x - coord.x);
  
  return polar_angle1 < polar_angle2;
}
  
std::ostream& operator<<(std::ostream& os, const Node& node) {
  os << "Node coordinate: " << node.coord << std::endl;
  if (node.nfaces() != 0) {
    os << "Node faces: ";
    for (int iface = 0; iface < node.nfaces(); iface++)
      os << node.faces[iface] << " ";
    os << std::endl;
  }
  if (node.ncells() != 0) {
    os << "Node cells: ";
    for (int icell = 0; icell < node.ncells(); icell++)
      os << node.cells[icell] << " ";
    os << std::endl;
  }
  
  return os;
}

std::ofstream& operator<<(std::ofstream& ofs, const Node& node) {
  ofs << std::scientific;
  ofs.precision(16);
  
  ofs << "# Node index:" << std::endl << node.index() << std::endl;
  ofs << "# Node coordinate:" << std::endl << node.coord.x << "\t" << node.coord.y << std::endl;
  if (node.nfaces() != 0) {
    ofs << "Indices of node faces: ";
    for (int iface = 0; iface < node.nfaces(); iface++)
      ofs << node.faces[iface] << " ";
    ofs << std::endl;
  }
  if (node.ncells() != 0) {
    ofs << "Indices of node cells: ";
    for (int icell = 0; icell < node.ncells(); icell++)
      ofs << node.cells[icell] << " ";
    ofs << std::endl;
  }
  ofs << "# Index of the parent node in the parent mesh:" << std::endl;
  ofs << node.iparent() << std::endl;
  
  return ofs;
}

}
