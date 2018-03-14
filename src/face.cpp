/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Los Alamos National Security, LLC.
 All rights reserved.
*/

#include "face.h"
#include "mesh.h"
#include "exception.h"
#include "tools.h"
#include "simple_vector.h"
#include <fstream>
#include <algorithm>

namespace XMOF2D {

Face::Face(const std::vector<int>& _nodes, int _face_ind, const MeshBC& _mesh) : mesh(_mesh) {
  XMOF2D_ASSERT_SIZE(_nodes.size(), 2);
  face_ind = _face_ind;
  nodes = _nodes;
  cell0 = cell1 = -1;
  parent_ind = -1;
  
  calculate_center();
  calculate_size();
  XMOF2D_ASSERT(is_not_equal(size_, 0.0), "Trying to create a zero length face!");
}

bool Face::is_boundary() const { return cell1 == -1; }

std::vector<double> Face::normal() const {
  return as_segment().normal();
}

bool Face::check_ndir() const {
  XMOF2D_ASSERT(cell0 != -1, "Cell for which the normal is outward is not set!");
  std::vector<double> ccen2fcen = center() - get_cell(false).center();
  
  return (ddot(ccen2fcen, normal()) > 0.0);
}

void Face::swap_nodes_order() {
  std::reverse(nodes.begin(), nodes.end());
}

int Face::get_other_node(int inode) const {
  XMOF2D_ASSERT((inode == nodes[0]) || (inode == nodes[1]), "Face does not contain the specified node!")
  return inode == nodes[0] ? nodes[1] : nodes[0];
}
  
int Face::get_neigh_cell(int icell) const {
  XMOF2D_ASSERT((icell == cell0) || (icell == cell1), "Face does not belong to the specified cell!")
  return icell == cell0 ? cell1 : cell0;
}

bool Face::normal_is_out(int icell) const {
  XMOF2D_ASSERT((icell == cell0) || (icell == cell1), "Face does not belong to the specified cell!")
  return icell == cell0;
}

int Face::get_node_index(int i) const {
  XMOF2D_ASSERT_SIZE_LESS(i, nnodes());
  return nodes[i];
}

int Face::get_cell_index(bool in_normal) const {
  return in_normal ? cell1 : cell0;
}

const Cell& Face::get_cell(bool in_normal) const {
  return mesh.get_cell(get_cell_index(in_normal));
}

const Point2D& Face::get_node_crd(int i) const {
  XMOF2D_ASSERT_SIZE_LESS(i, nnodes());
  return mesh.get_node_crd(nodes[i]);
}

void Face::set_parent(int parent) {
  XMOF2D_ASSERT(parent >= 0, "Negative parent face index: " << parent);
  parent_ind = parent;
}

void Face::set_node_index(int i4face, int i4mesh) {
  XMOF2D_ASSERT_SIZE_LESS(i4face, nnodes());
  XMOF2D_ASSERT_SIZE_LESS(i4mesh, mesh.nnodes());
  nodes[i4face] = i4mesh;
}

void Face::set_cell_index(bool in_normal, int icell, bool bypass_ndir_check) {
  if (in_normal)
    cell1 = icell;
  else {
    cell0 = icell;
    if (!bypass_ndir_check && !check_ndir()) swap_nodes_order();
  }
}

Segment Face::as_segment() const {
  return Segment(get_node_crd(0), get_node_crd(1));
}

void Face::calculate_size() {
  size_ = as_segment().size();
}

void Face::calculate_center() {
  center_ = 0.5*(get_node_crd(0) + get_node_crd(1));
}

std::ostream& operator<<(std::ostream& os, const Face& face) {
  os << "Face nodes: ";
  for (int iface = 0; iface < face.nnodes(); iface++)
    os << face.nodes[iface] << " ";
  os << std::endl;
  os << "Face cells(0, 1): " << face.cell0 << ", " << face.cell1 << std::endl;
  
  return os;
}

std::ofstream& operator<<(std::ofstream& ofs, const Face& face) {
  ofs << std::scientific;
  ofs.precision(16);
  
  ofs << "# Face index:" << std::endl << face.index() << std::endl;
  ofs << "# Indices of face nodes:" << std::endl;
  for (int inode = 0; inode < face.nnodes(); inode++)
    ofs << face.get_node_index(inode) << " ";
  ofs << std::endl;
  ofs << "# Indices of face cells:" << std::endl;
    ofs << face.get_cell_index(false) << " " << face.get_cell_index(true) << std::endl;
  ofs << "# Index of the parent face in the parent mesh:" << std::endl;
  ofs << face.iparent() << std::endl;
  ofs << "# Face size:" << std::endl;
  ofs << face.size() << std::endl;
  ofs << "# Face center:" << std::endl;
  ofs << face.center().x << "\t" << face.center().y << std::endl;
  
  return ofs;
}

}
