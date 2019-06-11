/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#include "cell.h"
#include "exception.h"
#include "simple_vector.h"
#include <map>
#include <fstream>
#include <numeric>
#include <algorithm>

namespace XMOF2D {

Cell::Cell(int cell_ind_, const MeshBC& mesh_) : mesh(mesh_) {
  cell_ind = cell_ind_;
  size_ = 0.0;
  center_ = BAD_POINT;
  minimesh = NULL;
}

Cell::Cell(const std::vector<int>& faces_, int cell_ind_, const MeshBC& mesh_) : mesh(mesh_) {
  faces = faces_;
  cell_ind = cell_ind_;
  nodes_from_faces_top();
  calculate_size();
  calculate_center();
  minimesh = NULL;
}

Cell::Cell(const std::vector<int>& nodes_, const std::vector<int>& faces_, int cell_ind_, const MeshBC& mesh_) : mesh(mesh_) {
  nodes = nodes_;
  faces = faces_;
  cell_ind = cell_ind_;
  calculate_size();
  calculate_center();
  minimesh = NULL;
}

Cell::Cell(const Cell& c) : mesh(c.mesh), cell_ind(c.cell_ind), parent_ind(c.parent_ind),
nodes(c.nodes), faces(c.faces), size_(c.size_), center_(c.center_) {
  if (c.has_minimesh())
    minimesh = new MiniMesh(c.get_minimesh());
  else
    minimesh = NULL;
}

Cell::~Cell() {
  delete_minimesh();
}

void Cell::delete_minimesh() {
  if (minimesh) {
    delete minimesh;
    minimesh = NULL;
  }
}

const Cell& Cell::operator=(const Cell& c) {
  if (this != &c) {
    this->~Cell();
    new (this) Cell(c);
  }
  
  return *this;
}

int Cell::get_node_index(int i) const {
  XMOF2D_ASSERT_SIZE_LESS(i, nodes.size());
  return nodes[i];
}

const Point2D& Cell::get_node_crd(int i) const {
  XMOF2D_ASSERT_SIZE_LESS(i, nodes.size());
  return mesh.get_node_crd(nodes[i]);
}

const Node& Cell::get_node(int i) const {
  XMOF2D_ASSERT_SIZE_LESS(i, nodes.size());
  return mesh.get_node(nodes[i]);
}

int Cell::get_face_index(int i) const {
  XMOF2D_ASSERT_SIZE_LESS(i, faces.size());
  return faces[i];
}

const Face& Cell::get_face(int i) const {
  XMOF2D_ASSERT_SIZE_LESS(i, faces.size());
  return mesh.get_face(faces[i]);
}

int Cell::get_material_index() const { return mesh.get_cells_material(cell_ind); }

const MiniMesh& Cell::get_minimesh() const {
  XMOF2D_ASSERT(minimesh != NULL, "No minimesh in this cell!");
  return *minimesh;
}

void Cell::add_face(int inew_face, int ineigh_face, bool after) {
  XMOF2D_ASSERT((inew_face < mesh.nfaces()) && (ineigh_face < mesh.nfaces()), "Unknown face indices!");
  std::vector<int>::iterator it = std::find (faces.begin(), faces.end(), ineigh_face);
  XMOF2D_ASSERT(it != faces.end(), "Given face is not a side of this cell!");
  if (after) it++;
  faces.insert(it, inew_face);
  nodes_from_faces_top();
}

SimpleConvex Cell::as_simpleconvex() const {
  std::vector<Point2D> scv(nfaces());
  
  for (int inode = 0; inode < nfaces(); inode++)
    scv[inode] = get_node_crd(inode);
  
  return SimpleConvex(scv, mesh.ddot_eps());
}

void Cell::calculate_size() {
  size_ = as_simpleconvex().size();
}

void Cell::calculate_center() {
  center_ = as_simpleconvex().center();
}

void Cell::nodes_from_faces_geo() {
  nodes.clear();
  
  if (is_ccw(center(), get_face(0).get_node_crd(0), get_face(0).get_node_crd(1), mesh.ddot_eps()))
    nodes.push_back(get_face(0).get_node_index(0));
  else
    nodes.push_back(get_face(0).get_node_index(1));
  
  for (int iside = 0; iside < nfaces() - 1; iside++) {
    if (is_ccw(center(), get_face(iside).get_node_crd(0), get_face(iside).get_node_crd(1), mesh.ddot_eps()))
      nodes.push_back(get_face(iside).get_node_index(1));
    else
      nodes.push_back(get_face(iside).get_node_index(0));
  }
}

void Cell::nodes_from_faces_top() {
  nodes.clear();
  
  if (get_face(0).normal_is_out(cell_ind))
    nodes.push_back(get_face(0).get_node_index(0));
  else
    nodes.push_back(get_face(0).get_node_index(1));
  
  for (int iside = 0; iside < nfaces() - 1; iside++) {
    if (get_face(iside).normal_is_out(cell_ind))
      nodes.push_back(get_face(iside).get_node_index(1));
    else
      nodes.push_back(get_face(iside).get_node_index(0));
  }
}

std::ostream& operator<<(std::ostream& os, const Cell& cell) {
  for (int iface = 0; iface < cell.faces.size(); iface++) {
    os << "Face #" << iface << ": " << cell.faces[iface] << std::endl;
    os << cell.get_face(iface);
  }
  return os;
}

std::ofstream& operator<<(std::ofstream& ofs, const Cell& cell) {
  ofs << std::scientific;
  ofs.precision(16);
  
  ofs << "# Cell index:" << std::endl << cell.index() << std::endl;
  ofs << "# Indices of cell nodes:" << std::endl;
  for (int inode = 0; inode < cell.nfaces(); inode++)
    ofs << cell.get_node_index(inode) << " ";
  ofs << std::endl;
  ofs << "# Indices of cell faces:" << std::endl;
  for (int iface = 0; iface < cell.nfaces(); iface++)
    ofs << cell.get_face_index(iface) << " ";
  ofs << std::endl;
  ofs << "# Index of the parent cell in the parent mesh:" << std::endl;
  ofs << cell.iparent() << std::endl;
  ofs << "# Cell size:" << std::endl;
  ofs << cell.size() << std::endl;
  ofs << "# Cell center:" << std::endl;
  ofs << cell.center().x << "\t" << cell.center().y << std::endl;
  ofs << "# Index of cell material:" << std::endl;
  ofs << cell.get_material_index() << std::endl;
  if (cell.get_material_index() == -1) {
    XMOF2D_ASSERT(cell.has_minimesh(), "Multi-material cell " << cell.index() << " does not have a minimesh!");
    ofs << "# Cell minimesh:" << std::endl;
    ofs << cell.get_minimesh() << std::endl;
  }
  
  return ofs;
}

void Cell::InitializeMinimesh() {
  XMOF2D_ASSERT(minimesh == NULL, "Minimesh is already initialized!");
  
  int nextfaces = nfaces();
  minimesh = new MiniMesh(*this, mesh.ddot_eps(), mesh.dist_eps());
  std::vector<int> first_cell_nodes(nextfaces);
  std::vector<int> first_cell_faces(nextfaces);
  std::iota(first_cell_nodes.begin(), first_cell_nodes.end(), 0);
  std::iota(first_cell_faces.begin(), first_cell_faces.end(), 0);
  for (int inode = 0; inode < nextfaces; inode++) {
    const Node& cnode = get_node(inode);
    minimesh->nodes.push_back(Node(cnode.get_crd(), inode, *minimesh));
    minimesh->nodes[inode].set_parent(get_node_index(inode));
  }
  for (int iface = 0; iface < nextfaces; iface++) {
    minimesh->faces.push_back(Face({ iface, (iface + 1)%nextfaces }, iface, *minimesh));
    minimesh->faces[iface].set_cell_index(false, 0, true);
    minimesh->faces[iface].set_parent(get_face_index(iface));
  }
  
  minimesh->cells.push_back(Cell(first_cell_nodes, first_cell_faces, 0, *minimesh));
}

void Cell::ConstructMinimesh(const std::vector<int>& imats, const std::vector<double>& vfracs,
                             const std::vector<Point2D>& centroids) {
  int nsubcells = (int) imats.size();
  XMOF2D_ASSERT((vfracs.size() == nsubcells) && (centroids.size() == nsubcells), "Incorrect subcell data!");
  
  InitializeMinimesh();
  SimpleConvex cell_SC = as_simpleconvex();
  std::vector<int> opt_mat_order(nsubcells);
  std::vector<double> opt_a2OX(nsubcells - 1);
  std::vector<double> opt_d2orgn(nsubcells - 1);
  
  cell_SC.compute_optimal_cuts(vfracs, centroids, mesh.ang_eps(), mesh.area_eps(),
                               mesh.ddot_eps(), mesh.dist_eps(), mesh.mof_max_iter(),
                               opt_mat_order, opt_a2OX, opt_d2orgn);
  
  std::vector<int> isc_mat(nsubcells);
  for (int isc = 0; isc < nsubcells; isc++)
    isc_mat[isc] = imats[opt_mat_order[isc]];

  minimesh->ConstructMinimesh(opt_a2OX, opt_d2orgn, isc_mat);
  minimesh->set_orig_icells(opt_mat_order);
}

}
