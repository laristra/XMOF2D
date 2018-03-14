/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Los Alamos National Security, LLC.
 All rights reserved.
*/

#ifndef single_mmc_wrapper_cpp
#define single_mmc_wrapper_cpp

#include "single_mmc_wrapper.h"

XMOF2D::XMOF_Reconstructor* xmof2d_rec = nullptr;
XMOF2D::IRTolerances xmof2d_tol = {
  .dist_eps = 1.0e-15,
  .div_eps = 1.0e-6,
  .ddot_eps = 1.0e-14,
  .vfrac_eps = 1.0e-14,
  .ang_eps = 1.0e-13,
  .mof_max_iter = 10000
};
XMOF2D::MeshConfig xmof2d_single_mmc_mconfig;
XMOF2D::CellsMatData xmof2d_mmc_mat_data;

void xmof2d_set_distance_tolerance(double dist_eps) {
  if (dist_eps <= 0.0) {
    throw XMOF2D::Exception("Distance tolerance should be positive!");
  }
  xmof2d_tol.dist_eps = dist_eps;
}

void xmof2d_get_distance_tolerance(double* dist_eps) {
  *dist_eps = xmof2d_tol.dist_eps;
}

void xmof2d_set_intersect_division_tolerance(double div_eps) {
  if (div_eps <= 0.0) {
    throw XMOF2D::Exception("Division tolerance should be positive!");
  }
  xmof2d_tol.div_eps = div_eps;
}

void xmof2d_get_intersect_division_tolerance(double* div_eps) {
  *div_eps =  xmof2d_tol.div_eps;
}

void xmof2d_set_dot_product_tolerance(double ddot_eps) {
  if (ddot_eps <= 0.0) {
    throw XMOF2D::Exception("Dot product tolerance should be positive!");
  }
  xmof2d_tol.ddot_eps = ddot_eps;
}

void xmof2d_get_dot_product_tolerance(double* ddot_eps) {
  *ddot_eps = xmof2d_tol.ddot_eps;
}

void xmof2d_set_vol_fraction_tolerance(double vfrac_eps) {
  if (vfrac_eps <= 0.0) {
    throw XMOF2D::Exception("Volume fraction tolerance should be positive!");
  }
  xmof2d_tol.vfrac_eps = vfrac_eps;
}

void xmof2d_get_vol_fraction_tolerance(double* vfrac_eps) {
  *vfrac_eps = xmof2d_tol.vfrac_eps;
}

void xmof2d_set_angle_tolerance(double ang_eps) {
  if (ang_eps <= 0.0) {
    throw XMOF2D::Exception("Angle tolerance should be positive!");
  }
  xmof2d_tol.ang_eps = ang_eps;
}

void xmof2d_get_angle_tolerance(double* ang_eps) {
  *ang_eps = xmof2d_tol.ang_eps;
}

void xmof2d_set_max_iter_num(int mof_max_iter) {
  if (mof_max_iter < 1) {
    throw XMOF2D::Exception("Max number of iterations should be positive!");
  }
  xmof2d_tol.mof_max_iter = mof_max_iter;
}

void xmof2d_get_max_iter_num(int* max_iter) {
  *max_iter = xmof2d_tol.mof_max_iter;
}

void xmof2d_set_mmc_vertices(int nvrts, double* crds) {
  if (nvrts < 3) {
    throw XMOF2D::Exception("A multi-material cell should have at least three vertices!");
  }
  xmof2d_single_mmc_mconfig.nodes_coords.resize(nvrts);
  xmof2d_single_mmc_mconfig.ifaces_nodes.resize(nvrts);
  xmof2d_single_mmc_mconfig.icells_faces.resize(1);
  xmof2d_single_mmc_mconfig.icells_faces[0].resize(nvrts);
  for (int ivrt = 0; ivrt < nvrts; ivrt++) {
    xmof2d_single_mmc_mconfig.nodes_coords[ivrt].x = crds[2*ivrt];
    xmof2d_single_mmc_mconfig.nodes_coords[ivrt].y = crds[2*ivrt + 1];
    xmof2d_single_mmc_mconfig.ifaces_nodes[ivrt] = { ivrt, (ivrt + 1)%nvrts };
    xmof2d_single_mmc_mconfig.icells_faces[0][ivrt] = ivrt;
  }
  xmof2d_single_mmc_mconfig.cells_material.clear();
  xmof2d_single_mmc_mconfig.cells_material.resize(1, -1);
}

void xmof2d_set_materials_data(int nmaterials, int* mat_ids, double* vol_fractions, double* centroids) {
  if (nmaterials < 2) {
    throw XMOF2D::Exception("A multi-material cell should contain at least two materials!");
  }
  xmof2d_mmc_mat_data.cells_materials.resize(1);
  xmof2d_mmc_mat_data.cells_vfracs.resize(1);
  xmof2d_mmc_mat_data.cells_centroids.resize(1);
  xmof2d_mmc_mat_data.cells_materials[0].resize(nmaterials);
  xmof2d_mmc_mat_data.cells_vfracs[0].resize(nmaterials);
  xmof2d_mmc_mat_data.cells_centroids[0].resize(nmaterials);
  for (int imat = 0; imat < nmaterials; imat++) {
    xmof2d_mmc_mat_data.cells_materials[0][imat] = mat_ids[imat];
    xmof2d_mmc_mat_data.cells_vfracs[0][imat] = vol_fractions[imat];
    xmof2d_mmc_mat_data.cells_centroids[0][imat].x = centroids[2*imat];
    xmof2d_mmc_mat_data.cells_centroids[0][imat].y = centroids[2*imat + 1];
  }
  if (xmof2d_rec != nullptr)
    xmof2d_rec->set_materials_data(xmof2d_mmc_mat_data);
}

void xmof2d_initialize_reconstructor() {
  if (xmof2d_single_mmc_mconfig.icells_faces.empty()) {
    throw XMOF2D::Exception("Vertices of a multi-material cell should be set first!");
  }
  
  delete xmof2d_rec;
  xmof2d_rec = new XMOF2D::XMOF_Reconstructor(xmof2d_single_mmc_mconfig, xmof2d_tol);
  if (!xmof2d_mmc_mat_data.cells_materials.empty())
    xmof2d_rec->set_materials_data(xmof2d_mmc_mat_data);
}

void xmof2d_free_reconstructor() {
  delete xmof2d_rec;
  xmof2d_rec = nullptr;
}

void xmof2d_perform_reconstruction() {
  if (xmof2d_rec->get_cell_materials(0).empty()) {
    throw XMOF2D::Exception("Materials data should be set first!");
  }
  xmof2d_rec->construct_minimesh(0);
}

void xmof2d_mesh_get_ncells(int* ncells) {
  *ncells = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().ncells();
}

void xmof2d_mesh_get_nfaces(int* nfaces) {
  *nfaces = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().nfaces();
}

void xmof2d_mesh_get_nnodes(int* nnodes) {
  *nnodes = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().nnodes();
}

void xmof2d_cell_get_mat_id(int icell, int* mat_id) {
  *mat_id = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_cell(icell).get_material_index();
}

void xmof2d_cell_get_nfaces(int icell, int* nfaces) {
  *nfaces = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_cell(icell).nfaces();
}

void xmof2d_cell_get_face_ids(int icell, int* ifaces) {
  const std::vector<int>& cell_ifaces = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_cell(icell).get_faces();
  std::copy(cell_ifaces.begin(), cell_ifaces.end(), ifaces);
}

void xmof2d_cell_get_node_ids(int icell, int* inodes) {
  const std::vector<int>& cell_inodes = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_cell(icell).get_nodes();
  std::copy(cell_inodes.begin(), cell_inodes.end(), inodes);
}

void xmof2d_cell_get_size(int icell, double* size) {
  *size = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_cell(icell).size();
}

void xmof2d_cell_get_center(int icell, double* center_coords) {
  XMOF2D::Point2D cen_crd = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_cell(icell).center();
  center_coords[0] = cen_crd.x;
  center_coords[1] = cen_crd.y;
}

void xmof2d_face_get_cell_ids(int iface, int* left_cell_id, int* right_cell_id) {
  *left_cell_id = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_face(iface).get_cell_index(false);
  *right_cell_id = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_face(iface).get_cell_index(true);
}

void xmof2d_face_get_node_ids(int iface, int* tail_id, int* head_id) {
  const std::vector<int>& face_inodes = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_face(iface).get_nodes();
  *tail_id = face_inodes[0];
  *head_id = face_inodes[1];
}

void xmof2d_face_is_boundary(int iface, bool* is_boundary) {
  *is_boundary = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_face(iface).is_boundary();
}

void xmof2d_face_get_parent_face_id(int iface, int* iparent) {
  *iparent = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_face(iface).iparent();
}

void xmof2d_face_get_size(int iface, double* size) {
  *size = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_face(iface).size();
}

void xmof2d_face_get_center(int iface, double* center_coords) {
  XMOF2D::Point2D cen_crd = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_face(iface).center();
  center_coords[0] = cen_crd.x;
  center_coords[1] = cen_crd.y;
}

void xmof2d_node_get_coords(int inode, double* coords) {
  const XMOF2D::Point2D& node_crd = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_node(inode).get_crd();
  coords[0] = node_crd.x;
  coords[1] = node_crd.y;
}

void xmof2d_node_get_ncells(int inode, int* ncells) {
  *ncells = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_node(inode).ncells();
}

void xmof2d_node_get_cell_ids(int inode, int* icells) {
  const std::vector<int>& node_icells = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_node(inode).get_cells();
  std::copy(node_icells.begin(), node_icells.end(), icells);
}

void xmof2d_node_get_nfaces(int inode, int* nfaces) {
  *nfaces = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_node(inode).nfaces();
}

void xmof2d_node_get_face_ids(int inode, int* ifaces) {
  const std::vector<int>& node_ifaces = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_node(inode).get_faces();
  std::copy(node_ifaces.begin(), node_ifaces.end(), ifaces);
}

void xmof2d_node_is_boundary(int inode, bool* is_boundary) {
  *is_boundary = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_node(inode).is_boundary();
}

void xmof2d_node_get_parent_node_id(int inode, int* iparent) {
  *iparent = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_node(inode).iparent();
}

void xmof2d_node_get_nparent_faces(int inode, int* nparent_faces) {
  *nparent_faces = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_node(inode).iparent_faces().size();
}

void xmof2d_node_get_parent_face_ids(int inode, int* iparent_faces) {
  std::vector<int> node_iparent_faces = xmof2d_rec->
    get_base_mesh().get_cell(0).get_minimesh().get_node(inode).iparent_faces();
  std::copy(node_iparent_faces.begin(), node_iparent_faces.end(), iparent_faces);
}

#endif
