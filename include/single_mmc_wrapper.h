/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Los Alamos National Security, LLC.
 All rights reserved.
*/

#ifndef single_mmc_wrapper_h
#define single_mmc_wrapper_h

#include "xmof2D.h"
#include "mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

void xmof2d_set_distance_tolerance(double dist_eps);
void xmof2d_get_distance_tolerance(double* dist_eps);
void xmof2d_set_intersect_division_tolerance(double div_eps);
void xmof2d_get_intersect_division_tolerance(double* div_eps);
void xmof2d_set_dot_product_tolerance(double ddot_eps);
void xmof2d_get_dot_product_tolerance(double* ddot_eps);
void xmof2d_set_vol_fraction_tolerance(double vfrac_eps);
void xmof2d_get_vol_fraction_tolerance(double* vfrac_eps);
void xmof2d_set_angle_tolerance(double ang_eps);
void xmof2d_get_angle_tolerance(double* ang_eps);
void xmof2d_set_max_iter_num(int mof_max_iter);
void xmof2d_get_max_iter_num(int* max_iter);
void xmof2d_set_mmc_vertices(int nvrts, double* crds);
void xmof2d_set_materials_data(int nmaterials, int* mat_ids, double* vol_fractions, double* centroids);
void xmof2d_initialize_reconstructor();
void xmof2d_free_reconstructor();
void xmof2d_perform_reconstruction();
void xmof2d_mesh_get_ncells(int* ncells);
void xmof2d_mesh_get_nfaces(int* nfaces);
void xmof2d_mesh_get_nnodes(int* nnodes);
void xmof2d_cell_get_mat_id(int icell, int* mat_id);
void xmof2d_cell_get_nfaces(int icell, int* nfaces);
void xmof2d_cell_get_face_ids(int icell, int* ifaces);
void xmof2d_cell_get_node_ids(int icell, int* inodes);
void xmof2d_cell_get_size(int icell, double* size);
void xmof2d_cell_get_center(int icell, double* center_coords);
void xmof2d_face_get_cell_ids(int iface, int* left_cell_id, int* right_cell_id);
void xmof2d_face_get_node_ids(int iface, int* tail_id, int* head_id);
void xmof2d_face_is_boundary(int iface, bool* is_boundary);
void xmof2d_face_get_parent_face_id(int iface, int* iparent);
void xmof2d_face_get_size(int iface, double* size);
void xmof2d_face_get_center(int iface, double* center_coords);
void xmof2d_node_get_coords(int inode, double* coords);
void xmof2d_node_get_ncells(int inode, int* ncells);
void xmof2d_node_get_cell_ids(int inode, int* icells);
void xmof2d_node_get_nfaces(int inode, int* nfaces);
void xmof2d_node_get_face_ids(int inode, int* ifaces);
void xmof2d_node_is_boundary(int inode, bool* is_boundary);
void xmof2d_node_get_parent_node_id(int inode, int* iparent);
void xmof2d_node_get_nparent_faces(int inode, int* nparent_faces);
void xmof2d_node_get_parent_face_ids(int inode, int* iparent_faces);

#ifdef __cplusplus
}
#endif

#endif
