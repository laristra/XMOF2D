/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#include <algorithm>
#include "mesh.h"
#include "tools.h"
#include "simple_vector.h"

namespace XMOF2D {

MiniMesh::MiniMesh(const Cell& base_cell_, double ddot_eps, double dist_eps) : base_cell(base_cell_) {
  ddot_eps_ = ddot_eps;
  dist_eps_ = dist_eps;
  nodes_shift = Point2D(0.0, 0.0);
}

const MeshBC& MiniMesh::get_parent_mesh() const {
  return get_base_cell().get_mesh();
}

void MiniMesh::translate_2origin() {
  XMOF2D_ASSERT(nodes_shift == Point2D(0.0, 0.0), "The mesh has been already translated to the origin!");
  
  nodes_shift = BAD_POINT;
  for (int inode = 0; inode < nnodes(); inode++) {
    if (get_node_crd(inode).x < nodes_shift.x)
      nodes_shift.x = get_node_crd(inode).x;
    if (get_node_crd(inode).y < nodes_shift.y)
      nodes_shift.y = get_node_crd(inode).y;
  }
  nodes_shift = -nodes_shift;
  for (int inode = 0; inode < nnodes(); inode++)
    nodes[inode].set_coord(nodes[inode].get_crd() + nodes_shift);
  
  for (int iface = 0; iface < nfaces(); iface++)
    faces[iface].calculate_center();
  for (int icell = 0; icell < ncells(); icell++)
    cells[icell].calculate_center();
}

void MiniMesh::translate_2orig_pos() {
  if (nodes_shift == Point2D(0.0, 0.0))
    return;
  
  nodes_shift = -nodes_shift;
  for (int inode = 0; inode < nnodes(); inode++)
    nodes[inode].set_coord(nodes[inode].get_crd() + nodes_shift);
  
  for (int iface = 0; iface < nfaces(); iface++)
    faces[iface].calculate_center();
  for (int icell = 0; icell < ncells(); icell++)
    cells[icell].calculate_center();
  
  nodes_shift = Point2D(0.0, 0.0);
}

void MiniMesh::Split(int cell_ind, const std::vector<double>& n, double d2orgn, int cell_mat) {
  XMOF2D_ASSERT((cell_ind >= 0) && (cell_ind < ncells()), "Unknown subcell " << cell_ind);
  
  translate_2origin();
  
  Cell& cur_cell = cells[cell_ind];
  Factory factory(*this, faces);
  
  int new_cell_ind[2] = {cell_ind, ncells()};
  int nsides = cur_cell.nfaces();
  std::vector< std::vector<int> > cells_faces(2);
  int int_point_ind[2] = {-1, -1};
  int nearest_lhp_ivrt[2] = {-1, -1};
  int iside = 0, ilastside = nsides - 1;
  while (iside <= ilastside) {
    int icurface = cur_cell.get_face_index(iside);
    Segment cur_side = get_face(icurface).as_segment();
    SegLine::Position side_pos = cur_side.PosWRT2Line(n, d2orgn, dist_eps());
    switch (side_pos) {
      case SegLine::Position::BELOW:
        cells_faces[0].push_back(icurface);
        iside++;
        continue;
        
      case SegLine::Position::ABOVE: {
        cells_faces[1].push_back(icurface);
        bool in_normal = !get_face(icurface).normal_is_out(cell_ind);
        faces[icurface].set_cell_index(in_normal, new_cell_ind[1], true);
        iside++;
        continue;
      }
      case SegLine::Position::TOUCHES: {
        XMOF2D_ASSERT(int_point_ind[1] == -1, "Mesh node was found to be the third intersection point!")
        int iv_on = (cur_side[0].PosWRT2Line(n, d2orgn, dist_eps()) == Point2DLine::Position::ON) ? 0 : 1;
        int iip = (int_point_ind[0] == -1) ? 0 : 1;
        int_point_ind[iip] = get_face(icurface).get_node_index(iv_on);

        int iv_off = (iv_on + 1)%2;
        int ifset = (cur_side[iv_off].PosWRT2Line(n, d2orgn, dist_eps()) == Point2DLine::Position::BELOW) 
                    ? 0 : 1;
        cells_faces[ifset].push_back(icurface);
        if (ifset == 1) {
          bool in_normal = !get_face(icurface).normal_is_out(cell_ind);
          faces[icurface].set_cell_index(in_normal, new_cell_ind[1], true);
        }

        XMOF2D_ASSERT(iside != ilastside, "Face lies on the cutting line!");

        int iadjface = cur_cell.get_face_index(iside + 1);
        if ((int_point_ind[iip] == get_face(iadjface).get_node_index(0)) ||
            (int_point_ind[iip] == get_face(iadjface).get_node_index(1)))
          iside += 2;
        else {
          XMOF2D_ASSERT(iside == 0, "Face lies on the cutting line!");
          iadjface = cur_cell.get_face_index(nsides - 1);
          iside++; ilastside--;
        }
        cells_faces[(ifset + 1)%2].push_back(iadjface);
        if (ifset == 0) {
          bool in_normal = !get_face(iadjface).normal_is_out(cell_ind);
          faces[iadjface].set_cell_index(in_normal, new_cell_ind[1], true);
        }
        continue;
      }
        
      case SegLine::Position::INTERSECTS: {
        Point2D int_point = cur_side.LineIntersect(n, d2orgn, ddot_eps(), dist_eps());
        XMOF2D_ASSERT(int_point_ind[1] == -1, "Extra intersection!");
        int iip = (int_point_ind[0] == -1) ? 0 : 1;

        for (int ivrt = 0; ivrt < 2; ivrt++) {
          if (cur_side[ivrt].PosWRT2Line(n, d2orgn, dist_eps()) == Point2DLine::Position::BELOW) {
            nearest_lhp_ivrt[iip] = get_face(icurface).get_node_index(ivrt);
            break;
          }
        }

        if (nearest_lhp_ivrt[0] == nearest_lhp_ivrt[1]) {
          XMOF2D_ASSERT(std::fabs(vpz(get_node_crd(nearest_lhp_ivrt[1]), int_point, 
                                      get_node_crd(int_point_ind[0]))) >= 2*epsilon,
                        "Area of the subcell below the line is smaller than the machine epsilon!");
        }
        
        int inew_node = nnodes();
        nodes.push_back(Node(int_point, inew_node, *this));
        int_point_ind[iip] = inew_node;
        
        int inewface_cell = -1;
        Point2DLine::Position v0_pos = cur_side[0].PosWRT2Line(n, d2orgn, dist_eps());
        bool in_normal = !get_face(icurface).normal_is_out(cell_ind);
        if (v0_pos == Point2DLine::Position::BELOW) {
          cells_faces[0].push_back(icurface);
          inewface_cell = new_cell_ind[1];
        }
        else {
          cells_faces[1].push_back(icurface);
          inewface_cell = new_cell_ind[0];
        }
        
        int inewface = factory.get_face(inewface_cell,
                                        { inew_node, get_face(icurface).get_node_index(1) });
        faces[icurface].set_node_index(1, inew_node);
        faces[icurface].calculate_size();
        faces[icurface].calculate_center();
        
        if (!get_face(icurface).is_boundary()) {
          faces[inewface].set_cell_index(in_normal, inewface_cell, true);
          int ineigh_cell = get_face(icurface).get_neigh_cell(cell_ind);
          faces[inewface].set_cell_index(!in_normal, ineigh_cell, true);
          cells[ineigh_cell].add_face(inewface, icurface, true);
        }
        else
          faces[inewface].set_parent(cur_cell.get_face(iside).iparent());
        
        if (v0_pos == Point2DLine::Position::ABOVE)
          faces[icurface].set_cell_index(in_normal, new_cell_ind[1], true);
        
        if (v0_pos == Point2DLine::Position::ABOVE)
          cells_faces[0].push_back(inewface);
        else
          cells_faces[1].push_back(inewface);
        iside++; 
        continue;
      }
        
      default:
        THROW_EXCEPTION("Face lies on the cutting line!");
        break;
    }
  }
  XMOF2D_ASSERT((int_point_ind[0] != -1) && (int_point_ind[1] != -1), "Less than two intersection points were found!");
  
  std::vector< std::vector<int> > iendfaces(2);
  for (int ic = 0; ic < 2; ic++) {
    iendfaces[ic] = { -1, -1 };
    if (cells_faces[ic].size() == 2)
      continue;
    for (int icf = 0; icf < cells_faces[ic].size(); icf++) {
      for (int ifn = 0; ifn < 2; ifn++)
        for (int iipt = 0; iipt < 2; iipt++)
          if (get_face(cells_faces[ic][icf]).get_node_index(ifn) == int_point_ind[iipt])
            iendfaces[ic][iipt] = icf;
      if ((iendfaces[ic][0] != -1) && (iendfaces[ic][1] != -1))
        break;
    }
  }
  for (int ic = 0; ic < 2; ic++) {
    int ncfaces = (int) cells_faces[ic].size();
    if (ncfaces == 2) {
      if (!XMOF2D::is_ccw(get_face(cells_faces[ic][0]).as_segment(),
                          get_face(cells_faces[ic][1]).as_segment(), dist_eps(), ddot_eps()))
        std::rotate(cells_faces[ic].begin(), cells_faces[ic].begin() + 1, cells_faces[ic].end());
      
      continue;
    }
    if (((iendfaces[ic][0] == 0) && (iendfaces[ic][1] == ncfaces - 1)) ||
        ((iendfaces[ic][1] == 0) && (iendfaces[ic][0] == ncfaces - 1)))
    continue;
    std::rotate(cells_faces[ic].begin(), cells_faces[ic].begin() + iendfaces[ic][1], cells_faces[ic].end());
  }
  
  std::vector< std::vector<int> > cells_nodes(2);
  for (int ic = 0; ic < 2; ic++) {
    cells_nodes[ic].resize(cells_faces[ic].size() + 1);
    for (int iside = 0; iside < cells_faces[ic].size() - 1; iside++) {
      const Face& cur_face = get_face(cells_faces[ic][iside]);
      cells_nodes[ic][iside] = (cur_face.normal_is_out(new_cell_ind[ic])) ? cur_face.get_node_index(0) : cur_face.get_node_index(1);
    }
    int ilastside = (int) cells_faces[ic].size() - 1;
    const Face& last_face = get_face(cells_faces[ic][ilastside]);
    int rev_dir = (last_face.normal_is_out(new_cell_ind[ic])) ? 0 : 1;
    cells_nodes[ic][ilastside] = last_face.get_node_index(rev_dir);
    cells_nodes[ic][ilastside + 1] = last_face.get_node_index((rev_dir + 1) % 2);
  }
  
  int iint_face = factory.get_face(new_cell_ind[0],
                                   { cells_nodes[0][cells_nodes[0].size() - 1], cells_nodes[0][0] });
  faces[iint_face].set_cell_index(true, new_cell_ind[1], true);
  for (int ic = 0; ic < 2; ic++)
  cells_faces[ic].push_back(iint_face);
  
  translate_2orig_pos();
  
  cells[cell_ind] = Cell(cells_nodes[0], cells_faces[0], cell_ind, *this);
  cells.push_back(Cell(cells_nodes[1], cells_faces[1], ncells(), *this));
  
  cells_material.resize(ncells());
  cells_material[cell_ind] = cell_mat;
}

void MiniMesh::ConstructMinimesh(const std::vector<double>& a2OX, const std::vector<double>& d2orgn, const std::vector<int>& cells_mat) {
  int nc = (int) cells_mat.size();
  
  std::vector<double> cur_line_n(2);
  for (int icell = 0; icell < nc - 1; icell++) {
    cur_line_n[0] = cos(a2OX[icell]);
    cur_line_n[1] = sin(a2OX[icell]);
    double n_scal = 1.0/dnrm2(cur_line_n);
    dscal(n_scal, cur_line_n);
    Split(icell, cur_line_n, d2orgn[icell], cells_mat[icell]);
  }
  cells_material[nc - 1] = cells_mat[nc - 1];
  
  set_nodes_ifaces();
  set_nodes_icells();
}

}
