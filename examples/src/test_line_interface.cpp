/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "xmof2D.h"
#include "simple_vector.h"

double mat_int_slope = 1.5;
double mat_int_shift = -0.25;
double deps = 1.0e-14;

XMOF2D::MeshConfig make_rectangular_mesh_data(int nx, int ny,
                                              const std::vector<double> xbnd,
                                              const std::vector<double> ybnd);
XMOF2D::CellsMatData read_mat_data(const std::string& matdata_fname);
int PosWRTLine(const XMOF2D::Point2D& p, double line_slope, double line_shift, double eps);

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << std::endl <<
      "Correct usage: test_line_interface <material_data_file> <nx> <ny>" << 
      std::endl << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string matdata_fname = argv[1];
  int nx = atoi(argv[2]);
  int ny = atoi(argv[3]);
  XMOF2D::MeshConfig mesh_cfg = make_rectangular_mesh_data(nx, ny, {0.0, 1.0}, {0.0, 1.0});
  XMOF2D::CellsMatData mat_data = read_mat_data(matdata_fname);
  if (nx*ny != mesh_cfg.icells_faces.size()) {
    std::ostringstream os;
    os << std::endl <<
      "Material data is provided for " << mesh_cfg.icells_faces.size() <<
      " cells, mesh has " << nx*ny << " cells!" << std::endl;
    throw XMOF2D::Exception(os.str());
  }

  XMOF2D::IRTolerances ir_tolerances;
  ir_tolerances.dist_eps = 1.0e-14;
  ir_tolerances.div_eps = 1.0e-6;
  ir_tolerances.ddot_eps = 1.0e-14;
  ir_tolerances.area_eps = 1.0e-14;
  ir_tolerances.ang_eps = 1.0e-13;
  ir_tolerances.mof_max_iter = 10000;
  
  mesh_cfg.cells_material.clear();
  mesh_cfg.cells_material.resize(mesh_cfg.icells_faces.size(), -1);
  
  XMOF2D::XMOF_Reconstructor xmof_ir(mesh_cfg, ir_tolerances);
  xmof_ir.set_materials_data(mat_data);
  
  const XMOF2D::BaseMesh& base_mesh = xmof_ir.get_base_mesh();
  std::vector<int> ncells_with_nmats(1, 0);
  for (int icell = 0; icell < base_mesh.ncells(); icell++) {
    std::cout << "\rProcessing cell " << icell << "/" << base_mesh.ncells() << "...     ";
    xmof_ir.construct_minimesh(icell);
    
    if (base_mesh.get_cell(icell).has_minimesh()) {
      int nmats = base_mesh.get_cell(icell).get_minimesh().ncells();
      if (nmats > ncells_with_nmats.size())
        ncells_with_nmats.resize(nmats, 0);
      ncells_with_nmats[nmats - 1]++;
    }
    else
      ncells_with_nmats[0]++;
  }
  std::cout << std::endl;
  
  std::cout << "Total number of cells: " << nx*ny << std::endl;
  std::cout << "Number of single-material cells: " << ncells_with_nmats[0] << std::endl;
  for (int inm = 0; inm < ncells_with_nmats.size() - 1; inm++)
    std::cout << "Number of cells with " << inm + 2 << " materials: " <<
      ncells_with_nmats[inm + 1] << std::endl;
  std::cout << std::endl;
  

  std::vector<XMOF2D::Point2D> ref_line = {
    XMOF2D::Point2D(0.0, mat_int_shift),
    XMOF2D::Point2D(1.0, mat_int_slope + mat_int_shift)};

  double max_hausdorff = 0.0;
  int ncmp_int = 0;
  for (int icell = 0; icell < base_mesh.ncells(); icell++) {
    const XMOF2D::Cell& cur_cell = base_mesh.get_cell(icell);
    int nsides = cur_cell.nfaces();
    bool should_be_mmc = false;
    std::vector<XMOF2D::Point2D> ref_int_pts;
    for (int iside = 0; iside < nsides; iside++) {
      std::vector<XMOF2D::Point2D> side_vrts = cur_cell.get_face(iside).as_segment().vrts();
      std::vector<int> vrts_pos(2);
      for (int ivrt = 0; ivrt < 2; ivrt++)
        vrts_pos[ivrt] = PosWRTLine(side_vrts[ivrt], mat_int_slope, mat_int_shift, deps);
      
      if (!vrts_pos[0] || !vrts_pos[1]) {
        XMOF2D::Point2D ref_int = vrts_pos[0] ? side_vrts[1] : side_vrts[0];
        bool new_int_pt = true;
        for (int iip = 0; iip < ref_int_pts.size(); iip++)
          if (XMOF2D::distance(ref_int, ref_int_pts[iip]) < deps) {
            new_int_pt = false;
            break;
          }
        if (new_int_pt) {
          ref_int_pts.push_back(ref_int);
          continue;
        }
      }
      else if (vrts_pos[0] != vrts_pos[1]) {
        should_be_mmc = true;
        
        XMOF2D::Point2D ref_int = LineLineIntersect(ref_line, side_vrts);
        ref_int_pts.push_back(ref_int);
      }
    }
    
    if (should_be_mmc) {
      if (ref_int_pts.size() != 2) {
        std::ostringstream os;
        os << std::endl <<
          "Interfce line should have two intersections with a multi-material cell!" <<
          std::endl;
        throw XMOF2D::Exception(os.str());
      }

      const XMOF2D::MiniMesh& cur_mmesh = cur_cell.get_minimesh();
      std::vector<int> iintfaces;
      for (int iface = 0; iface < cur_mmesh.nfaces(); iface++)
        if (!cur_mmesh.get_face(iface).is_boundary())
          iintfaces.push_back(iface);
      if (iintfaces.size() != 1)
        std::cout << "Cell #" << icell <<
          " should have one material interface, but contains " << iintfaces.size() <<
          std::endl;
      else {
        ncmp_int++;
        std::vector<XMOF2D::Segment> cmp_segs = {
          XMOF2D::Segment(ref_int_pts),
          cur_mmesh.get_face(iintfaces[0]).as_segment() };

        for (int iseg = 0; iseg < 2; iseg++)
          for (int ivrt = 0; ivrt < 2; ivrt++) {
            double cur_dist = cmp_segs[iseg].dist(cmp_segs[(iseg + 1)%2][ivrt]);
            if (cur_dist > max_hausdorff)
              max_hausdorff = cur_dist;
          }
      }
    }
    else {
      if (cur_cell.has_minimesh())
        std::cout << "Cell #" << icell << " should be single-material, but contains " <<
          cur_cell.get_minimesh().ncells() << " materials!" << std::endl;
    }
  }
  std::cout << "Checked Hausdorff distance for " << ncmp_int <<
    " in-cell material interfaces" << std::endl;
  std::cout << "Max Hausdorff distance between actual and reference " <<
    "material interface segments -> " << max_hausdorff << std::endl;
  
  return 0;
}

XMOF2D::MeshConfig make_rectangular_mesh_data(int nx, int ny,
                                              const std::vector<double> xbnd,
                                              const std::vector<double> ybnd) {
  XMOF2D::MeshConfig mesh_cfg;
  int nnodes = (nx + 1)*(ny + 1);
  double Lx = xbnd[1] - xbnd[0], Ly = ybnd[1] - ybnd[0];
  double hx = Lx/nx, hy = Ly/ny;
  mesh_cfg.nodes_coords.resize(nnodes);
  
  for (int iy = 0; iy <= ny; iy++)
    for (int ix = 0; ix <= nx; ix++) {
      int ind = iy*(nx + 1) + ix;
      mesh_cfg.nodes_coords[ind].x = xbnd[0] + ix*hx;
      mesh_cfg.nodes_coords[ind].y = ybnd[0] + iy*hy;
    }
  
  for (int ix = 0; ix < nx; ix++) 
    mesh_cfg.ifaces_nodes.push_back({ix, ix + 1});
  
  for (int iy = 0; iy < ny; iy++) {
    int offsetn = (iy + 1)*(nx + 1);
    int offsetc = iy*nx;
    for (int ix = 0; ix < nx; ix++)
      mesh_cfg.ifaces_nodes.push_back({offsetn + ix + 1, offsetn + ix});
  }
  for (int iy = 0; iy < ny; iy++) 
    mesh_cfg.ifaces_nodes.push_back({(iy + 1)*(nx + 1), iy*(nx + 1)});
  
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++) {
      int offsetn = iy*(nx + 1) + ix + 1;
      mesh_cfg.ifaces_nodes.push_back({offsetn, offsetn + nx + 1});
    }
  
  mesh_cfg.icells_faces.resize(nx*ny);
  int vfaces_offset = (ny + 1)*nx;
  for (int iy = 0; iy < ny; iy++)
    for (int ix = 0; ix < nx; ix++)
      mesh_cfg.icells_faces[iy*nx + ix] = {
        iy*nx + ix,
        vfaces_offset + (ix + 1)*ny + iy,
        (iy + 1)*nx + ix,
        vfaces_offset + ix*ny + iy
      };
  
  return mesh_cfg;
}

XMOF2D::CellsMatData read_mat_data(const std::string& matdata_fname) {
  std::ifstream os(matdata_fname.c_str(), std::ifstream::binary);
  if (!os.good()) {
    std::ostringstream os;
    os << std::endl <<
      "Cannot open " << matdata_fname << " for binary input" << std::endl;
    throw XMOF2D::Exception(os.str());
  }
  
  XMOF2D::CellsMatData mat_data;
  
  int data_dim;
  os.read(reinterpret_cast<char *>(&data_dim), sizeof(int));
  if(data_dim != 2) {
    std::ostringstream os;
    os << std::endl << "Material data should be for a 2D mesh!" << std::endl;
    throw XMOF2D::Exception(os.str());
  }

  int ncells;
  os.read(reinterpret_cast<char *>(&ncells), sizeof(int));

  mat_data.cells_materials.resize(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    int nmats;
    os.read(reinterpret_cast<char *>(&nmats), sizeof(int));
    mat_data.cells_materials[icell].resize(nmats);
    for (int im = 0; im < nmats; im++)
      os.read(reinterpret_cast<char *>(&mat_data.cells_materials[icell][im]), sizeof(int));
  }
  mat_data.cells_vfracs.resize(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    int nmats = mat_data.cells_materials[icell].size();
    mat_data.cells_vfracs[icell].resize(nmats);
    if (nmats == 1) {
      mat_data.cells_vfracs[icell][0] = 1.0;
      continue;
    }
    for (int im = 0; im < nmats; im++)
      os.read(reinterpret_cast<char *>(&mat_data.cells_vfracs[icell][im]), sizeof(double));
  }
  mat_data.cells_centroids.resize(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    int nmats = mat_data.cells_materials[icell].size();
    mat_data.cells_centroids[icell].resize(nmats);
    if (nmats == 1) {
      mat_data.cells_centroids[icell][0] = XMOF2D::BAD_POINT;
      continue;
    }
    for (int im = 0; im < nmats; im++) {
      double cen_x, cen_y;
      os.read(reinterpret_cast<char *>(&cen_x), sizeof(double));
      os.read(reinterpret_cast<char *>(&cen_y), sizeof(double));
      mat_data.cells_centroids[icell][im] = XMOF2D::Point2D(cen_x, cen_y);
    }
  }

  os.close();
  
  return mat_data;
}

int PosWRTLine(const XMOF2D::Point2D& p, double line_slope, double line_shift, double eps) {
  std::vector<double> lort_vec = {line_slope, -1.0};
  std::vector<double> p2l_vec = {p.x, p.y - line_shift};
  double prj = XMOF2D::ddot(lort_vec, p2l_vec);
  if (std::fabs(prj) < eps)
    return 0;
  return std::signbit(prj) ? -1 : 1;
}
