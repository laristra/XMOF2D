/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "single_mmc_wrapper.h"

int main(int argc, char** argv) {

  std::string gnuplot_mmc_fname = "gnuplot_script_mmc.gnu";
  std::string gnuplot_minimesh_fname = "gnuplot_script_minimesh.gnu";
  
  std::vector<XMOF2D::Point2D> base_cell_vrts = {
    XMOF2D::Point2D(8.517743e-01, 3.781231e-01),
    XMOF2D::Point2D(8.506274e-01, 3.789720e-01),
    XMOF2D::Point2D(8.486559e-01, 3.770629e-01),
    XMOF2D::Point2D(8.490196e-01, 3.748427e-01),
    XMOF2D::Point2D(8.505602e-01, 3.748087e-01),
    XMOF2D::Point2D(8.519694e-01, 3.762168e-01)
  };
  
  std::vector<int> material_ids = {2, 3, 5, 7};
  std::vector<double> mat_volume_fractions = {
    1.862369e-01, 4.742934e-01, 8.186527e-02, 2.576044e-01};
  
  std::vector<XMOF2D::Point2D> material_centroids = {
    XMOF2D::Point2D(8.501062e-01, 3.779714e-01),
    XMOF2D::Point2D(8.501696e-01, 3.767483e-01),
    XMOF2D::Point2D(8.493403e-01, 3.752849e-01),
    XMOF2D::Point2D(8.511479e-01, 3.761662e-01)
  };
  
  int nbase_vrts = base_cell_vrts.size();
  int nmaterials = material_ids.size();
  
  double* vrt_crds = new double[2*nbase_vrts];
  for (int ivrt = 0; ivrt < nbase_vrts; ivrt++) {
    vrt_crds[2*ivrt] = base_cell_vrts[ivrt].x;
    vrt_crds[2*ivrt + 1] = base_cell_vrts[ivrt].y;
  }
  
  double* mat_centroids = new double[2*nmaterials];
  for (int imat = 0; imat < nmaterials; imat++) {
    mat_centroids[2*imat] = material_centroids[imat].x;
    mat_centroids[2*imat + 1] = material_centroids[imat].y;
  }

  std::cout << "****************************************" << std::endl;
  std::cout << "Input data:" << std::endl;
  std::cout << "  Vertices of the base cell" << std::endl;
  for (int ivrt = 0; ivrt < nbase_vrts; ivrt++) 
    std::cout << "    Vertex #" << ivrt << " at (" << vrt_crds[2*ivrt] 
              << ", " << vrt_crds[2*ivrt + 1] << ")" << std::endl;
  std::cout << std::endl;
  std::cout << "  IDs of materials inside the base cell: ";
  for (int imat = 0; imat < nmaterials; imat++)
    std::cout << material_ids[imat] << " ";
  std::cout << std::endl;
  std::cout << "  Volume fractions of materials inside the base cell" << std::endl;
  for (int imat = 0; imat < nmaterials; imat++)
    std::cout << "    Material #" << material_ids[imat] << ": " 
              << mat_volume_fractions[imat] << std::endl;
  std::cout << std::endl;
  std::cout << "  Reference centroids of materials inside the base cell" << std::endl;
  for (int imat = 0; imat < nmaterials; imat++)
    std::cout << "    Material #" << material_ids[imat] << ": ("
              << mat_centroids[2*imat] << ", " 
              <<  mat_centroids[2*imat + 1] << ")" << std::endl;
  std::cout << std::endl;
  std::cout << "****************************************" << std::endl;

  std::string mmc_nodes_fname = "mmc_node_crds.txt";
  std::ofstream ofs;
  ofs.open(mmc_nodes_fname.c_str(), std::ofstream::out);
  if (!ofs.good()) {
    std::ostringstream error_msg_stream;
    error_msg_stream << std::endl << "Could not open " << mmc_nodes_fname
                     << " for writing!" << std::endl;
    throw XMOF2D::Exception(error_msg_stream.str());
  }
  for (int ivrt = 0; ivrt <= nbase_vrts; ivrt++)
    ofs << std::setprecision(16) << base_cell_vrts[ivrt%nbase_vrts].x << " "
        << base_cell_vrts[ivrt%nbase_vrts].y << std::endl;
  ofs.close();
  
  std::ofstream script_ofs;
  script_ofs.open(gnuplot_mmc_fname.c_str(), std::ofstream::out);
  if (!script_ofs.good()) {
    std::ostringstream error_msg_stream;
    error_msg_stream << std::endl << "Could not open " << gnuplot_mmc_fname
                     << " for writing!" << std::endl;
    throw XMOF2D::Exception(error_msg_stream.str());
  }
  
  for (int ivrt = 0; ivrt < nbase_vrts; ivrt++) {
    script_ofs << "set label \"" << ivrt << "\" at " << base_cell_vrts[ivrt].x
               << "," <<  base_cell_vrts[ivrt].y << " tc lt 4" << std::endl;
    XMOF2D::Point2D fcen = 0.5*(base_cell_vrts[ivrt] +
                                base_cell_vrts[(ivrt + 1)%nbase_vrts]);
    script_ofs << "set label \"" << ivrt << "\" at " << fcen.x
               << "," <<  fcen.y << " right tc lt 3" << std::endl;
  }
  script_ofs << "set size square" << std::endl;
  script_ofs << "plot \'" << mmc_nodes_fname << "\' w lp notitle" << std::endl;
  script_ofs.close();
  
  std::cout << "Setting coordinates of multi-material cell's vertices..." << std::endl;
  xmof2d_set_mmc_vertices(nbase_vrts, vrt_crds);
  std::cout << "Setting IDs, volume fractions, and centroids for materials..." << std::endl;
  xmof2d_set_materials_data(nmaterials, &material_ids[0], &mat_volume_fractions[0],
                            mat_centroids);
  delete [] vrt_crds;
  delete [] mat_centroids;
  
  double dist_tol, ddot_tol, area_tol, ang_tol;
  int max_iter_num;
  xmof2d_get_distance_tolerance(&dist_tol);
  std::cout << "Default distance tolerance: " << dist_tol << std::endl;
  xmof2d_get_area_tolerance(&area_tol);
  std::cout << "Default area tolerance: " << area_tol << std::endl;
  xmof2d_get_angle_tolerance(&ang_tol);
  std::cout << "Default angle tolerance: " << ang_tol << std::endl;
  xmof2d_get_max_iter_num(&max_iter_num);
  std::cout << "Default max number of iterations: " << max_iter_num << std::endl;
  
  area_tol = 1.0e-14;
  std::cout << "Relaxing area tolerance to " << area_tol << std::endl;
  xmof2d_set_area_tolerance(area_tol);
  
  std::cout << "Creating a reconstructor instance..." << std::endl;
  xmof2d_initialize_reconstructor();
  std::cout << "Constructing a minimesh..." << std::endl;
  bool rec_success;
  xmof2d_perform_reconstruction(&rec_success);
  if (!rec_success) {
    throw XMOF2D::Exception("Reconstruction has failed, exiting...");
  }
  std::cout << "****************************************" << std::endl;
  
  int ncells, nfaces, nnodes;
  xmof2d_mesh_get_ncells(&ncells);
  xmof2d_mesh_get_nfaces(&nfaces);
  xmof2d_mesh_get_nnodes(&nnodes);
  std::cout << "Resulting minimesh has " << ncells << " cells, "
            << nfaces << " faces, and " << nnodes << " nodes." << std::endl;
  std::cout << "****************************************" << std::endl;
  
  script_ofs.open(gnuplot_minimesh_fname.c_str(), std::ofstream::out);
  if (!script_ofs.good()) {
    std::ostringstream error_msg_stream;
    error_msg_stream << std::endl << "Could not open " << gnuplot_minimesh_fname
                     << " for writing!" << std::endl;
    throw XMOF2D::Exception(error_msg_stream.str());
  }
  
  std::vector<std::string> cells_fnames;
  
  std::cout << "Cells data:" << std::endl;
  for (int icell = 0; icell < ncells; icell++) {
    std::cout << "  Cell #" << icell << std::endl;
    
    int cell_mat_id;
    xmof2d_cell_get_mat_id(icell, &cell_mat_id);
    std::cout << "    Material ID: " << cell_mat_id << std::endl;
    
    double cell_area;
    xmof2d_cell_get_size(icell, &cell_area);
    std::cout << "    Area: " << cell_area << std::endl;
    
    double cell_center[2];
    xmof2d_cell_get_center(icell, &cell_center[0]);
    std::cout << "    Actual Centroid: " << "(" << cell_center[0] << ", "
              << cell_center[1] << ")" << std::endl;
    int iref = (int) (std::find(material_ids.begin(), material_ids.end(), cell_mat_id) -
                      material_ids.begin());
    double cen_diff = sqrt(pow(cell_center[0] - material_centroids[iref].x, 2) +
                           pow(cell_center[1] - material_centroids[iref].y, 2));
    std::cout << "    Discrepancy in centroids: " << cen_diff << std::endl;
    
    int cell_nfaces;
    xmof2d_cell_get_nfaces(icell, &cell_nfaces);
    std::cout << "    Number of faces: " << cell_nfaces << std::endl;
    
    int* cell_ifaces = new int[cell_nfaces];
    xmof2d_cell_get_face_ids(icell, cell_ifaces);
    std::cout << "    Indices of faces:";
    for (int iside = 0; iside < cell_nfaces; iside++)
      std::cout << " " << cell_ifaces[iside];
    std::cout << std::endl;
    delete [] cell_ifaces;
      
    int* cell_inodes = new int[cell_nfaces];
    xmof2d_cell_get_node_ids(icell, cell_inodes);
    std::cout << "    Indices of nodes:";
    for (int inode = 0; inode < cell_nfaces; inode++)
      std::cout << " " << cell_inodes[inode];
    std::cout << std::endl;
    
    script_ofs << "set label \"#" << icell << ",MID " << cell_mat_id << "\" at "
               << cell_center[0] << "," <<  cell_center[1]
               << " center tc lt 1" << std::endl;
    
    std::ostringstream cell_nodes_fname_stream;
    cell_nodes_fname_stream << "cell_" << icell << "_node_crds.txt";
    std::string cell_nodes_fname = cell_nodes_fname_stream.str();
    std::ofstream ofs;
    ofs.open(cell_nodes_fname.c_str(), std::ofstream::out);
    if (!ofs.good()) {
      std::ostringstream error_msg_stream;
      error_msg_stream << std::endl << "Could not open " << cell_nodes_fname
                       << " for writing!" << std::endl;
      throw XMOF2D::Exception(error_msg_stream.str());
    }
    for (int icn = 0; icn <= cell_nfaces; icn++) {
      double node_crd[2];
      xmof2d_node_get_coords(cell_inodes[icn%cell_nfaces], &node_crd[0]);
      ofs << std::setprecision(16) << node_crd[0] << " " << node_crd[1] << std::endl;
    }
    ofs.close();
    delete [] cell_inodes;
    
    cells_fnames.push_back(cell_nodes_fname);
    std::cout << std::endl;
  }
  std::cout << "****************************************" << std::endl;
  std::cout << "Faces data:" << std::endl;
  for (int iface = 0; iface < nfaces; iface++) {
    std::cout << "  Face #" << iface << std::endl;
    
    int itail, ihead;
    xmof2d_face_get_node_ids(iface, &itail, &ihead);
    std::cout << "    Index of the tail node: " << itail << std::endl;
    std::cout << "    Index of the head node: " << ihead << std::endl;
    
    int ileft_cell, iright_cell;
    xmof2d_face_get_cell_ids(iface, &ileft_cell, &iright_cell);
    std::cout << "    Index of the cell to the left: " << ileft_cell << std::endl;
    std::cout << "    Index of the cell to the right: ";
    if (iright_cell == -1)
      std::cout << "NONE";
    else
      std::cout << iright_cell;
    std::cout << std::endl;
    
    double face_length;
    xmof2d_face_get_size(iface, &face_length);
    std::cout << "    Length: " << face_length << std::endl;
    
    double face_center[2];
    xmof2d_face_get_center(iface, &face_center[0]);
    std::cout << "    Center: " << "(" << face_center[0] << ", "
              << face_center[1] << ")" << std::endl;
    
    bool is_boundary;
    xmof2d_face_is_boundary(iface, &is_boundary);
    std::cout << "    This face is ";
    if (!is_boundary)
      std::cout << "NOT ";
    std::cout << "boundary" << std::endl;
    
    int iparent;
    xmof2d_face_get_parent_face_id(iface, &iparent);
    std::cout << "    Parent face ID: ";
    if (iparent == -1)
      std::cout << "NONE";
    else
      std::cout << iparent;
    std::cout << std::endl << std::endl;
    
    script_ofs << "set label \"" << iface << "\" at "
               << face_center[0] << "," <<  face_center[1]
               << " right tc lt 3" << std::endl;
  }
  
  std::cout << "****************************************" << std::endl;
  std::cout << "Nodes data:" << std::endl;
  for (int inode = 0; inode < nnodes; inode++) {
    std::cout << "  Node #" << inode << std::endl;
    
    double node_crd[2];
    xmof2d_node_get_coords(inode, &node_crd[0]);
    std::cout << "    Coordinates: (" << node_crd[0] << ", "
                                    << node_crd[1] << ")" << std::endl;
    
    int node_ncells;
    xmof2d_node_get_ncells(inode, &node_ncells);
    std::cout << "    This is a node of " << node_ncells << " cell(s)" << std::endl;
    
    int* node_icells = new int[node_ncells];
    xmof2d_node_get_cell_ids(inode, node_icells);
    std::cout << "    Indices of cells:";
    for (int icell = 0; icell < node_ncells; icell++)
      std::cout << " " << node_icells[icell];
    std::cout << std::endl;
    delete [] node_icells;
    
    int node_nfaces;
    xmof2d_node_get_nfaces(inode, &node_nfaces);
    std::cout << "    This is a node of " << node_nfaces << " face(s)" << std::endl;
    
    int* node_ifaces = new int[node_nfaces];
    xmof2d_node_get_face_ids(inode, node_ifaces);
    std::cout << "    Indices of faces:";
    for (int iface = 0; iface < node_nfaces; iface++)
      std::cout << " " << node_ifaces[iface];
    std::cout << std::endl;
    delete [] node_ifaces;
    
    bool is_boundary;
    xmof2d_node_is_boundary(inode, &is_boundary);
    std::cout << "    This node is ";
    if (!is_boundary)
      std::cout << "NOT ";
    std::cout << "boundary" << std::endl;
    
    int iparent_node;
    xmof2d_node_get_parent_node_id(inode, &iparent_node);
    std::cout << "    Parent node ID: ";
    if (iparent_node == -1)
      std::cout << "NONE";
    else
      std::cout << iparent_node;
    std::cout << std::endl;
    
    int node_nparent_faces;
    xmof2d_node_get_nparent_faces(inode, &node_nparent_faces);
    std::cout << "    This node belongs to ";
    if (node_nparent_faces == 0)
      std::cout << "NO";
    else
      std::cout << node_nparent_faces;
    std::cout << " base face(s)" << std::endl;
    
    if (node_nparent_faces) {
      int* node_iparent_faces = new int[node_nparent_faces];
      xmof2d_node_get_parent_face_ids(inode, node_iparent_faces);
      std::cout << "    Indices of those base faces:";
      for (int iface = 0; iface < node_nparent_faces; iface++)
        std::cout << " " << node_iparent_faces[iface];
      std::cout << std::endl << std::endl;
      delete [] node_iparent_faces;
    }
    
    script_ofs << "set label \"" << inode << "\" at "
               << node_crd[0] << "," <<  node_crd[1] << " tc lt 4" << std::endl;
  }
  
  script_ofs << "set size square" << std::endl;
  //std::string cell_fnames = cell_fnames_stream.str();
  //cell_fnames.resize(cell_fnames.size() - 1);
  //script_ofs << "filenames = \"" << cell_fnames << "\"" << std::endl;
  //script_ofs << "plot for [ifile=1:words(filenames)] word(filenames,ifile) w l notitle" << std::endl;
  script_ofs << "plot";
  for (int ifile = 0; ifile < cells_fnames.size() - 1; ifile++)
    script_ofs << "  \"" << cells_fnames[ifile] << "\" w l notitle, \\" << std::endl;
  script_ofs << "  \"" << cells_fnames[cells_fnames.size() - 1] << "\" w l notitle" << std::endl;
  script_ofs.close();
  
  std::cout << "****************************************" << std::endl;
  std::cout << "Removing the reconstructor instance..." << std::endl;
  xmof2d_free_reconstructor();
  std::cout << "All done!" << std::endl;
  
  return 0;
}

