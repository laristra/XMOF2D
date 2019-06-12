/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#include "mesh.h"
#include "tools.h"
#include "simple_vector.h"
#include <fstream>
#include <algorithm>
#include <functional>

namespace XMOF2D {

int MeshBC::nmaterials() const {
  return *std::max_element(cells_material.begin(), cells_material.end());
}

const Point2D& MeshBC::get_node_crd(int inode) const {
  XMOF2D_ASSERT_SIZE_LESS(inode, nodes.size());
  return nodes[inode].get_crd();
}
const Node& MeshBC::get_node(int inode) const {
  XMOF2D_ASSERT_SIZE_LESS(inode, nodes.size());
  return nodes[inode];
}
const Face& MeshBC::get_face(int iface) const {
  XMOF2D_ASSERT_SIZE_LESS(iface, faces.size());
  return faces[iface];
}
const Cell& MeshBC::get_cell(int icell) const {
  XMOF2D_ASSERT_SIZE_LESS(icell, cells.size());
  return cells[icell];
}
int MeshBC::get_cells_material(int icell) const {
  XMOF2D_ASSERT_SIZE_LESS(icell, ncells());
  return cells_material[icell];
}

void MeshBC::set_node_crd(int inode, const Point2D& crd) {
  XMOF2D_ASSERT_SIZE_LESS(inode, nodes.size());
  nodes[inode].set_coord(crd);
  
  const std::vector<int>& inode_cells = get_node(inode).get_cells();
  for (int inc = 0; inc < inode_cells.size(); inc++)
    cells[inode_cells[inc]].delete_minimesh();
}
  
void MeshBC::set_nodes_ifaces() {
  std::vector< std::vector<int> > nodes_ifaces(nnodes());
  for (int iface = 0; iface < nfaces(); iface++)
    for (int ifn = 0; ifn < 2; ifn++) {
      int inode = get_face(iface).get_node_index(ifn);
      nodes_ifaces[inode].push_back(iface);
    }
  
  for (int inode = 0; inode < nnodes(); inode++) {
    std::function<bool(int, int)> cmp_fun = [&](int iface1, int iface2) {
      return this->get_node(inode).faces_in_ccw_order(iface1, iface2);
    };
    std::sort(nodes_ifaces[inode].begin(), nodes_ifaces[inode].end(),
              cmp_fun);
    nodes[inode].set_faces(nodes_ifaces[inode]);
  }
}

void MeshBC::set_nodes_icells() {
  std::vector< std::vector<int> > nodes_icells(nnodes());
  for (int icell = 0; icell < ncells(); icell++)
    for (int icn = 0; icn < get_cell(icell).nfaces(); icn++) {
      int inode = get_cell(icell).get_node_index(icn);
      nodes_icells[inode].push_back(icell);
    }
  
  for (int inode = 0; inode < nnodes(); inode++) {
    std::function<bool(int, int)> cmp_fun = [&](int icell1, int icell2) {
      return this->get_node(inode).cells_in_ccw_order(icell1, icell2);
    };
    std::sort(nodes_icells[inode].begin(), nodes_icells[inode].end(),
              cmp_fun);
    nodes[inode].set_cells(nodes_icells[inode]);
  }
}

std::vector<int> MeshBC::get_bad_cells() const {
  std::vector<int> ibadcells;
  for (int icell = 0; icell < ncells(); icell++) {
    if (!is_ccw(get_cell(icell).get_nodes()))
      ibadcells.push_back(icell);
  }
  
  return ibadcells;
}

bool MeshBC::is_ccw(const std::vector<int> inodes) const {
  int np = (int) inodes.size();
  XMOF2D_ASSERT(np > 2, "Sequence should contain at least three nodes");
  
  for (int i = 0; i < np; i++)  {
    if (!XMOF2D::is_ccw(nodes[inodes[i]].get_crd(), nodes[inodes[(i + 1)%np]].get_crd(),
                        nodes[inodes[(i + 2)%np]].get_crd(), dist_eps(), ddot_eps()))
      return false;
  }
  
  return true;
}
  
bool MeshBC::on_same_line(const std::vector<int> inodes) const {
  int np = (int) inodes.size();
  XMOF2D_ASSERT(np > 2, "Sequence should contain at least three nodes");
  
  for (int i = 0; i < np - 2; i++)  {
    if (!XMOF2D::on_same_line(nodes[inodes[0]].get_crd(), nodes[inodes[i + 1]].get_crd(),
                              nodes[inodes[i + 2]].get_crd(), dist_eps(), ddot_eps()))
      return false;
  }
  
  return true;
}
  
std::pair< std::vector<int>, std::vector<bool> > MeshBC::get_ccw_faces_ordering(const std::vector<int>& ifaces) {
  int nsides = (int) ifaces.size();
  XMOF2D_ASSERT(nsides > 2, "Does not make sense for less than three faces!");
  
  std::vector< std::vector< std::pair<int, int> > > sconn(nsides,
                                                          {std::make_pair(-1, -1), std::make_pair(-1, -1)});
  for (int i = 0; i < nsides; i++) {
    const std::vector<int> &icur_side_nodes = get_face(ifaces[i]).get_nodes();
    int nconn = 0;
    for (int j = i + 1; j < nsides; j++) {
      const std::vector<int> &icmp_side_nodes = get_face(ifaces[j]).get_nodes();
      std::vector< std::pair<int, int> > isame_vrts = isamei(icur_side_nodes, icmp_side_nodes);
      XMOF2D_ASSERT_SIZE_LESS(isame_vrts.size(), 2);
      if (!isame_vrts.empty()) {
        int ivrt_curf = isame_vrts[0].first;
        int ivrt_cmpf = isame_vrts[0].second;
        sconn[i][ivrt_curf].first = j;
        sconn[i][ivrt_curf].second = ivrt_cmpf;
        sconn[j][ivrt_cmpf].first = i;
        sconn[j][ivrt_cmpf].second = ivrt_curf;
        nconn++;
      }
      if (nconn > 1)
        break;
    }
  }
  int istarts = 0;
  int istart_svrt = 0;
  int nhanging_nodes = 0;
  for (int iside = 0; iside < nsides; iside++)
    for (int ivrt = 0; ivrt < 2; ivrt++)
      if (sconn[iside][ivrt].first == -1) {
        istarts = iside;
        istart_svrt = ivrt;
        nhanging_nodes++;
      }
  XMOF2D_ASSERT(nhanging_nodes < 3, "Sequence of faces should be connected!");
  if (nhanging_nodes > 0) {
    XMOF2D_ASSERT(nhanging_nodes > 1, "Sequence of faces should be open on both ends!");
  }
  
  std::vector<int> inode_seq(nsides + 1);
  inode_seq[0] = get_face(ifaces[istarts]).get_nodes()[istart_svrt];
  std::vector<int> iface_seq(nsides);
  int inexts = istarts;
  int inext_svrt = (istart_svrt + 1)%2;
  std::vector<bool> side_is_ccw(nsides, true);
  for (int iside = 0; iside < nsides; iside++) {
    iface_seq[iside] = ifaces[inexts];
    inode_seq[iside + 1] = get_face(ifaces[inexts]).get_nodes()[inext_svrt];
    if (iside > 0) {
      int iprev_svrt = (inext_svrt + 1)%2;
      side_is_ccw[iside] = (iprev_svrt == sconn[inexts][iprev_svrt].second) ? !side_is_ccw[iside - 1] : side_is_ccw[iside - 1];
    }
    int iprevs = inexts;
    inexts = sconn[iprevs][inext_svrt].first;
    inext_svrt = (sconn[iprevs][inext_svrt].second + 1)%2;
    
  }
  if (nhanging_nodes == 0) {
    XMOF2D_ASSERT(inode_seq[0] == inode_seq[nsides], "Indices of end nodes for closed sequence do not match!")
  }
  
  std::vector<int> orig2red(nsides, 0);
  std::vector<int> inode_redseq = {inode_seq[0]};
  for (int iside = 0; iside < nsides - 1; iside++) {
    std::vector<int> inodes_2sides = {inode_seq[iside], inode_seq[iside + 1], inode_seq[iside + 2]};
    if (on_same_line(inodes_2sides))
      orig2red[iside + 1] = orig2red[iside];
    else {
      orig2red[iside + 1] = (int) inode_redseq.size();
      inode_redseq.push_back(inode_seq[iside + 1]);
    }
  }
  if (nhanging_nodes > 0)
    inode_redseq.push_back(inode_seq[nsides]);
  
  int nsides_red = (int) inode_redseq.size() - 1;
  if (nsides_red < nsides) {
    int ifirst_subside = 0;
    for (int iside = 0; iside < nsides; iside++) {
      if ((iside == nsides - 1) || (orig2red[iside + 1] != orig2red[iside])){
        if (ifirst_subside == iside)
          ifirst_subside++;
        else {
          const Point2D& first_node = get_node_crd(inode_seq[ifirst_subside]);
          std::vector<double> side_vec = get_node_crd(inode_seq[iside + 1]) - first_node;
          double side_vec_size = dnrm2(side_vec);
          for (int issnode = 0; issnode < iside - ifirst_subside; issnode++) {
            int icur_node = inode_seq[ifirst_subside + issnode + 1];
            double d2first = dnrm2(get_node_crd(icur_node) - first_node);
            Point2D corrected_crd(first_node.x + d2first*side_vec[0]/side_vec_size,
                                  first_node.y + d2first*side_vec[1]/side_vec_size);
            nodes[icur_node].set_coord(corrected_crd);
            ifirst_subside = iside + 1;
          }
        }
      }
    }
  }
  
  if (!is_ccw(inode_redseq)) {
    std::reverse(iface_seq.begin(), iface_seq.end());
    std::reverse(side_is_ccw.begin(), side_is_ccw.end());
    for (int iside = 0; iside < nsides; iside++)
      side_is_ccw[iside] = !side_is_ccw[iside];
  }
  
  return std::make_pair(iface_seq, side_is_ccw);
}
  
std::ostream& operator<<(std::ostream& os, const MeshBC& mesh) {
  os << "Number of nodes: " << mesh.nnodes() << std::endl;
  os << "Number of faces: " << mesh.nfaces() << std::endl;
  os << "Number of cells: " << mesh.ncells() << std::endl;

  for (int icell = 0; icell < mesh.cells.size(); icell++) {
    os << "Cell #" << icell << std::endl;
    os << mesh.cells[icell];
  }

  return os;
}

std::ofstream& operator<<(std::ofstream& ofs, const MeshBC& mesh) {
  ofs << std::scientific;
  ofs.precision(16);
  
  ofs << "# Number of nodes:" << std::endl << mesh.nnodes() << std::endl;
  ofs << "# Number of faces:" << std::endl << mesh.nfaces() << std::endl;
  ofs << "# Number of cells:" << std::endl << mesh.ncells() << std::endl;
  
  ofs << std::endl << "# Nodes:" << std::endl;
  for (int inode = 0; inode < mesh.nnodes(); inode++)
    ofs << inode << "\t" << mesh.get_node_crd(inode).x << "\t" << mesh.get_node_crd(inode).y << std::endl;
  
  ofs << std::endl << "# Faces:" << std::endl;
  for (int iface = 0; iface < mesh.nfaces(); iface++)
    ofs << mesh.get_face(iface) << std::endl;
  
  ofs << std::endl << "# Cells:" << std::endl;
  for (int icell = 0; icell < mesh.ncells(); icell++)
    ofs << mesh.get_cell(icell) << std::endl;
  
  return ofs;
}

void MeshBC::dump_gmv(const std::string& filename) const {
  std::ofstream os(filename.c_str());
  if (!os.good())
  THROW_EXCEPTION("Cannot open " << filename << " for writing");
  
  os << "gmvinput ascii" << std::endl;
  
  os << "nodev " << nnodes() << std::endl;
  for (int inode = 0; inode < nnodes(); inode++)
    os << nodes[inode].get_crd().x << " " << nodes[inode].get_crd().y << " " << 0 << std::endl;
  os << std::endl;
  
  os << "cells " << ncells() << std::endl;
  for (int icell = 0; icell < ncells(); icell++) {
    const int nfaces = get_cell(icell).nfaces();
    os << "general " << 1 << std::endl;
    const std::vector<int>& cfaces = get_cell(icell).get_faces();
    int nnodes_out = 0;
    for(int iface = 0; iface < nfaces; iface++)
      nnodes_out += get_face(cfaces[iface]).nnodes() - 1;
    os << nnodes_out << std::endl;
    for(int iface = 0; iface < nfaces; iface++) {
      int nfnodes = get_face(cfaces[iface]).nnodes();
      if (get_face(cfaces[iface]).normal_is_out(icell))
        for(int j = 0; j < nfnodes - 1; j++)
          os << get_face(cfaces[iface]).get_node_index(j) + 1 << " ";
      else
        for(int j = 0; j < nfnodes - 1; j++)
          os << get_face(cfaces[iface]).get_node_index(nfnodes - 1 - j) + 1 << " ";
    }
    os << std::endl;
  }
  
  int nmats = nmaterials();
  if (nmats != 1) {
    os << std::endl;
    os << "material " << nmats << " 0" << std::endl;
    for (int i = 0; i < nmats; i++)
      os << "mat" << i + 1 << " ";
    os << std::endl;
    for (int i = 0; i < cells_material.size(); i++)
      os << cells_material[i] + 1 << " ";
    os << std::endl;
  }
  os << std::endl;
  
  os << "tracers " << ncells() << std::endl;
  for (int icell = 0; icell < ncells(); icell++)
    os << cells[icell].center().x << " ";
  os << std::endl;
  for (int icell = 0; icell < ncells(); icell++)
    os << cells[icell].center().y << " ";
  os << std::endl;
  for (int icell = 0; icell < ncells(); icell++)
    os << 0 << " ";
  os << std::endl;
  os << "icells" << std::endl;
  for (int icell = 0; icell < ncells(); icell++)
    os << icell + 1 << " ";
  os << std::endl << "endtrace" << std::endl << std::endl;
  
  os << std::endl << "endgmv" << std::endl;
  
  os.close();
}

}
