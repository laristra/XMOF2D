/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#ifndef mesh_h
#define mesh_h

#include "geometry.h"
#include "exception.h"
#include "node.h"
#include "face.h"
#include "cell.h"

namespace XMOF2D {

struct MeshConfig {
  std::vector<Point2D>            nodes_coords;     //coordinates of base mesh nodes
  std::vector< std::vector<int> > ifaces_nodes;     //indices of their nodes for base mesh faces
  std::vector< std::vector<int> > icells_faces;     //indices of their faces for base mesh cells
  std::vector<int>                cells_material;   //indices of their materials for base mesh cells 
                                                    //(for MMCs it's -1)
};

struct IRTolerances {
  IRTolerances() : dist_eps(1.0e-15), div_eps(1.0e-08), ddot_eps(1.0e-14),
                   area_eps(1.0e-15), ang_eps(1.0e-14), mof_max_iter(10000) {}

  double dist_eps;      //distance tolerance: two points within this distance are considered coincidental
  double div_eps;       //the min value of denominator before switching to bisections when finding line-segment intersections
  double ddot_eps;      //dot product tolerance: for values below this, vectors are considered to have the same direction
  double area_eps;      //area tolerance for nested disections
  double ang_eps;       //angle tolerance for nested disections
  int    mof_max_iter;  //maximum number of iterations when finding the optimal distance or the optimal angle in nested disections
};

class MeshBC {
protected:
  std::vector<Node>  nodes;
  std::vector<Face>  faces;
  std::vector<Cell>  cells;
  std::vector<int>   nodes_global_ind;
  std::vector<int>   faces_global_ind;
  std::vector<int>   cells_global_ind;
  std::vector<int>   cells_material;
  
  double dist_eps_;
  double div_eps_;
  double ddot_eps_;
  double area_eps_;
  double ang_eps_;
  int    mof_max_iter_;

  void set_nodes_ifaces();
  void set_nodes_icells();
  
public:
  virtual ~MeshBC() {}
  
  int nnodes() const { return (int) nodes.size(); }
  int nfaces() const { return (int) faces.size(); }
  int ncells() const { return (int) cells.size(); }

  double dist_eps() const { return dist_eps_; }
  double div_eps() const { return div_eps_; }
  double ddot_eps() const { return ddot_eps_; }
  double area_eps() const { return area_eps_; }
  double ang_eps() const { return ang_eps_; }
  int    mof_max_iter() const { return mof_max_iter_; }
  
  int nmaterials () const;
  
  const Point2D& get_node_crd(int inode) const;
  const Node& get_node(int inode) const;
  const Face& get_face(int iface) const;
  const Cell& get_cell(int icell) const;
  int get_cells_material(int icell) const;
  
  void set_node_crd(int inode, const Point2D& crd);
  std::vector<int> get_bad_cells() const;
  
  bool is_ccw(const std::vector<int> inodes, double ddot_eps = 1.0e-14) const;
  bool on_same_line(const std::vector<int> inodes, double ddot_eps = 1.0e-14) const;
  std::pair< std::vector<int>, std::vector<bool> > get_ccw_faces_ordering(
    const std::vector<int>& ifaces);
  
  void dump_gmv(const std::string& filename) const;
  friend std::ostream& operator<<(std::ostream& os, const MeshBC& mesh);
  friend std::ofstream& operator<<(std::ofstream& os, const MeshBC& mesh);
};

class MiniMesh : public MeshBC {
private:
  const Cell& base_cell;
  Point2D nodes_shift;
  std::vector<int> orig_icells;
  
  void translate_2origin();
  void translate_2orig_pos();
  
  void Split(int cell_ind, const std::vector<double>& n, double d2orgn, int cell_mat);
  void ConstructMinimesh(const std::vector<double>& a2OX, const std::vector<double>& d2orgn, const std::vector<int>& cells_mat);  
public:
  MiniMesh(const Cell& base_cell_, double div_eps, double ddot_eps, double dist_eps);
  const Cell& get_base_cell() const { return base_cell; }
  const MeshBC& get_parent_mesh() const;

  void set_orig_icells(const std::vector<int>& opt_mat_order) { orig_icells = opt_mat_order; }
  const std::vector<int>& get_orig_icells() const { return orig_icells; }

  friend class Cell;
};

class BaseMesh : public MeshBC {
private:

public:
  BaseMesh(const MeshConfig& config, const IRTolerances& tolerances);
  void ConstructMinimesh(int icell, const std::vector<int>& imats, const std::vector<double>& vfracs, const std::vector<Point2D>& centroids);
};

class Factory {
private:
  const MeshBC& mesh;
  std::vector<Face>& faces;
  std::vector< std::vector< std::pair<std::vector<int>, int> > > sorted;
  
public:
  Factory(const MeshBC& _mesh, std::vector<Face>& _faces);
  
  int get_face(int cell_ind, const std::vector<int>& face);
};

}
#endif /* mesh_h */
