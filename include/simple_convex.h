/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#ifndef simple_convex_h
#define simple_convex_h

#include "gobject.h"

namespace XMOF2D {

class SimpleConvex : public GObject2D {
protected:
  std::vector<Point2D> v;
  Point2D              shift;
  
  bool vrts_are_ccw(double dist_eps = 1.0e-15, double ddot_eps = 1.0e-14);

public:
  SimpleConvex();
  SimpleConvex(const std::vector<Point2D>& v, double ddot_eps = 1.0e-14);
  
  Point2D                     	vertex(int i) const;
  const std::vector<Point2D>&	vertices() const;
  Point2D                     	center() const;
  std::unique_ptr<GObject1D>  	face(int i) const;
  int                       	nfaces() const;
  virtual double            	size() const;
  virtual bool              	contains(const Point2D& p, double eps = 1.0e-14) const;
  bool                      	contains(const SimpleConvex& sc, double eps = 1.0e-14) const;
  double                    	dist(const Point2D& p) const;
  
  void translate2origin();
  void translate2orig_pos();
  bool is_boundary(const Point2D& p, double eps = 1.0e-15) const;
  bool is_interior(const Point2D& p, double dist_eps = 1.0e-15, double ddot_eps = 1.0e-14) const;
  std::vector<SimpleConvex> SimpleConvexCutByLine(double a2OX, double d2orgn,
                                                  double ddot_eps, double dist_eps,
                                                  bool ipts_check_on);
  double compute_cutting_dist(double a2OX, double cut_area,
                              double area_eps, double ddot_eps, double dist_eps, int max_iter);
  double centroid_error(Point2D ref_centroid);
  double compute_optimal_angle(double cut_area, Point2D ref_centroid,
                               double ang_eps, double area_eps,
                               double ddot_eps, double dist_eps, int max_iter);
  void compute_optimal_cuts(const std::vector<double>& cut_vol_fracs,
                            const std::vector<Point2D>& ref_centroids,
                            double ang_eps, double area_eps,
                            double ddot_eps, double dist_eps, int max_iter,
                            std::vector<int>& opt_mat_order,
                            std::vector<double>& opt_a2OX, std::vector<double>& opt_d2orgn);
};

}
#endif /* simple_convex_h */
