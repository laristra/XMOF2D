/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#include <numeric>
#include <algorithm>
#include <iostream>
#include "simple_convex.h"
#include "exception.h"
#include "tools.h"
#include "logger.h"
#include "simple_vector.h"
#include "segment.h"
#include <iomanip>
namespace XMOF2D {

static inline double tri_size(const Point2D& a, const Point2D& b, const Point2D& c) {
  double v1[2] = { b.x - a.x, b.y - a.y };
  double v2[2] = { c.x - a.x, c.y - a.y };
  
  return 0.5*std::fabs(v1[0]*v2[1] - v1[1]*v2[0]);
}

SimpleConvex::SimpleConvex() {
  shift = BAD_POINT;
}

SimpleConvex::SimpleConvex(const std::vector<Point2D>& p, double ddot_eps) {
  v = p;
  shift = Point2D(0.0, 0.0);
  XMOF2D_ASSERT(vrts_are_ccw(ddot_eps), "Tried to construct a non-convex polygon!");
}

bool SimpleConvex::vrts_are_ccw(double ddot_eps) {
  return is_ccw(v, ddot_eps);
}

Point2D SimpleConvex::vertex(int i) const {
  XMOF2D_ASSERT_SIZE_LESS(i, nfaces());
  return v[i];
}

const std::vector<Point2D>& SimpleConvex::vertices() const {
  return v;
}

int SimpleConvex::nfaces() const {
  return (int) v.size();
}

std::unique_ptr<GObject1D> SimpleConvex::face(int i) const {
  const int n = nfaces();
  XMOF2D_ASSERT_SIZE_LESS(i, n);
  return std::unique_ptr<GObject1D>(new Segment(v[i], v[(i+1)%n]));
}

double SimpleConvex::size() const {
  const int n = nfaces();
  double size = 0;
  for (int i = 0; i < n; i++) {
    int ifv = i, isv = (i + 1) % n;
    size += v[ifv].x * v[isv].y - v[isv].x * v[ifv].y;
  }
  size *= 0.5;
  
  if (is_equal(size, 0.0))
  return 0.0;
  
  XMOF2D_ASSERT(size > 0.0, "Computed size is negative!");
  
  return size;
}

Point2D SimpleConvex::center() const {
  const int n = nfaces();
  double size = 0.0;
  Point2D centroid(0.0, 0.0);
  for (int i = 0; i < n; i++) {
    int ifv = i, isv = (i + 1)%n;
    double cur_size_term = v[ifv].x * v[isv].y - v[isv].x * v[ifv].y;
    size += cur_size_term;
    centroid.x += cur_size_term * (v[ifv].x + v[isv].x);
    centroid.y += cur_size_term * (v[ifv].y + v[isv].y);
  }
  size *= 0.5;
  
  if (is_equal(size, 0.0))
  return BAD_POINT;
  
  XMOF2D_ASSERT(size > 0.0, "Computed size is negative!");
  
  return centroid/(6*size);
}

bool SimpleConvex::contains(const Point2D& p, double eps) const {
  const int n = nfaces();
  double new_size = 0;
  for (int i = 0; i < n; i++)
    new_size += tri_size(p, v[i], v[(i + 1)%n]);

  return std::fabs(size() - new_size) < eps;
}

bool SimpleConvex::contains(const SimpleConvex& sc, double eps) const {
  const int m = sc.nfaces();
  for (int i = 0; i < m; i++)
    if (contains(sc.v[i], eps) == false)
      return false;
  return true;
}

double SimpleConvex::dist(const Point2D& p) const {
  if (contains(p))
    return 0.0;
  
  const int n = nfaces();
  std::vector<double> d(n);
  for (int i = 0; i < n; i++)
    d[i] = Segment(v[i], v[(i + 1)%n]).dist(p);
  
  return *std::min_element(d.begin(), d.end());
}

void SimpleConvex::translate2origin() {
  XMOF2D_ASSERT(shift == Point2D(0.0, 0.0), "This convex polygon has been already translated to the origin!");
  
  shift = BAD_POINT;
  for (int inode = 0; inode < nfaces(); inode++) {
    if (vertex(inode).x < shift.x)
      shift.x = vertex(inode).x;
    if (vertex(inode).y < shift.y)
      shift.y = vertex(inode).y;
  }
  shift = -shift;
  for (int inode = 0; inode < nfaces(); inode++)
    v[inode] += shift;
}

void SimpleConvex::translate2orig_pos() {
  shift = -shift;
  for (int inode = 0; inode < nfaces(); inode++)
    v[inode] += shift;
  
  shift = Point2D(0.0, 0.0);
}

bool SimpleConvex::is_boundary(const Point2D& p, double eps) const {
  int nnodes = nfaces();
  for (int inode = 0; inode < nnodes; inode++) {
    Segment cur_side(v[inode], v[(inode + 1)%nnodes]);
    if (cur_side.contains(p, eps))
      return true;
  }
  return false;
}

bool SimpleConvex::is_interior(const Point2D& p, double ddot_eps) const {
  int wn = 0;
  int nnodes = nfaces();
  for (int inode = 0; inode < nnodes; inode++) {
    int inextnode = (inode + 1)%nnodes;
    if (v[inode].y <= p.y) {
      if (v[inextnode].y  > p.y)
        if (is_ccw(v[inode], v[inextnode], p, ddot_eps))
          ++wn;
    }
    else {
      if (v[inextnode].y <= p.y)
        if (is_cw(v[inode], v[inextnode], p, ddot_eps))
          --wn;
    }
  }
  return wn;
}

std::vector<SimpleConvex> SimpleConvex::SimpleConvexCutByLine(double a2OX, const double d2orgn,
                                                              const double ddot_eps,
                                                              const double dist_eps,
                                                              bool area_check_on) {
  std::vector<double> n(2);
  n[0] = cos(a2OX);
  n[1] = sin(a2OX);
  double n_scal = 1.0/dnrm2(n);
  dscal(n_scal, n);
  
  int nsides = nfaces();
  std::vector<Point2DLine::Position> vrts_pos(nsides);
  for (int ivrt = 0; ivrt < nsides; ivrt++)
    vrts_pos[ivrt] = vertex(ivrt).PosWRT2Line(n, d2orgn, dist_eps);

  std::vector<Point2D> v_below;
  std::vector<Point2D> v_above;
  int int_point_ind[2] = {-1, -1};
  int nearest_lhp_ivrt[2] = {-1, -1};

  for (int iside = 0; iside < nsides; iside++) {
    int iv0 = iside, iv1 = (iside + 1)%nsides;
    switch(vrts_pos[iv0]) {
      case Point2DLine::Position::BELOW: 
        v_below.push_back(vertex(iv0)); break;
      case Point2DLine::Position::ABOVE: 
        v_above.push_back(vertex(iv0)); break;
      case Point2DLine::Position::ON:
        if (vrts_pos[iv0] == vrts_pos[iv1]) {
          THROW_EXCEPTION("Side lies on the cutting line!");
        }
        else {
          //Line intersects the side at the first node
          XMOF2D_ASSERT(int_point_ind[1] == -1, "Extra intersection!");
          int iip = (int_point_ind[0] == -1) ? 0 : 1;
          int_point_ind[iip] = (int) v_below.size();
          v_below.push_back(vertex(iv0));
          v_above.push_back(vertex(iv0));
          continue;
        }
      default: THROW_EXCEPTION("Unknown point with respect to line position!");
    }

    //Check if line intersects the side at the second node
    if (vrts_pos[iv1] == Point2DLine::Position::ON) continue;
    
    //Check if line intersects the interior of the side
    if (vrts_pos[iv0] != vrts_pos[iv1]) {
      XMOF2D_ASSERT(int_point_ind[1] == -1, "Extra intersection!");

      Segment cur_side = Segment(vertex(iv0), vertex(iv1));
      Point2D int_point = cur_side.LineIntersect(n, d2orgn);
      int iip = (int_point_ind[0] == -1) ? 0 : 1;

      if (area_check_on) {          
        nearest_lhp_ivrt[iip] = 
          (vrts_pos[iv0] == Point2DLine::Position::BELOW) ? iv0 : iv1;

        if (nearest_lhp_ivrt[0] == nearest_lhp_ivrt[1]) {
          XMOF2D_ASSERT(std::fabs(vpz(vertex(nearest_lhp_ivrt[1]), int_point, 
                                      v_below[int_point_ind[0]])) >= 2*epsilon,
                        "Area of the cutoff is below the machine epsilon!");
        }
      }
      
      int_point_ind[iip] = (int) v_below.size();
      v_below.push_back(int_point);
      v_above.push_back(int_point);
    }
  }
  
  XMOF2D_ASSERT(int_point_ind[1] != -1, 
                "Cutting line does not pass through the interior!");
  
  std::vector<SimpleConvex> SC_cuts(2);
  SC_cuts[0] = SimpleConvex(v_below, ddot_eps);
  SC_cuts[1] = SimpleConvex(v_above, ddot_eps);
  
  return SC_cuts;
}

double SimpleConvex::compute_cutting_dist(double a2OX, double cut_area,
                                          double area_eps, double ddot_eps, double dist_eps,
                                          int max_iter) {
  std::vector<double> n(2);
  n[0] = cos(a2OX);
  n[1] = sin(a2OX);
  double n_scal = 1.0/dnrm2(n);
  dscal(n_scal, n);
  
  double full_size = size();

  std::vector<double> dbnd = { DBL_MAX, -DBL_MAX };
  for (int i = 0; i < nfaces(); i++) {
    double cur_d = ddot(n, v[i].vec());
    if (cur_d < dbnd[0]) dbnd[0] = cur_d;
    if (cur_d > dbnd[1]) dbnd[1] = cur_d;
  }
 
  if (cut_area < area_eps) return dbnd[0];
  if (cut_area > full_size - area_eps) return dbnd[1];

  double d2orgn = 0.5*(dbnd[0] + dbnd[1]);
  std::vector<double> area_bnd = { 0.0, full_size };
  double area_err = full_size, darea = full_size;
  std::vector<double> cur_dbnd = dbnd;

  bool use_secant = true;
  int iter_count = 0;
  
  while (area_err > area_eps) {
    double sec_coef;
    if (use_secant) {
      // Both the change in area on the previous step and the max possible change
      // in area on this step should be above volume tolerance to keep using 
      // the secant method
      use_secant = (std::fabs(darea) > area_eps) && 
                   (std::fabs(area_bnd[1] - area_bnd[0]) > area_eps);
      if (use_secant) {
        sec_coef = (cur_dbnd[1] - cur_dbnd[0])/(area_bnd[1] - area_bnd[0]);
        d2orgn = cur_dbnd[1] + sec_coef*(cut_area - area_bnd[1]);
        use_secant = ((d2orgn > dbnd[0]) == (d2orgn < dbnd[1]));          
      }

      if(!use_secant) {
        // Secant method failed: fall back to the bisection algorithm
        bool valid_lower_bnd = (area_bnd[0] < cut_area);
        if ( valid_lower_bnd != (area_bnd[1] > cut_area) ) {
          // Target volume fraction is no longer within the current bounds: 
          // reset current bounds
          cur_dbnd = dbnd;
        }
        else if (!valid_lower_bnd) {
          // frac_bnd[0] is the upper bound: swap the current bounds
          std::reverse(cur_dbnd.begin(), cur_dbnd.end());
        }

        d2orgn = 0.5*(cur_dbnd[0] + cur_dbnd[1]);
      }
    }
    else
      d2orgn = 0.5*(cur_dbnd[0] + cur_dbnd[1]);
    
    std::vector<SimpleConvex> SC_cuts;
    double cur_area;
    try {
      SC_cuts = SimpleConvexCutByLine(a2OX, d2orgn, ddot_eps, dist_eps, true);
    }
    catch (Exception e) {
      if (e.is("Side lies on the cutting line!") ||
          e.is("Area of the cutoff is below the machine epsilon!") ||
          e.is("Cutting line does not pass through the interior!")) {
        //Distance check with the original bounds: 
        //we either sliced off nothing or everything
        cur_area = (std::fabs(d2orgn - dbnd[0]) < std::fabs(d2orgn - dbnd[1])) ? 
          0.0 : full_size;
        darea = cur_area - area_bnd[1];

        cur_dbnd[0] = cur_dbnd[1];
        cur_dbnd[1] = d2orgn;
        area_bnd[0] = area_bnd[1];
        area_bnd[1] = cur_area;

        continue;                    
      } 
      else THROW_EXCEPTION(e.what());
    }

    cur_area = SC_cuts[0].size();
    area_err = std::fabs(cur_area - cut_area);     
    
    iter_count++;
    if (iter_count > max_iter) {
      std::clog << "Max iterations count exceeded when finding the cutting distance, " <<
        "achieved error in area is " << area_err << std::endl;
      break;
    }

    if (!use_secant) {
      if (cur_area > cut_area)
        cur_dbnd[1] = d2orgn;
      else
        cur_dbnd[0] = d2orgn;
    }
    else {
      darea = cur_area - area_bnd[1];

      cur_dbnd[0] = cur_dbnd[1];
      cur_dbnd[1] = d2orgn;
      area_bnd[0] = area_bnd[1];
      area_bnd[1] = cur_area;
    }
  }
  
  return d2orgn;
}

double SimpleConvex::centroid_error(Point2D ref_centroid) {
  return pow(distance(center(), ref_centroid), 2)/size();
}

double SimpleConvex::compute_optimal_angle(double cut_area, Point2D ref_centroid, double ang_eps, double area_eps, double ddot_eps, double dist_eps, int max_iter) {
  Point2D cur_centroid = center();
  std::vector<double> dcen = cur_centroid - ref_centroid;
  double a2OX = atan2(dcen[1], dcen[0]);
  double dangle = PI / 360.0;
  std::vector<double> a_guess = { a2OX - dangle, a2OX, a2OX + dangle };
  std::vector<double> guess_err(3);

  for (int ia = 0; ia < 3; ia++) {
    double d2orgn = compute_cutting_dist(a_guess[ia], cut_area, area_eps, ddot_eps, dist_eps, max_iter);
    std::vector<SimpleConvex> SC_cuts = SimpleConvexCutByLine(a_guess[ia], d2orgn, ddot_eps, dist_eps, true);
    guess_err[ia] = SC_cuts[0].centroid_error(ref_centroid);
    if (guess_err[ia] < dist_eps) return a_guess[ia];
  }
  
  int iter_count = 0;
  bool inside_bracket = false;
  while (!inside_bracket) {
    XMOF2D_ASSERT(iter_count < max_iter, "Max iterations count exceeded while bracketing the optimal angle!");
    
    if ( (std::fabs(guess_err[1] - guess_err[0]) <= dist_eps) ||
         (std::fabs(guess_err[1] - guess_err[2]) <= dist_eps) ) {
      a_guess[0] -= a_guess[1] - a_guess[0];
      a_guess[2] += a_guess[2] - a_guess[1];

      double d2orgn = compute_cutting_dist(a_guess[0], cut_area, area_eps, ddot_eps, dist_eps, max_iter);
      std::vector<SimpleConvex> SC_cuts = SimpleConvexCutByLine(a_guess[0], d2orgn, ddot_eps, dist_eps, true);
      guess_err[0] = SC_cuts[0].centroid_error(ref_centroid);
      if (guess_err[0] < dist_eps) return a_guess[0];

      d2orgn = compute_cutting_dist(a_guess[2], cut_area, area_eps, ddot_eps, dist_eps, max_iter);
      SC_cuts = SimpleConvexCutByLine(a_guess[2], d2orgn, ddot_eps, dist_eps, true);
      guess_err[2] = SC_cuts[0].centroid_error(ref_centroid);
      if (guess_err[2] < dist_eps) return a_guess[2];
    }
    else if ( (guess_err[0] < guess_err[1]) ||
              (guess_err[2] < guess_err[1]) ) {
      if (guess_err[0] < guess_err[2]) {
        std::rotate(a_guess.begin(), a_guess.end() - 1, a_guess.end());
        a_guess[0] = a_guess[1] - 2*(a_guess[2] - a_guess[1]);
        std::rotate(guess_err.begin(), guess_err.end() - 1, guess_err.end());
        
        double d2orgn = compute_cutting_dist(a_guess[0], cut_area, area_eps, ddot_eps, dist_eps, max_iter);
        std::vector<SimpleConvex> SC_cuts = SimpleConvexCutByLine(a_guess[0], d2orgn, ddot_eps, dist_eps, true);
        guess_err[0] = SC_cuts[0].centroid_error(ref_centroid);
        if (guess_err[0] < dist_eps) return a_guess[0];
      }
      else {
        std::rotate(a_guess.begin(), a_guess.begin() + 1, a_guess.end());
        a_guess[2] = a_guess[1] + 2*(a_guess[1] - a_guess[0]);
        std::rotate(guess_err.begin(), guess_err.begin() + 1, guess_err.end());
        
        double d2orgn = compute_cutting_dist(a_guess[2], cut_area, area_eps, ddot_eps, dist_eps, max_iter);
        std::vector<SimpleConvex> SC_cuts = SimpleConvexCutByLine(a_guess[2], d2orgn, ddot_eps, dist_eps, true);
        guess_err[2] = SC_cuts[0].centroid_error(ref_centroid);
        if (guess_err[2] < dist_eps) return a_guess[2];        
      }
    }
    else inside_bracket = true;

    if (a_guess[2] - a_guess[0] > 2*PI) inside_bracket = true;

    iter_count++;
  }

  double gratio = 0.5*(sqrt(5) - 1.0);
  double daguess = a_guess[2] - a_guess[0];
  std::vector<double> gr_guess = { a_guess[0] + (1.0 - gratio)*daguess, a_guess[0] + gratio*daguess };
  std::vector<double> gr_error(2);
  for (int ia = 0; ia < 2; ia++) {
    double d2orgn = compute_cutting_dist(gr_guess[ia], cut_area, area_eps, ddot_eps, dist_eps, max_iter);
    std::vector<SimpleConvex> SC_cuts = SimpleConvexCutByLine(gr_guess[ia], d2orgn, ddot_eps, dist_eps, true);
    gr_error[ia] = SC_cuts[0].centroid_error(ref_centroid);
  }
  for (iter_count = 0; iter_count < max_iter; iter_count++) {
    int imin_err = (gr_error[0] < gr_error[1]) ? 0 : 1;
    if ((gr_error[imin_err] < dist_eps) || (daguess < ang_eps))
      return gr_guess[imin_err];
    
    if (!imin_err) {
      a_guess[2] = gr_guess[1];
      daguess = a_guess[2] - a_guess[0];
      gr_guess[1] = gr_guess[0];
      gr_error[1] = gr_error[0];
      gr_guess[0] = a_guess[0] + (1.0 - gratio)*daguess;
      
      double d2orgn = compute_cutting_dist(gr_guess[0], cut_area, area_eps, ddot_eps, dist_eps, max_iter);
      std::vector<SimpleConvex> SC_cuts = SimpleConvexCutByLine(gr_guess[0], d2orgn, ddot_eps, dist_eps, true);
      gr_error[0] = SC_cuts[0].centroid_error(ref_centroid);
    }
    else {
      a_guess[0] = gr_guess[0];
      daguess = a_guess[2] - a_guess[0];
      gr_guess[0] = gr_guess[1];
      gr_error[0] = gr_error[1];
      gr_guess[1] = a_guess[0] + gratio*daguess;
      
      double d2orgn = compute_cutting_dist(gr_guess[1], cut_area, area_eps, ddot_eps, dist_eps, max_iter);
      std::vector<SimpleConvex> SC_cuts = SimpleConvexCutByLine(gr_guess[1], d2orgn, ddot_eps, dist_eps, true);
      gr_error[1] = SC_cuts[0].centroid_error(ref_centroid);
    }
  }
  return 0.5*(a_guess[0] + a_guess[2]);
}

void SimpleConvex::compute_optimal_cuts(const std::vector<double>& cut_vol_fracs, const std::vector<Point2D>& ref_centroids, double ang_eps, double area_eps, double ddot_eps, double dist_eps, int max_iter, std::vector<int>& opt_mat_order, std::vector<double>& opt_a2OX, std::vector<double>& opt_d2orgn) {
  translate2origin();
  
  int nmat = (int) cut_vol_fracs.size();
  std::vector<int> cur_mat_order(nmat);
  std::iota(cur_mat_order.begin(), cur_mat_order.end(), 0);
  
  double min_err = DBL_MAX;
  std::vector<double> cur_a2OX(nmat - 1);
  std::vector<double> cur_d2orgn(nmat - 1);
  do {
    double aggr_err = 0.0;
    SimpleConvex cur_SC = SimpleConvex(vertices());
    for(int imat = 0; imat < nmat - 1; imat++) {
      int cur_mat = cur_mat_order[imat];
      double cut_area = cut_vol_fracs[cur_mat] * size();
      
      cur_a2OX[imat] = cur_SC.compute_optimal_angle(cut_area, ref_centroids[cur_mat] + shift, ang_eps, area_eps, ddot_eps, dist_eps, max_iter);
      cur_d2orgn[imat] = cur_SC.compute_cutting_dist(cur_a2OX[imat], cut_area, area_eps, ddot_eps, dist_eps, max_iter);
      
      std::vector<SimpleConvex> SC_cuts = cur_SC.SimpleConvexCutByLine(cur_a2OX[imat], cur_d2orgn[imat], ddot_eps, dist_eps, true);
      
      aggr_err += SC_cuts[0].centroid_error(ref_centroids[cur_mat] + shift);
      cur_SC = SC_cuts[1];
    }
    aggr_err += cur_SC.centroid_error(ref_centroids[cur_mat_order[nmat - 1]] + shift);
    
    if (aggr_err < min_err) {
      min_err = aggr_err;
      opt_mat_order = cur_mat_order;
      opt_a2OX = cur_a2OX;
      opt_d2orgn = cur_d2orgn;
    }
  } while (std::next_permutation(cur_mat_order.begin(), cur_mat_order.end()));
  
  translate2orig_pos();
}

}
