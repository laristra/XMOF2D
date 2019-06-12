/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#include "segment.h"
#include "exception.h"
#include "simple_vector.h"
#include "tools.h"

namespace XMOF2D {

Segment::Segment(const Point2D& a, const Point2D& b) {
  v[0] = a; v[1] = b;
}

Segment::Segment(const std::vector<Point2D>& points) {
  XMOF2D_ASSERT_SIZE(points.size(), 2);
  v[0] = points[0]; v[1] = points[1];
}

const Segment& Segment::operator=(const Segment& seg) {
  if (this == &seg) return *this;
  v[0] = seg[0]; v[1] = seg[1];
  
  return *this;
}

Point2D Segment::operator[](int i) const {
  XMOF2D_ASSERT_SIZE_LESS(i, 2);
  return v[i];
}

double Segment::size() const {
  return distance(v[0], v[1]);
}

Point2D Segment::middle() const {
  return 0.5*(v[0] + v[1]);
}

std::vector<double> Segment::direction() const {
  std::vector<double> direction = v[1] - v[0];
  dscal(1.0/dnrm2(direction), direction);
  return direction;
}

std::vector<double> Segment::normal() const {
  double len = distance(v[0], v[1]);
  std::vector<double> normal = {
    (v[1].y - v[0].y)/len, -(v[1].x - v[0].x)/len };
  return normal;
}

double Segment::dist(const Point2D& p) const {
  std::vector<double> dirv = direction();
  double prj_shift = ddot(dirv, p - v[0]);
  if ((prj_shift < 0.0) || is_equal(prj_shift, 0.0))
    return distance(p, v[0]);
  if ((prj_shift > size()) || is_equal(prj_shift, size()))
    return distance(p, v[1]);
  dscal(prj_shift, dirv);
  return distance(p, v[0] + dirv);
}

bool Segment::contains(const Point2D& p, double eps) const {
  return (std::fabs(distance(v[0], v[1]) - (distance(v[0], p) + distance(p, v[1]))) < eps);
}

SegLine::Position Segment::PosWRT2Line(std::vector<double> n, double d2orgn, double eps) {
  Point2DLine::Position v_pos[2];
  for (int inode = 0; inode < 2; inode++)
    v_pos[inode] = v[inode].PosWRT2Line(n, d2orgn, eps);
  
  SegLine::Position pos;
  if (v_pos[0] == Point2DLine::ON) {
    if (v_pos[1] == Point2DLine::ON)
      pos = SegLine::Position::ON;
    else
      pos = SegLine::Position::TOUCHES;
  }
  else {
    if (v_pos[1] == Point2DLine::Position::ON)
      pos = SegLine::Position::TOUCHES;
    else {
      if (v_pos[0] != v_pos[1])
        pos = SegLine::Position::INTERSECTS;
      else {
        if (v_pos[0] == Point2DLine::Position::ABOVE)
          pos = SegLine::Position::ABOVE;
        else
          pos = SegLine::Position::BELOW;
      }
    }
  }
  
  return pos;
}

Point2D Segment::LineIntersect(std::vector<double> n, double d2orgn, double denom_eps, double eps) {
  SegLine::Position pos = PosWRT2Line(n, d2orgn, eps);
  
  XMOF2D_ASSERT((pos == SegLine::Position::TOUCHES) || (pos == SegLine::Position::INTERSECTS),
         "Segment is above, below, or on the line!");
  
  Point2D pint = BAD_POINT;
  if (pos == SegLine::Position::TOUCHES) {
    double lv[2];
    for (int inode = 0; inode < 2; inode++) {
      std::vector<double> v2l = v[inode].vec();
      daxpy(-d2orgn, n, v2l);
      lv[inode] = ddot(v2l, n);
    }
    pint = (std::fabs(lv[0]) < eps) ? v[0] : v[1];
  }
  else {
    std::vector<double> sdv = v[1] - v[0];
    double cpz = std::fabs(ddot(sdv, n));
    
    if(std::fabs(cpz) > denom_eps) {
      double t = std::fabs(-d2orgn + ddot(v[0].vec(), n)) / cpz;
      XMOF2D_ASSERT((t > 0.0) && (t < 1.0), "Intersection point should belong to the interior!");
      dscal(t, sdv);
      pint = v[0] + sdv;
    }
    else {
      Segment bisected_seg = *this;
      SegLine::Position bs_pos = SegLine::Position::INTERSECTS;
      do {
        Segment test_seg = Segment(bisected_seg[0], bisected_seg.middle());
        SegLine::Position cur_pos = test_seg.PosWRT2Line(n, d2orgn, eps);
        if ((cur_pos == SegLine::Position::BELOW) || (cur_pos == SegLine::Position::ABOVE))
          bisected_seg = Segment(bisected_seg.middle(), bisected_seg[1]);
        else
          bisected_seg = test_seg;
        bs_pos = bisected_seg.PosWRT2Line(n, d2orgn, eps);
      } while (bs_pos == SegLine::Position::INTERSECTS);
      pint = bisected_seg.LineIntersect(n, d2orgn, denom_eps, eps);
    }
  }
  
  return pint;
}

std::ostream& operator<<(std::ostream& os, const Segment& seg) {
  os << "[" << seg.v[0] << ", " << seg.v[1] << "]";
  return os;
}

std::vector<int> isame_vrt(const Segment& a, const Segment& b) {
  std::vector<int> isamev(2, -1);
  for (int iva = 0; iva < 2; iva++)
    for (int ivb = 0; ivb < 2; ivb++)
      if (a[iva] == b[ivb]) {
        isamev[0] = iva;
        isamev[1] = ivb;
        return isamev;
      }
  return isamev;
}

bool is_ccw(const Segment& a, const Segment& b, double dist_eps, double ddot_eps) {
  std::vector<int> isamev = isame_vrt(a, b);
  if (isamev[0] == -1)
    return false;
  std::vector<Point2D> p = { a[(isamev[0] + 1)%2], a[isamev[0]], b[(isamev[1] + 1)%2] };
  return is_ccw(p[0], p[1], p[2], dist_eps, ddot_eps);
}

bool is_ccw(const std::vector<Segment>& segs, double dist_eps, double ddot_eps) {
  int nseg = (int) segs.size();
  XMOF2D_ASSERT(nseg > 1, "Sequence should contain at least two segments");
  
  for (int i = 0; i < nseg - 1; i++)  {
    if (!is_ccw(segs[i], segs[i + 1], dist_eps, ddot_eps))
      return false;
  }
  
  return true;
}

}
