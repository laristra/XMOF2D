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

const Point2D& Segment::operator[](int i) const {
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
  double seg_size = size();
  std::vector<double> normal = {
    (v[1].y - v[0].y)/seg_size, -(v[1].x - v[0].x)/seg_size };
  return normal;
}

double Segment::dist_to_seg(const Point2D& p) const {
  std::vector<double> dirv = direction();
  double prj_shift = ddot(dirv, p - v[0]);
  if ((prj_shift < 0.0) || is_equal(prj_shift, 0.0))
    return distance(p, v[0]);
  double seg_size = size();
  if ( (prj_shift > seg_size) || is_equal(prj_shift, seg_size) )
    return distance(p, v[1]);

  return std::fabs(vpz(v[0], p, v[1]))/seg_size;
}

double Segment::dist_to_line(const Point2D& p) const {
  return std::fabs(vpz(v[0], p, v[1]))/size();
}

bool Segment::contains(const Point2D& p, double dist_eps) const {
  double seg_size = size();

  for (int iv = 0; iv < 2; iv++) {
    double dpv = distance(v[iv], p);

    if (dpv < dist_eps) return true;
    if (dpv > seg_size) return false;
  }

  return (std::fabs(vpz(v[0], p, v[1]))/seg_size < dist_eps);
}

bool Segment::is_interior(const Point2D& p, double dist_eps) const {
  double seg_size = size();

  for (int iv = 0; iv < 2; iv++) {
    double dpv = distance(v[iv], p);

    if ( (dpv < dist_eps) || (dpv >= sqrt(pow2(seg_size) + pow2(dist_eps))) )
      return false;
  }

  return (std::fabs(vpz(v[0], p, v[1]))/seg_size < dist_eps);
}

SegLine::Position Segment::PosWRT2Line(std::vector<double> n, double d2orgn, 
                                       double dist_eps) const{
  Point2DLine::Position v_pos[2];
  for (int inode = 0; inode < 2; inode++)
    v_pos[inode] = v[inode].PosWRT2Line(n, d2orgn, dist_eps);
  
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

Point2D Segment::LineIntersect(std::vector<double> n, double d2orgn, 
                               double dist_eps) const {
  SegLine::Position pos = PosWRT2Line(n, d2orgn, dist_eps);
  
  XMOF2D_ASSERT((pos == SegLine::Position::TOUCHES) || (pos == SegLine::Position::INTERSECTS),
         "Segment is above, below, or on the line!");
  
  Point2D pint;
  if (pos == SegLine::Position::TOUCHES)
    pint = (v[0].PosWRT2Line(n, d2orgn, dist_eps) == Point2DLine::ON) ? v[0] : v[1];
  else {
    std::vector<double> sdv = v[1] - v[0];
    double cpz = std::fabs(ddot(sdv, n));
    
    bool use_bisections = false;
    double t_eps = dist_eps/dnrm2(sdv);
    double t = std::fabs(-d2orgn + ddot(v[0].vec(), n)) / cpz;
    if ( (t <= t_eps) || (t >= 1.0 - t_eps) )
      use_bisections = true;
    else {
      dscal(t, sdv);
      pint = v[0] + sdv;
      if ( !is_interior(pint, dist_eps) || 
           (pint.PosWRT2Line(n, d2orgn, dist_eps) != Point2DLine::Position::ON))
        use_bisections = true;
    }

    if (use_bisections) {
      Segment bisected_seg = *this;
      bool int_pt_inside = true;
      do {
        Segment test_seg = Segment(bisected_seg[0], bisected_seg.middle());
        SegLine::Position cur_pos = test_seg.PosWRT2Line(n, d2orgn, dist_eps);
        if ((cur_pos == SegLine::Position::BELOW) || (cur_pos == SegLine::Position::ABOVE))
          bisected_seg = Segment(bisected_seg.middle(), bisected_seg[1]);
        else {
          bisected_seg = test_seg;
          if (cur_pos == SegLine::Position::TOUCHES)
            int_pt_inside = false;
        }
      } while (int_pt_inside);
      pint = bisected_seg[1];
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

bool is_ccw(const Segment& a, const Segment& b, double dist_eps) {
  std::vector<int> isamev = isame_vrt(a, b);
  if (isamev[0] == -1)
    return false;

  return is_ccw(a[(isamev[0] + 1)%2], a[isamev[0]], b[(isamev[1] + 1)%2], dist_eps);
}

bool is_ccw(const std::vector<Segment>& segs, double dist_eps) {
  int nseg = (int) segs.size();
  XMOF2D_ASSERT(nseg > 1, "Sequence should contain at least two segments");
  
  std::vector<Point2D> pts(nseg);
  std::vector<int> isamev;
  for (int iseg = 0; iseg < nseg - 1; iseg++) {
    isamev = isame_vrt(segs[iseg], segs[iseg + 1]);
    if (isamev[0] == -1)
      return false;

    pts[iseg] = segs[iseg][(isamev[0] + 1)%2];
  }
  pts[nseg - 1] = segs[nseg - 1][isamev[1]];
  if (pts[0] != segs[nseg - 1][(isamev[1] + 1)%2])
    pts.push_back(segs[nseg - 1][(isamev[1] + 1)%2]);

  return is_ccw(pts, dist_eps);
}

Point2D LineLineIntersect(const Segment& l1, const Segment& l2, double dist_eps) {
  double seg_sizes[] = {l1.size(), l2.size()};
  const Segment* s[2];
  if (seg_sizes[0] > seg_sizes[1]) { s[0] = &l2; s[1] = &l1; }
  else { s[0] = &l1; s[1] = &l2; }

  std::vector<double> s1_n = s[1]->normal();
  double s1_d2orgn = ddot((*s[1])[0].vec(), s1_n);
  SegLine::Position s0_pos = s[0]->PosWRT2Line(s1_n, s1_d2orgn, dist_eps);
  if (s0_pos == SegLine::Position::TOUCHES)
    return ((*s[0])[0].PosWRT2Line(s1_n, s1_d2orgn, dist_eps) == Point2DLine::ON) ? 
      (*s[0])[0] : (*s[0])[1];

  if (s0_pos == SegLine::Position::INTERSECTS)
    return s[0]->LineIntersect(s1_n, s1_d2orgn, dist_eps);

  double d[2] = { s[0]->dist_to_line((*s[1])[0]), s[0]->dist_to_line((*s[1])[1]) };
  double dd = d[0] - d[1];
  if (std::fabs(dd) < dist_eps) return BAD_POINT;

  Point2D pint;
  std::vector<double> seg_vec = (*s[1])[1] - (*s[1])[0];
  if (dd > 0) {
    dscal(d[1]/dd, seg_vec);
    pint = (*s[1])[1] + seg_vec;
  }
  else {
    dscal(d[0]/dd, seg_vec);
    pint = (*s[1])[0] + seg_vec;
  }

  if (s[0]->dist_to_line(pint) >= dist_eps) return BAD_POINT;

  return pint;
}

}
