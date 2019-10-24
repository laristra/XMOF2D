/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#include "point2D.h"
#include "segment.h"
#include "tools.h"
#include "exception.h"
#include "simple_vector.h"

namespace XMOF2D {

Point2D::Point2D() : x(0), y(0) {
}

Point2D::Point2D(double _x, double _y) : x(_x), y(_y) {
}

Point2D::Point2D(const Point2D& p) : x(p.x), y(p.y) {
}

bool Point2D::operator==(const Point2D& p) const {
  return is_equal(x, p.x) && is_equal(y, p.y);
}

bool Point2D::operator!=(const Point2D& p) const {
  return !(*this == p);
}

const Point2D& Point2D::operator=(const Point2D& p) {
  if (this == &p) return *this;
  x = p.x; y = p.y;
  return *this;
}

Point2D Point2D::operator-() const {
  return Point2D(-x, -y);
}

Point2D Point2D::operator+(const Point2D& p) const {
  return Point2D(x + p.x, y + p.y);
}

Point2D Point2D::operator+(const std::vector<double>& v) const {
  XMOF2D_ASSERT_SIZE(v.size(), 2);
  return Point2D(x + v[0], y + v[1]);
}

std::vector<double> Point2D::operator-(const Point2D& p) const {
  return {x - p.x, y - p.y};
}

Point2D Point2D::operator*(double f) const {
  return Point2D(f*x, f*y);
}

Point2D Point2D::operator/(double f) const {
  XMOF2D_ASSERT(is_not_equal(f, 0.0), "Division by zero");
  return Point2D(x/f, y/f);
}

Point2D operator*(double f, const Point2D& p) {
  return p*f;
}

void Point2D::operator+=(const Point2D& p) {
  x += p.x;
  y += p.y;
}

void Point2D::operator/=(double f) {
  XMOF2D_ASSERT(is_not_equal(f, 0.0), "Division by zero");
  x /= f;
  y /= f;
}

std::vector<double> Point2D::vec() const {
  return {x, y};
}

double distance(const Point2D& p1, const Point2D& p2) {
  return sqrt(pow2(p1.x - p2.x) + pow2(p1.y - p2.y));
}

double vpz(const Point2D& a, const Point2D& b, const Point2D& c) {
  return ( (b.x - a.x)*(c.y - a.y) - (c.x -  a.x)*(b.y - a.y) );
}

int vpsign(const Point2D& a, const Point2D& b, const Point2D& c, 
           double dist_eps) {
  std::vector<double> vec[2] = {b - a, c - a};
  for (int iv = 0; iv < 2; iv++) {
    double vsize = dnrm2(vec[iv]);
    XMOF2D_ASSERT(vsize >= dist_eps, 
      "When computing the sign of a cross product, one of the vectors is of length " <<
      vsize << ", which is below the distance tolerance of " << dist_eps);
    dscal(1.0/vsize, vec[iv]);
  }

  double cpz = vec[0][0]*vec[1][1] - vec[1][0]*vec[0][1];

  if (is_equal(cpz, 0.0)) return 0;
  else return (std::signbit(cpz)) ? -1 : 1;
}

bool is_cw(const Point2D& a, const Point2D& b, const Point2D& c, double dist_eps) {
  return (vpsign(a, b, c, dist_eps) < 0);
}
  
bool is_ccw(const Point2D& a, const Point2D& b, const Point2D& c, double dist_eps) {
  return (vpsign(a, b, c, dist_eps) > 0);  
}

bool is_ccw(const std::vector<Point2D>& p, double dist_eps) {
  int np = (int) p.size();
  XMOF2D_ASSERT(np > 2, "Sequence should contain at least three points");

  std::vector<bool> is_hanging(np, false);
  int ifp = 0, isp;
  for (int i = 0; i < np - 1; i++) {
    isp = (i + 2)%np;
    if (ifp == isp) return false;
    if (Segment(p[ifp], p[isp]).contains(p[i + 1], dist_eps))
      is_hanging[i + 1] = true;
    else ifp = i + 1;
  }
  isp = 1;
  while (is_hanging[isp])
    isp = (isp + 1)%np;
  if (isp == ifp) return false;
  is_hanging[0] = Segment(p[ifp], p[isp]).contains(p[0], dist_eps);

  for (int i = 0; i < np; i++)  {
    int imp = (i + 1)%np;
    if (is_hanging[imp]) continue;

    int ifp = i, isp = (i + 2)%np;
    while (is_hanging[ifp])
      ifp = (np + ifp - 1)%np;
    while (is_hanging[isp])
      isp = (isp + 1)%np;

    if (!is_ccw(p[ifp], p[imp], p[isp], dist_eps))
      return false;
  }
  
  return true;
}
  
Point2DLine::Position Point2D::PosWRT2Line(std::vector<double> n, double d2orgn, 
                                           double dist_eps) const {
  std::vector<double> v2l = vec();
  daxpy(-d2orgn, n, v2l);
  double lv = ddot(v2l, n);
    
  Point2DLine::Position pos;
  if (std::fabs(lv) < dist_eps)
    pos = Point2DLine::Position::ON;
    else {
      if (std::signbit(lv))
        pos = Point2DLine::Position::BELOW;
      else
        pos = Point2DLine::Position::ABOVE;
  }
  
  return pos;
}

std::ostream& operator<<(std::ostream& os, const Point2D& p) {
  return os << "(" << p.x << ", " << p.y << ")";
}

}
