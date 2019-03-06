/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#include "point2D.h"
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

double vpz(const Point2D& a, const Point2D& b, const Point2D& c)
{
  return ( (b.x - a.x)*(c.y - a.y) - (c.x -  a.x)*(b.y - a.y) );
}

bool is_cw(const Point2D& a, const Point2D& b, const Point2D& c, double ddot_eps) {
  double pts_vpz = vpz(a, b, c);
  
  return (pts_vpz < ddot_eps);
}

bool on_same_line(const Point2D& a, const Point2D& b, const Point2D& c, double ddot_eps) {
  double pts_vpz = vpz(a, b, c);
  
  return (std::fabs(pts_vpz) < ddot_eps);
}
  
bool is_ccw(const Point2D& a, const Point2D& b, const Point2D& c, double ddot_eps) {
  double pts_vpz = vpz(a, b, c);
  
  return (pts_vpz > -ddot_eps);
}

bool is_ccw(const std::vector<Point2D>& p, double ddot_eps) {
  int np = (int) p.size();
  XMOF2D_ASSERT(np > 2, "Sequence should contain at least three points");
  
  for (int i = 0; i < np; i++)  {
    if (!is_ccw(p[i], p[(i + 1)%np], p[(i + 2)%np], ddot_eps))
      return false;
  }
  
  return true;
}
  
Point2DLine::Position Point2D::PosWRT2Line(std::vector<double> n, double d2orgn, double eps) {
  std::vector<double> v2l = vec();
  daxpy(-d2orgn, n, v2l);
  double lv = ddot(v2l, n);
  
  
  Point2DLine::Position pos;
  if (std::fabs(lv) < eps)
  pos = Point2DLine::Position::ON;
  else {
    if (std::signbit(lv))
    pos = Point2DLine::Position::BELOW;
    else
    pos = Point2DLine::Position::ABOVE;
  }
  
  return pos;
}

Point2D LineLineIntersect(const std::vector<Point2D>& l1p, const std::vector<Point2D>& l2p, double eps) {
  XMOF2D_ASSERT_SIZE(l1p.size(), 2);
  XMOF2D_ASSERT_SIZE(l2p.size(), 2);
  double denom = (l1p[0].x - l1p[1].x)*(l2p[0].y - l2p[1].y) -
                 (l1p[0].y - l1p[1].y)*(l2p[0].x - l2p[1].x);
  if (std::fabs(denom) < eps)
    return BAD_POINT;
  
  Point2D pint;
  pint.x = (l1p[0].x*l1p[1].y - l1p[0].y*l1p[1].x)*(l2p[0].x - l2p[1].x) -
           (l2p[0].x*l2p[1].y - l2p[0].y*l2p[1].x)*(l1p[0].x - l1p[1].x);
  pint.y = (l1p[0].x*l1p[1].y - l1p[0].y*l1p[1].x)*(l2p[0].y - l2p[1].y) -
           (l2p[0].x*l2p[1].y - l2p[0].y*l2p[1].x)*(l1p[0].y - l1p[1].y);
  
  return pint/denom;
}

std::ostream& operator<<(std::ostream& os, const Point2D& p) {
  return os << "(" << p.x << ", " << p.y << ")";
}

}
