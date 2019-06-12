/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#ifndef point_h
#define point_h

#include <cfloat>
#include <vector>
#include <iostream>

namespace XMOF2D {

namespace Point2DLine {
  enum Position {
    BELOW, ON, ABOVE
  };
}

struct Point2D {
  double x, y;
  
  Point2D();
  Point2D(double _x, double _y);
  Point2D(const Point2D& _p);
  
  bool operator==(const Point2D& p) const;
  bool operator!=(const Point2D& p) const;
  
  const Point2D& operator=(const Point2D& p);
  
  Point2D  operator+(const Point2D& p) const;
  Point2D  operator+(const std::vector<double>& v) const;
  Point2D  operator-() const;
  std::vector<double> operator-(const Point2D& p) const;
  Point2D  operator*(double f) const;
  Point2D  operator/(double f) const;
  
  void operator+=(const Point2D& p);
  void operator/=(double f);
  
  std::vector<double> vec() const;
  Point2DLine::Position PosWRT2Line(std::vector<double> n, double d2orgn, double eps = 1.0e-15);
};

Point2D operator*(double f, const Point2D& p);
std::ostream& operator<<(std::ostream& os, const Point2D& p);
double distance(const Point2D& p1, const Point2D& p2);

double vpz(const Point2D& a, const Point2D& b, const Point2D& c);
int vpsign(const Point2D& a, const Point2D& b, const Point2D& c, 
           double dist_eps = 1.0e-15, double ddot_eps = 1.0e-14);
int dpsign(const Point2D& a, const Point2D& b, const Point2D& c, 
           double dist_eps = 1.0e-15, double ddot_eps = 1.0e-14);
bool is_cw(const Point2D& a, const Point2D& b, const Point2D& c, 
           double dist_eps = 1.0e-15, double ddot_eps = 1.0e-14);
bool is_ccw(const Point2D& a, const Point2D& b, const Point2D& c, 
            double dist_eps = 1.0e-15, double ddot_eps = 1.0e-14);
bool is_ccw(const std::vector<Point2D>& p, double dist_eps = 1.0e-15, double ddot_eps = 1.0e-14);
bool on_same_line(const Point2D& a, const Point2D& b, const Point2D& c, 
                  double dist_eps = 1.0e-15, double ddot_eps = 1.0e-14);

Point2D LineLineIntersect(const std::vector<Point2D>& l1p, const std::vector<Point2D>& l2p, double eps = 1e-15);

static const Point2D BAD_POINT(DBL_MAX,DBL_MAX);

}
#endif /* point_h */
