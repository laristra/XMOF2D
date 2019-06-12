/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#ifndef segment_h
#define segment_h

#include "gobject.h"

namespace XMOF2D {

namespace SegLine {
  enum Position {
    BELOW, ON, TOUCHES, INTERSECTS, ABOVE
  };
}

class Segment : public GObject1D {
private:
  Point2D v[2];
public:
  Segment(const Point2D& a, const Point2D& b);
  Segment(const std::vector<Point2D>& points);
  
  const Segment& operator=(const Segment& seg);
  bool  operator==(const Segment& p) const;
  bool  operator!=(const Segment& p) const;
  Point2D operator[](int i) const;
  
  std::vector<Point2D> vrts() const { return {v[0], v[1]}; }
  
  int     npoints() const { return 2; }
  bool    contains(const Point2D& p, double eps = 1.0e-15) const;
  double  dist(const Point2D& p) const;
  
  double              size() const;
  Point2D               middle() const;
  std::vector<double> normal() const;
  std::vector<double> direction() const;
  
  SegLine::Position PosWRT2Line(std::vector<double> n, double d2orgn, double eps = 1.0e-15);
  Point2D LineIntersect(std::vector<double> n, double d2orgn, double denom_eps = 1.0e-14, double eps = 1.0e-15);
  
  friend std::ostream& operator<<(std::ostream& os, const Segment& s);
};

bool is_ccw(const Segment& a, const Segment& b, double dist_eps = 1.0e-15, double ddot_eps = 1.0e-14);
bool is_ccw(const std::vector<Segment>& segs, double dist_eps = 1.0e-15, double ddot_eps = 1.0e-14);
std::vector<int> isame_vrt(const Segment& a, const Segment& b);

}

#endif /* segment_h */
