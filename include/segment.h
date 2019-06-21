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
  const Point2D& operator[](int i) const;
  
  std::vector<Point2D> vrts() const { return {v[0], v[1]}; }
  
  int     npoints() const { return 2; }
  bool    contains(const Point2D& p, double dist_eps) const;
  bool    is_interior(const Point2D& p, double dist_eps) const;
  double  dist_to_seg(const Point2D& p) const;
  double  dist_to_line(const Point2D& p) const;
  
  double              size() const;
  Point2D             middle() const;
  std::vector<double> normal() const;
  std::vector<double> direction() const;
  
  SegLine::Position PosWRT2Line(std::vector<double> n, double d2orgn, double dist_eps) const;
  Point2D LineIntersect(std::vector<double> n, double d2orgn, double dist_eps) const;
  
  friend std::ostream& operator<<(std::ostream& os, const Segment& s);
};

bool is_ccw(const Segment& a, const Segment& b, double dist_eps);
bool is_ccw(const std::vector<Segment>& segs, double dist_eps);
std::vector<int> isame_vrt(const Segment& a, const Segment& b);

Point2D LineLineIntersect(const Segment& l1, const Segment& l2, double dist_eps);

}

#endif /* segment_h */
