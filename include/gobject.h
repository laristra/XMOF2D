/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Los Alamos National Security, LLC.
 All rights reserved.
*/

#ifndef gobject_h
#define gobject_h

#include "point2D.h"
#include <iostream>
#include <memory>

namespace XMOF2D {

class GObject {
public:
  virtual ~GObject() {};
  virtual double size() const = 0;
};

class GObject1D : public GObject {
};

class GObject2D : public GObject {
public:
  virtual std::unique_ptr<GObject1D> face(int i) const = 0;
  virtual int nfaces() const = 0;
};

}
#endif /* gobject_h */
