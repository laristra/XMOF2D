/*
 This file is part of the Ristra XMOF2D project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/XMOF2D/blob/master/LICENSE
 
 Created by Evgeny Kikinzon.
 Copyright © 2018, Triad National Security, LLC.
 All rights reserved.
*/

#ifndef exception_h
#define exception_h

#include <string>
#include <exception>
#include <iostream>
#include <sstream>

namespace XMOF2D {

class Exception : public std::exception {
protected:
  std::string msg;
public:
  Exception(const std::string& _msg = "Unspecified error occurred") throw() : msg(_msg) { }
  virtual ~Exception() throw() { }
  
  const char* what() const throw () { return msg.c_str(); }
  bool is(std::string expected_msg) { 
    std::size_t found = msg.find(expected_msg);
    return found != std::string::npos; 
  }
};

#define THROW_EXCEPTION(what) \
{ \
  std::ostringstream os; \
  os << "EXCEPTION (" << __FILE__ << ":" << __LINE__ << "): " << what << std::endl; \
  throw Exception(os.str()); \
}

#define XMOF2D_ASSERT(expression, what) \
{ \
  if (!(expression)) \
    THROW_EXCEPTION(what); \
}

#define XMOF2D_ASSERT_SIZE(got, expected) \
  XMOF2D_ASSERT(got == expected, "Wrong size: " << got << " (expected " << expected << ")")
#define XMOF2D_ASSERT_SIZE_LESS(got, expected) \
  XMOF2D_ASSERT(int(got) < int(expected), "Index is out of boundaries: " << got << \
    " (should be less than " << expected << ")")
}
#endif /* exception_h */
