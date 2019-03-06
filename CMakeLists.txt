#[[
This file is part of the Ristra XMOF2D project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/XMOF2D/blob/master/LICENSE

Created by Evgeny Kikinzon.
Copyright © 2018, Triad National Security, LLC. 
All rights reserved.
]]

# -*- mode: cmake -*-
#
# Top-Level CMake file for XMOF2D

# Require cmake 3.0.2 or higher
cmake_minimum_required(VERSION 3.0.2)
set(CMAKE_CTEST_COMMAND ${CMAKE_CTEST_COMMAND})

project (XMOF2D CXX)

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

if (PREFER_STATIC_LIBRARIES)
   set(CMAKE_FIND_LIBRARY_SUFFIXES ".a; .so")
endif()

##############################################################################
# XMOF2D LIBRARY BUILD CONFIGURATION
##############################################################################

option(ENABLE_FORTRAN "Build Fortran module and example" OFF)

# Default build type is Release (Optimized)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()

# include files and source files

include_directories(${PROJECT_SOURCE_DIR}/include)
set (xmof2Dfiles src/basemesh.cpp
                 src/cell.cpp
                 src/face.cpp
                 src/factory.cpp
                 src/meshbc.cpp
                 src/minimesh.cpp
                 src/node.cpp
                 src/point2D.cpp
                 src/segment.cpp
                 src/simple_convex.cpp
                 src/simple_vector.cpp
                 src/single_mmc_wrapper.cpp
                 src/xmof2D.cpp)

# Extra library suffix for Debug builds

if (CMAKE_BUILD_TYPE STREQUAL Debug)
  if (INSTALL_ADD_DEBUG_SUFFIX)
    set (CONFIG_SUFFIX -d)
  endif ()
  add_definitions(-DDEBUG)
endif ()

# Set library name

set (XMOF2DLIB xmof2D${CONFIG_SUFFIX})

add_library(${XMOF2DLIB} ${xmof2Dfiles})
message(STATUS "Building library " ${XMOF2DLIB})

# Default INSTALL Directory

if (NOT INSTALL_DIR)
  set (INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
endif ()

# where to install library
install(TARGETS ${XMOF2DLIB}
  EXPORT xmof2D
  ARCHIVE
  DESTINATION ${INSTALL_DIR}/lib
  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
)

# where to install include files

if (${INSTALL_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
else ()
file (GLOB FILES ${PROJECT_SOURCE_DIR}/include/*.h)
foreach (FILE ${FILES})
  set (hfiles ${hfiles} ${FILE})
endforeach ()
INSTALL(FILES 
        ${hfiles}
        DESTINATION ${INSTALL_DIR}/include)
endif ()

# Fortran module

if (ENABLE_FORTRAN)
  enable_language (Fortran)

  get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
  add_library(fxmof2D_single_mmc OBJECT fortran/xmof2d_single_mmc.f)
  set_target_properties(fxmof2D_single_mmc PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/fortran)

  # Install the Fortran module
  install(FILES ${CMAKE_SOURCE_DIR}/fortran/xmof2d_single_mmc.f 
          DESTINATION ${INSTALL_DIR}/fortran)
  install(FILES ${PROJECT_BINARY_DIR}/fortran/xmof2d_single_mmc.mod 
          DESTINATION ${INSTALL_DIR}/fortran)

  #Create Makefile_standalone
  set(makefile_standalone "${CMAKE_SOURCE_DIR}/fortran/src/Makefile" 
      CACHE STRING "Makefile for linking the library to a single .f file")
  FILE(WRITE ${makefile_standalone} "FORTRAN_APP_NAME=single_cell_fortran\n")
  FILE(APPEND ${makefile_standalone} "XMOF2D_INSTALL_DIR=${INSTALL_DIR}\n")
  FILE(APPEND ${makefile_standalone} "FC=${CMAKE_Fortran_COMPILER}\n")
  FILE(APPEND ${makefile_standalone} "\n\n\n################################################################################\n\n")
  FILE(APPEND ${makefile_standalone} "FC_FLAGS=-lstdc++\nMOD_DIR=\${XMOF2D_INSTALL_DIR}/fortran\nLIB_DIR=\${XMOF2D_INSTALL_DIR}/lib\n")
  FILE(APPEND ${makefile_standalone} "SRCS=\${FORTRAN_APP_NAME}.f\nOBJS=\$(SRCS:%.f=%.o)\n")
  FILE(APPEND ${makefile_standalone} "all:\n\t\${FC} -c \${SRCS} -I\${MOD_DIR}\n\t\${FC} -o \${FORTRAN_APP_NAME} \${OBJS} \${LIB_DIR}/libxmof2D.a \${FC_FLAGS}\n")
  FILE(APPEND ${makefile_standalone} "clean:\n\trm -f ./*.o\n\trm -f \${FORTRAN_APP_NAME}\n")
endif ()
###############################################################################
# Example
###############################################################################

add_subdirectory(examples)
if (ENABLE_FORTRAN)
  add_subdirectory(fortran)
endif()
############################################################################## 
# Write a configuration file from template replacing only variables enclosed
# by the @ sign. This will let other programs built on XMOF2D to discover how
# it was built
#############################################################################

configure_file(${PROJECT_SOURCE_DIR}/cmake/XMOF2DConfig.cmake.in 
               ${PROJECT_BINARY_DIR}/XMOF2DConfig.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/XMOF2DConfig.cmake 
        DESTINATION ${INSTALL_DIR}/share/cmake/)
install(FILES ${PROJECT_BINARY_DIR}/XMOF2DConfig.cmake
        DESTINATION ${INSTALL_DIR}/lib)
