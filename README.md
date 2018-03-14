# XMOF2D

The XMOF2D library provides an implementation of the X-MOF interface reconstruction method on meshes with convex polygonal cells. The X-MOF method is based on the Moment-of-Fluid method,

Vadim Dyadechko, Mikhail Shashkov,
*Reconstruction of multi-material interfaces from moment data,*
Journal of Computational Physics 227(11), Pages 5361-5384, 2008,

and establishes the topology of meshes inside multi-material cells in a robust and efficient way.
This code operates on the so-called base mesh, which can contain single-material and multi-material cells. Given the material data for each cell, which includes reference volume fractions and centroids of contained  materials, interface reconstruction can be performed in multi-material cells, resulting in so-called minimeshes: local meshes representing the material distribution inside multi-material cells. Each minimesh has full local topology and, additionally, parentage information, i.e. for every mesh entity in a minimesh, the index of respective parent entity in the base mesh is known.
An optional Fortran module that provides a Fortran-friendly interface for interface reconstruction on a single multi-material cell is also provided.

## Getting Started

To obtain a copy of XMOF2D from GitHub, simply clone:
```sh
git clone https://github.com/laristra/XMOF2D.git
```

## Prerequisites

XMOF2D uses standard C++11 features, therefore a compatible compiler is required. The build system requires CMake version 2.8.8 or later.
Compiling the Fortran module requires a Fortran compiler compatible with the 2008 standard.

## Installing

The recommended way to build the library is by executing the `config/do-configure` script from the XMOF2D directory after specifying the desired configuration parameters.
Some important parameters are described below.

- `XMOF2D_DIR`: full path to the XMOF2D directory in which `config/do-configure` is located.
- `CC`, `CXX`, and `FC`: C, C++, and Fortran compilers to be used.
- `INSTALL_DIR`: directory where the library will be installed, by default the standard system directories will be used (e.g. `/usr/local` on OS X).
- `BUILD_TYPE`: controls the optimization flags, the default option is Release.

## Installation using modules

On HPC machines or workstations with modules,
```sh
module list
```
provides the list of currently loaded modules and
```sh
module avail
```
provides the list of modules that can be loaded.
Before running `config/do-configure`, the modules for compilers and cmake need to be loaded, for example
```sh
module load gcc/6.1.0
module load cmake/3.7.2
```
Then, if GNU compiler is to be used as above, ``which gcc``, ``which g++``, and ``which gfortran`` can be specified as `CC`, `CXX`, and `FC` parameters in `config/do-configure`.

Please note that the same module for the compiler need to be loaded when running applications linked to XMOF2D, otherwise incorrect shared libraries might be used.

# License

This project is licensed under a modified 3-clause BSD license - see
the [LICENSE](https://github.com/larista/XMOF2D/blob/master/LICENSE)
file for details.

# Release

This software has been approved for open source release and has been
assigned **LA-CC-18-010**.

## Examples

Along with the XMOF2D library, examples are built and installed in the `examples` subdirectory of `XMOF2D_DIR`. This directory includes two C++ examples:
- `single_cell`: this application demonstrates the use of the wrapper to perform interface reconstruction on a single multi-material cell. It invokes most of the available wrapper functions and outputs the resulting minimesh information. Additionally, the gnuplot files for plotting the base cell and the respective minimesh are written to the directory from which `example/single_cell` is executed.
To plot the base cell in gnuplot, use `load "gnuplot_script_mmc.gnu"` and to plot the corresponding minimesh, use `load "gnuplot_script_minimesh.gnu"`
The functions provided by the wrapper can be found in `INSTALL_DIR/include/single_mmc_wrapper.h`
- `test_line_interface`: this application demonstrates interface reconstruction on a unit square containing two materials separated by a linear material interface. The base mesh is a rectangular mesh of dimensions specified in the command line. The material data is loaded from an external file, sample data files for square meshes of different resolutions can be found in `examples/test_line_interface_data`. If the sample data file for, say, 64x64 mesh is to be used, then this example can be executed from `XMOF2D_DIR` with
```sh
examples/test_line_interface examples/test_line_interface_data/line_64x64_matdata.dat 64 64
```
Note that the provided material data was obtained using sampling with the similar number of points regardless of the mesh resolution. It illustrates the dependence of accuracy of interface reconstruction on the combination of sampling and mesh resolutions.
The output is the max Hausdorff distance between reconstructed and reference material interfaces over all multi-material cells. While theoretically zero, in practice it depends on tolerances and sampling resolution. For provided sample data, the resulting value should be on the order of 1e-6 or less.

## Fortran module

If `ENABLE_FORTRAN` parameter is set to `ON` in `config/do-configure`,
```sh
-D ENABLE_FORTRAN:BOOL=ON \
```
the `xmof2d_single_mmc` module will be built and installed to `INSTALL_DIR/fortran`. This module provides a Fortran interface to the wrapper in `INSTALL_DIR/include/single_mmc_wrapper.h`.
This wrapper exposes a list of procedures that allow to work with a single multi-material cell without any knowledge of XMOF2D's internal data structures.
The `single_cell_fortran` example will be built along with the XMOF2D library and placed in the `fortran` subdirectory. This example is a Fortran variant of the C++ example in `examples/single_cell`, which was described above.

Note that during the compilation of Fortran applications that use `xmof2d_single_mmc` module, `INSTALL_DIR/fortran` should be provided to the compiler as the include directory.

To illustrate linking of the XMOF2D library to a standalone Fortran application, a makefile template is given in `fortran/src/Makefile`. It allows to build an application from a single Fortran file that uses the `xmof2d_single_mmc` module.

The makefile in `fortran/src/Makefile` can be placed in any directory with the Fortran code that uses the Fortran module.
The parameters in the makefile that can be modified by the user are:

- `XMOF2D_INSTALL_DIR`: Path to the XMOF2D installation directory, usually `install-Debug` or `install-Release` in the XMOF2D directory.

- `FORTRAN_APP_NAME`: Name of the file with the Fortran code that uses XMOF2D module. To build an application from `single_cell_fortran.f`, it should be `FORTRAN_APP_NAME=single_cell_fortran`. The name of the executable will be `$FORTRAN_APP_NAME`.

- `FC`: Fortran compiler, ``which gfortran`` will use the default GNU compiler (or a GNU compiler from a loaded module).

After these parameters are specified, typing `make` will build an executable. `make clean` can be used to remove the executable and the object file.

When this makefile is created, `XMOF2D_INSTALL_DIR` and `FC` parameters are set to be the same as in `config/do-configure` script used to build the library.

Note that it is advised to use the same Fortran compiler that was used to build the XMOF2D library. If modules are used, it might be necessary to have the compiler module loaded when running the application: otherwise, incorrect shared libraries might be used.

As provided, the makefile can be used to build a copy of `fortran/single_cell_fortran` in `fortran/src` directory without rebuilding the XMOF2D library.
