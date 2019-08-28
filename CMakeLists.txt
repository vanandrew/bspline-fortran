#
# Basic CMake configuration file for bspline-fortran
#

cmake_minimum_required (VERSION 2.8)

set(PROJECT_DESCRIPTION "Multidimensional B-Spline Interpolation of Data on a Regular Grid")
set(PROJECT_URL "https://github.com/jacobwilliams/bspline-fortran")

project (bspline-fortran NONE)
enable_language (Fortran)

add_library(${PROJECT_NAME} STATIC
    src/bspline_interp.f90
    src/bspline_kinds_module.f90
    src/bspline_module.f90
    src/bspline_sub_module.f90)
target_compile_options(${PROJECT_NAME} PRIVATE -fopenmp -ftree-vectorize -O3 -march=native)
