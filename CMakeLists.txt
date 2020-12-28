cmake_minimum_required (VERSION 3.0)
project (flyingDiscSimulator)
enable_language (Fortran)

file(GLOB_RECURSE sources  src/*.f90 src/*.h)

add_executable(flyingDiscSimulator ${sources})
