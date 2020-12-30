cmake_minimum_required (VERSION 3.15)
project (flyingDiscSimulator)
enable_language (Fortran)

set(CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -fimplicit-none")
    set(bounds "-fbounds-check")
endif()

file(GLOB sources  ${CMAKE_SOURCE_DIR}/*.f90 ${CMAKE_SOURCE_DIR}/*.h)
add_executable(flyingDiscSimulator ${sources})
