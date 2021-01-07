cmake_minimum_required (VERSION 3.15)
enable_language (Fortran)

# if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#     set(dialect "-ffree-form -fimplicit-none -fcheck=all")
#     set(bounds "-fbounds-check -fcheck=all")
#     set(debug "-fcheck=all -Wall -g")
# endif()

# set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

# message("Fortran Flags" :${CMAKE_Fortran_FLAGS})
message("CMAKE_SOURCE_DIR" :${CMAKE_SOURCE_DIR})

file(GLOB sources  ${CMAKE_SOURCE_DIR}/extern/*.f90)
# add_executable(matinv66 ${sources})
add_library(matinv66 ${sources})

