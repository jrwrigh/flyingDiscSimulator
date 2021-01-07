cmake_minimum_required (VERSION 3.15)
project (flyingDiscSimulator)
enable_language (Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -fimplicit-none -fcheck=all")
    set(bounds "-fbounds-check -fcheck=all")
    set(debug "-fcheck=all -Wall -g")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

message("Fortran Flags" :${CMAKE_Fortran_FLAGS})

set(CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

add_subdirectory(extern)

file(GLOB sources  ${CMAKE_SOURCE_DIR}/*.f90 ${CMAKE_SOURCE_DIR}/*.h)
add_executable(flyingDiscSimulator ${sources})

target_link_libraries(flyingDiscSimulator matinv66)

# add_custom_command(TARGET flyingDiscSimulator
#                    PRE_BUILD
#                    COMMAND touch testfile && echo '!cmake command worked'>>${CMAKE_SOURCE_DIR}/main.f90)
