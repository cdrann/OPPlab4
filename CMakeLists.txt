cmake_minimum_required(VERSION 3.14)
project(OPPlab4)

set(CMAKE_CXX_STANDARD 14)

add_executable(OPPlab4 main_mpi.cpp Test.h Test.cpp JacobiMethod.h JacobiMethod.cpp)
# Require MPI for this project:

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(OPPlab4 ${MPI_LIBRARIES})