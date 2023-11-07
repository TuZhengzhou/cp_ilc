cmake_minimum_required(VERSION 3.0)

project(ZKP_{ALC-AC})

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_BUILD_TYPE "Release")
# set(CMAKE_BUILD_TYPE "Debug")

add_definitions(-DPROFILE_OP_COUNTS)

add_subdirectory(depends)
add_subdirectory(cp_ilc)
