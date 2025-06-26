cmake_minimum_required (VERSION 3.8)

project(at_lat_ea)
set(CMAKE_CXX_FLAGS "-Wall -pedantic -O3 -Wextra")

add_library(mv_block OBJECT mv_block.cpp)
set_property(TARGET mv_block PROPERTY CXX_STANDARD 11)

add_library(lattice_ea OBJECT lattice_ea.cpp)
set_property(TARGET lattice_ea PROPERTY CXX_STANDARD 11)

add_library(lattice_ea_ops OBJECT lattice_ea_ops.cpp)
set_property(TARGET lattice_ea_ops PROPERTY CXX_STANDARD 11)

add_library(at_lat_ea $<TARGET_OBJECTS:mv_block> $<TARGET_OBJECTS:lattice_ea> $<TARGET_OBJECTS:lattice_ea_ops>)
