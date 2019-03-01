cmake_minimum_required (VERSION 2.6)

project(lea)
set(CMAKE_CXX_FLAGS "-Wall -pedantic -O3 -Wextra")
add_subdirectory(at_lat_ea)

add_executable(lea lea.cpp)
set_property(TARGET lea PROPERTY CXX_STANDARD 14)

add_executable(lattice_ea_tests lattice_ea_tests.cpp)
set_property(TARGET lattice_ea_tests PROPERTY CXX_STANDARD 11)
#target_link_libraries(lattice_ea_tests lattice_ea_ops lattice_ea mv_block)
target_link_libraries(lattice_ea_tests at_lat_ea)

add_executable(mv_block_tests mv_block_tests.cpp)
set_property(TARGET mv_block_tests PROPERTY CXX_STANDARD 11)
#target_link_libraries(mv_block_tests $<TARGET_OBJECTS:at_lat_ea/mv_block>)
target_link_libraries(mv_block_tests $<TARGET_OBJECTS:mv_block>)


