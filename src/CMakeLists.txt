find_package(CGAL QUIET)

link_libraries(CGAL::CGAL)
# Required for GCC 8
link_libraries(stdc++fs)

add_executable(edge_cardinality edge_cardinality.cpp)
add_executable(slice_5_cube_c slice_5_cube.c)
add_executable(slice_5_cube_cpp slice_5_cube.cpp)
add_executable(slice_cube_degree_two slice_cube_degree_two.cpp)
add_executable(slice_cube_one_weight slice_cube_one_weight.cpp)
add_executable(smallest_equivalent_low_weight smallest_equivalent_low_weight.cpp)
add_executable(stats stats.cpp)
add_executable(storage storage.cpp)
add_executable(write_hyperplanes write_hyperplanes.cpp)

include_directories(../include ../extern)
