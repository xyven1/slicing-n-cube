find_package(CGAL QUIET)

find_package(Threads REQUIRED)

link_libraries(CGAL::CGAL Threads::Threads)
# Required for GCC 8
link_libraries(stdc++fs)

add_executable(slice_cube_low_weight slice_cube_low_weight.cpp)

add_executable(slice_cube_degree_two slice_cube_degree_two.cpp)

add_executable(hyperplanes_low_weight hyperplanes_low_weight.cpp)

add_executable(write_hyperplanes write_hyperplanes.cpp)

add_executable(storage storage.cpp)

add_executable(stats stats.cpp)

add_executable(slice_cube slice_cube.cpp)

add_executable(final final.c)

include_directories(../include ../extern)
