find_package(CGAL QUIET)

find_package(Threads REQUIRED)

add_executable(slice_cube_low_weight slice_cube_low_weight.cpp)
target_link_libraries(slice_cube_low_weight Threads::Threads)

add_executable(slice_cube_degree_two slice_cube_degree_two.cpp)
# Required for GCC 8
target_link_libraries(slice_cube_degree_two stdc++fs)

add_executable(hyperplanes_low_weight hyperplanes_low_weight.cpp)
target_link_libraries(hyperplanes_low_weight CGAL::CGAL Threads::Threads)

add_executable(write_hyperplanes write_hyperplanes.cpp)
# Required for GCC 8
target_link_libraries(write_hyperplanes stdc++fs)

add_executable(storage storage.cpp)
target_link_libraries(storage CGAL::CGAL Threads::Threads)
# Required for GCC 8
target_link_libraries(storage stdc++fs)

add_executable(stats stats.cpp)
target_link_libraries(stats CGAL::CGAL Threads::Threads)

add_executable(slice_cube slice_cube.cpp)
target_link_libraries(slice_cube CGAL::CGAL Threads::Threads)

add_executable(final final.c)

include_directories(../include ../extern)
