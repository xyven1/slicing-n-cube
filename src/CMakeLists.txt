find_package(CGAL QUIET)

find_package (Threads REQUIRED)

add_executable(main main.cpp lp.cpp)

target_link_libraries(main CGAL::CGAL Threads::Threads)

add_executable(slice_cube_low_weight slice_cube_low_weight.cpp)

add_executable(storage storage.cpp lp.cpp)

target_link_libraries(storage CGAL::CGAL)

add_executable(stats stats.cpp lp.cpp)

target_link_libraries(stats CGAL::CGAL)

add_executable(slice_cube slice_cube.cpp lp.cpp)

target_link_libraries(slice_cube CGAL::CGAL)

add_executable(final final.c)

include_directories(../include ../extern)
