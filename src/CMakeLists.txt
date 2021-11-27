find_package(CGAL QUIET)

find_package(Threads REQUIRED)

add_library(lp lp.cpp)
target_link_libraries(lp CGAL::CGAL)

add_executable(main main.cpp)
target_link_libraries(main lp)

add_executable(slice_cube_low_weight slice_cube_low_weight.cpp)

add_executable(storage storage.cpp)
target_link_libraries(storage lp)

add_executable(stats stats.cpp)
target_link_libraries(stats lp)

add_executable(slice_cube slice_cube.cpp)
target_link_libraries(slice_cube lp)

add_executable(final final.c)

include_directories(../include ../extern)
