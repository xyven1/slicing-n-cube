find_package(CGAL QUIET)

find_package (Threads REQUIRED)

add_executable(main main.cpp edge.cpp lp.cpp symmetry.cpp)

target_link_libraries(main CGAL::CGAL Threads::Threads)

add_executable(slice_cube_low_weight slice_cube_low_weight.cpp edge.cpp symmetry.cpp)

add_executable(storage storage.cpp edge.cpp lp.cpp symmetry.cpp)

target_link_libraries(storage CGAL::CGAL)

add_executable(final final.c)

include_directories(../include ../extern)
