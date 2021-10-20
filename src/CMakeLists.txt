find_package(CGAL QUIET)

add_executable(main main.cpp complex.cpp edge.cpp lp.cpp sliceable_set.cpp symmetry.cpp)

target_link_libraries(main CGAL::CGAL)

add_executable(final final.c)

include_directories(../include ../extern)
