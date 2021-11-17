find_package(CGAL QUIET)

add_executable(main main.cpp edge.cpp lp.cpp symmetry.cpp)

target_link_libraries(main CGAL::CGAL)

add_executable(final final.c)

include_directories(../include ../extern)
