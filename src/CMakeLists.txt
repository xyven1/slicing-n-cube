find_package(CGAL QUIET)

find_package (Threads REQUIRED)

add_executable(main main.cpp edge.cpp lp.cpp symmetry.cpp)

target_link_libraries(main CGAL::CGAL Threads::Threads)

add_executable(final final.c)

include_directories(../include ../extern)
