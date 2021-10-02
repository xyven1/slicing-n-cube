find_package(CGAL QUIET)

add_executable(main main.cpp complex.cpp lp.cpp symmetry.cpp)

target_link_libraries(main CGAL::CGAL)

include_directories(../include)
