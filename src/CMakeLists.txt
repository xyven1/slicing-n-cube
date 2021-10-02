find_package(CGAL QUIET)
add_executable(main main.cpp)
target_link_libraries(main CGAL::CGAL)
