find_package(CGAL QUIET)

find_package(Threads REQUIRED)

add_library(lp lp.cpp)
target_link_libraries(lp CGAL::CGAL)

add_executable(main main.cpp)
target_link_libraries(main lp Threads::Threads)

add_executable(slice_cube_low_weight slice_cube_low_weight.cpp)
target_link_libraries(slice_cube_low_weight Threads::Threads)

add_executable(hyperplanes_low_weight hyperplanes_low_weight.cpp)
target_link_libraries(hyperplanes_low_weight lp Threads::Threads)

add_executable(write_hyperplanes write_hyperplanes.cpp)
target_link_libraries(write_hyperplanes)

add_executable(storage storage.cpp)
target_link_libraries(storage lp Threads::Threads)

add_executable(stats stats.cpp)
target_link_libraries(stats lp Threads::Threads)

add_executable(slice_cube slice_cube.cpp)
target_link_libraries(slice_cube lp Threads::Threads)

add_executable(final final.c)

include_directories(../include ../extern)
