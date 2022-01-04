cmake_minimum_required(VERSION 3.20)

project(slicing_n_cube)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wfatal-errors")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wfatal-errors")

add_compile_definitions(N_CUBE_OUT_DIR="${N_CUBE_OUT_DIR}")

add_subdirectory(src)
