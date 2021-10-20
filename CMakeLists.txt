cmake_minimum_required(VERSION 3.20)

project(slicing_n_cube)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wfatal-errors")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wfatal-errors")

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wbad-function-cast -Wc++-compat -Wjump-misses-init -Wmissing-prototypes -Wnested-externs -Wold-style-definition -Wstrict-prototypes -Wtraditional -Wtraditional-conversion -Wunsuffixed-float-constants")

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wabi -Wabi-tag -Waggregate-return -Walloc-zero -Walloca -Wc++11-compat -Wc++14-compat -Wc++17-compat -Wcast-align -Wcast-align=strict -Wcast-qual -Wclass-memaccess -Wconditionally-supported -Wctor-dtor-privacy -Wdate-time -Wdelete-non-virtual-dtor -Wdeprecated-copy -Wdeprecated-copy-dtor -Wdisabled-optimization -Wdouble-promotion -Wduplicated-branches -Wduplicated-cond -Weffc++ -Wextra-semi -Wfloat-conversion -Wfloat-equal -Wformat-nonliteral -Wformat-security -Wformat-signedness -Wformat-y2k -Winit-self -Winline -Winvalid-pch -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wmultichar -Wmultiple-inheritance -Wnamespaces -Wnoexcept -Wnoexcept-type -Wnon-virtual-dtor -Wnull-dereference -Wold-style-cast -Woverloaded-virtual -Wpacked -Wpadded -Wpessimizing-move -Wredundant-decls -Wredundant-move -Wregister -Wreorder -Wsign-conversion -Wsign-promo -Wstack-protector -Wstrict-null-sentinel -Wsuggest-attribute=cold -Wsuggest-attribute=const -Wsuggest-attribute=format -Wsuggest-attribute=malloc -Wsuggest-attribute=noreturn -Wsuggest-attribute=pure -Wsuggest-final-methods -Wsuggest-final-types -Wsuggest-override -Wswitch-default -Wswitch-enum -Wsynth -Wsystem-headers -Wtemplates -Wtrampolines -Wundef -Wunused-macros -Wuseless-cast -Wvector-operation-performance -Wvirtual-inheritance -Wwrite-strings -Wzero-as-null-pointer-constant")

add_compile_definitions(NCUBE_DIR="${CMAKE_SOURCE_DIR}/../Results/")

add_subdirectory(src)
