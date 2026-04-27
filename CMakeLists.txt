cmake_minimum_required(VERSION 3.10)

project(NUBSTrajectoryTests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

find_package(Eigen3 REQUIRED)

add_library(nubs_trajectory INTERFACE)
target_include_directories(nubs_trajectory INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(nubs_trajectory INTERFACE Eigen3::Eigen)

enable_testing()

function(add_nubs_test target)
    add_executable(${target} ${ARGN})
    target_link_libraries(${target} PRIVATE nubs_trajectory)
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    add_test(NAME ${target} COMMAND $<TARGET_FILE:${target}>)
endfunction()

add_nubs_test(test_basic_s2_d2 src/test_basic_s2_d2.cpp)
add_nubs_test(test_basic_s3_d3 src/test_basic_s3_d3.cpp)
add_nubs_test(test_basic_s4_d3 src/test_basic_s4_d3.cpp)
add_nubs_test(test_usage_example src/test_usage_example.cpp)

add_nubs_test(test_minco_s2_3d src/test_minco_s2_3d.cpp)
add_nubs_test(test_minco_s3_3d src/test_minco_s3_3d.cpp)
add_nubs_test(test_minco_s4_3d src/test_minco_s4_3d.cpp)

add_nubs_test(test_minco_gradient_s2_3d src/test_minco_gradient_s2_3d.cpp)
add_nubs_test(test_minco_gradient_s3_3d src/test_minco_gradient_s3_3d.cpp)
add_nubs_test(test_minco_gradient_s4_3d src/test_minco_gradient_s4_3d.cpp)

add_nubs_test(test_centered_gradient_s2_d2 src/test_centered_gradient_s2_d2.cpp)
add_nubs_test(test_centered_gradient_s3_d3 src/test_centered_gradient_s3_d3.cpp)
add_nubs_test(test_centered_gradient_s4_d3 src/test_centered_gradient_s4_d3.cpp)

add_nubs_test(test_boundary_grid src/test_boundary_grid.cpp)
add_nubs_test(test_generic_specialized_equivalence src/test_generic_specialized_equivalence.cpp)
add_nubs_test(test_external_gradient src/test_external_gradient.cpp)
add_nubs_test(test_local_vs_full_fd src/test_local_vs_full_fd.cpp)
add_nubs_test(test_small_duration src/test_small_duration.cpp)
add_nubs_test(bench_finite_diff_s3_d3 src/bench_finite_diff_s3_d3.cpp)
