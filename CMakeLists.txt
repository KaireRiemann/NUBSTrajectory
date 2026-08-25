cmake_minimum_required(VERSION 3.10)

project(NUBSTrajectoryTests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

find_package(Eigen3 REQUIRED)

option(NUBS_ENABLE_REFERENCE_FD
       "Build finite-difference validation paths and tests" ON)
option(NUBS_ENABLE_DENSE_ANALYTIC
       "Build dense analytic timing-gradient validation paths" ON)
option(NUBS_ENABLE_BENCHMARKS
       "Build reproducible benchmark executables" ON)
option(NUBS_ENABLE_PLANNER_TESTS
       "Enable external planner integration tests when supplied by a parent project" OFF)
option(NUBS_STRICT_NUMERICS
       "Treat compiler warnings as errors for this header-only library's tests" OFF)

add_library(nubs_trajectory INTERFACE)
target_include_directories(nubs_trajectory INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(nubs_trajectory INTERFACE Eigen3::Eigen)
target_compile_definitions(nubs_trajectory INTERFACE
    $<$<BOOL:${NUBS_ENABLE_REFERENCE_FD}>:NUBS_ENABLE_REFERENCE_FD=1>
    $<$<BOOL:${NUBS_ENABLE_DENSE_ANALYTIC}>:NUBS_ENABLE_DENSE_ANALYTIC=1>
)

enable_testing()

function(add_nubs_test target)
    add_executable(${target} ${ARGN})
    target_link_libraries(${target} PRIVATE nubs_trajectory)
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    if(NUBS_STRICT_NUMERICS)
        target_compile_options(${target} PRIVATE -Werror)
    endif()
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
add_nubs_test(test_fixed_ratio_total_duration src/test_fixed_ratio_total_duration.cpp)
add_nubs_test(test_uniform_time src/test_uniform_time.cpp)
add_nubs_test(test_large_scale_uniform_compare src/test_large_scale_uniform_compare.cpp)
add_nubs_test(test_small_duration src/test_small_duration.cpp)
add_nubs_test(bench_finite_diff_s3_d3 src/bench_finite_diff_s3_d3.cpp)
add_nubs_test(test_dual_scalar src/test_dual_scalar.cpp)
add_nubs_test(test_timing_stencil src/test_timing_stencil.cpp)
add_nubs_test(test_nonuniform_robustness src/test_nonuniform_robustness.cpp)

if(NUBS_ENABLE_DENSE_ANALYTIC AND NUBS_ENABLE_REFERENCE_FD)
    add_nubs_test(test_local_ad_gradient src/test_local_ad_gradient.cpp)
endif()

if(NUBS_ENABLE_BENCHMARKS)
    add_executable(bench_large_scale_construction src/bench_large_scale_construction.cpp)
    target_link_libraries(bench_large_scale_construction PRIVATE nubs_trajectory)
    target_compile_options(bench_large_scale_construction PRIVATE -Wall -Wextra -Wpedantic)

    add_executable(bench_nonuniform_construction src/bench_nonuniform_construction.cpp)
    target_link_libraries(bench_nonuniform_construction PRIVATE nubs_trajectory)
    target_compile_options(bench_nonuniform_construction PRIVATE -Wall -Wextra -Wpedantic)

    add_executable(bench_energy_gradient_validation src/bench_energy_gradient_validation.cpp)
    target_link_libraries(bench_energy_gradient_validation PRIVATE nubs_trajectory)
    target_compile_options(bench_energy_gradient_validation PRIVATE -Wall -Wextra -Wpedantic)

    add_executable(bench_nonuniform_robustness src/bench_nonuniform_robustness.cpp)
    target_link_libraries(bench_nonuniform_robustness PRIVATE nubs_trajectory)
    target_compile_options(bench_nonuniform_robustness PRIVATE -Wall -Wextra -Wpedantic)

    add_executable(export_equivalence_figure_data src/export_equivalence_figure_data.cpp)
    target_link_libraries(export_equivalence_figure_data PRIVATE nubs_trajectory)
    target_compile_options(export_equivalence_figure_data PRIVATE -Wall -Wextra -Wpedantic)
endif()
