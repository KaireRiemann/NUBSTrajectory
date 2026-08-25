#!/usr/bin/env bash
set -euo pipefail

# Rebuild and run the complete deterministic validation suite, then export the
# two paper-facing numerical tables. Usage:
#   scripts/run_validation.sh [build-directory] [benchmark-repetitions]
build_dir="${1:-build}"
runs="${2:-20}"
results_dir="${build_dir}/validation"

cmake -S . -B "${build_dir}" -DCMAKE_BUILD_TYPE=Release \
  -DNUBS_STRICT_NUMERICS=ON \
  -DNUBS_ENABLE_REFERENCE_FD=ON \
  -DNUBS_ENABLE_DENSE_ANALYTIC=ON \
  -DNUBS_ENABLE_BENCHMARKS=ON
cmake --build "${build_dir}" --parallel
ctest --test-dir "${build_dir}" --output-on-failure

mkdir -p "${results_dir}"
./bin/bench_nonuniform_construction "${runs}" \
  "${results_dir}/nonuniform_construction.csv"
./bin/bench_energy_gradient_validation "${runs}" \
  "${results_dir}/energy_gradient_validation.csv"
./bin/bench_nonuniform_robustness "${runs}" \
  "${results_dir}/nonuniform_robustness.csv"
python3 scripts/plot_nonuniform_robustness.py \
  "${results_dir}/nonuniform_robustness.csv" "${results_dir}/figures"
./bin/bench_gradient_propagation "${runs}" \
  "${results_dir}/gradient_propagation_benchmark.csv"
python3 scripts/plot_gradient_propagation_speed.py \
  "${results_dir}/gradient_propagation_benchmark.csv" \
  "${results_dir}/gradient_propagation_speed.svg"

echo "Validation outputs: ${results_dir}"
