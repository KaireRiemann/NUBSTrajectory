#include "NUBSTrajectory.hpp"
#include "tools/minco_adapter.hpp"
#include "tools/optimization_test_cases.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr int Dim = 3;

struct BenchmarkRow
{
    int order = 0;
    int piece_num = 0;
    int runs = 0;
    int point_variables = 0;
    int time_variables = 0;
    int nubs_control_variables = 0;
    int minco_coefficient_variables = 0;
    double nubs_direct_gradient_us = 0.0;
    double minco_direct_gradient_us = 0.0;
    double minco_vs_nubs_direct_speedup = 0.0;
    double nubs_full_propagation_us = 0.0;
    double minco_full_propagation_us = 0.0;
    double minco_vs_nubs_full_speedup = 0.0;
    double cost_abs_error = 0.0;
    double point_gradient_max_abs_error = 0.0;
    double time_gradient_max_abs_error = 0.0;
};

inline void keepAlive(const double value)
{
    volatile double sink = value;
    (void)sink;
}

inline void printRow(const BenchmarkRow &row)
{
    std::cout << "S=" << row.order << ", M=" << row.piece_num
              << ", direct: nubs=" << std::fixed << std::setprecision(3)
              << row.nubs_direct_gradient_us << " us, minco="
              << row.minco_direct_gradient_us << " us, minco/nubs="
              << row.minco_vs_nubs_direct_speedup << "x"
              << "; full: nubs=" << row.nubs_full_propagation_us
              << " us, minco=" << row.minco_full_propagation_us
              << " us, minco/nubs=" << row.minco_vs_nubs_full_speedup << "x"
              << std::scientific << std::setprecision(3)
              << "; point_err=" << row.point_gradient_max_abs_error
              << ", time_err=" << row.time_gradient_max_abs_error
              << std::endl;
}

inline void writeCsv(const std::string &path,
                     const std::vector<BenchmarkRow> &rows)
{
    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error("failed to open CSV output: " + path);
    }
    out << "order,piece_count,runs,point_variables,time_variables,"
           "nubs_control_variables,minco_coefficient_variables,"
           "nubs_direct_gradient_us,minco_direct_gradient_us,"
           "minco_vs_nubs_direct_speedup,nubs_full_propagation_us,"
           "minco_full_propagation_us,minco_vs_nubs_full_speedup,"
           "cost_abs_error,point_gradient_max_abs_error,"
           "time_gradient_max_abs_error\n";
    out << std::setprecision(17);
    for (const BenchmarkRow &row : rows)
    {
        out << row.order << ',' << row.piece_num << ',' << row.runs << ','
            << row.point_variables << ',' << row.time_variables << ','
            << row.nubs_control_variables << ','
            << row.minco_coefficient_variables << ','
            << row.nubs_direct_gradient_us << ','
            << row.minco_direct_gradient_us << ','
            << row.minco_vs_nubs_direct_speedup << ','
            << row.nubs_full_propagation_us << ','
            << row.minco_full_propagation_us << ','
            << row.minco_vs_nubs_full_speedup << ','
            << row.cost_abs_error << ','
            << row.point_gradient_max_abs_error << ','
            << row.time_gradient_max_abs_error << '\n';
    }
}

template <int S>
BenchmarkRow runOne(const int piece_num, const int runs)
{
    const auto data = nubs_test::makeRandomProblem<Dim>(
        S, piece_num, static_cast<unsigned int>(970000 + 1000 * S + piece_num));

    nubs::NUBSTrajectoryT<Dim, S> nubs_trajectory;
    Eigen::MatrixXd nubs_controls;
    nubs_trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                             data.durations, nubs_controls);

    nubs_test::SplineMincoAdapter<Dim> minco_trajectory;
    minco_trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                              data.durations, S);

    // The direct terms are deliberately timed in each implementation's native
    // internal coordinates (B-spline controls versus polynomial coefficients).
    // Only the fully propagated waypoint/time gradients are coordinate
    // invariant, so those are the numerical cross-method validation target.
    double nubs_direct_cost = 0.0;
    Eigen::MatrixXd nubs_direct_coefficients;
    Eigen::VectorXd nubs_direct_times;
    nubs_trajectory.getEnergyPartialGradLocalAD(
        nubs_direct_cost, nubs_direct_coefficients, nubs_direct_times);

    double minco_direct_cost = minco_trajectory.getEnergy();
    Eigen::MatrixXd minco_direct_coefficients;
    Eigen::VectorXd minco_direct_times;
    minco_trajectory.getEnergyPartialGradByCoeffs(minco_direct_coefficients);
    minco_trajectory.getEnergyPartialGradByTimes(minco_direct_times);

    double nubs_full_cost = 0.0;
    Eigen::MatrixXd nubs_point_gradient;
    Eigen::VectorXd nubs_time_gradient;
    nubs_trajectory.getEnergyAndGrad(nubs_full_cost, nubs_point_gradient,
                                     nubs_time_gradient);

    double minco_full_cost = 0.0;
    Eigen::MatrixXd minco_full_coefficients;
    Eigen::Matrix<double, Eigen::Dynamic, Dim> minco_point_gradient;
    Eigen::VectorXd minco_time_gradient;
    Eigen::VectorXd minco_partial_time_gradient;
    minco_trajectory.getEnergyAndGrad(
        minco_full_cost, minco_full_coefficients, minco_point_gradient,
        minco_time_gradient, minco_partial_time_gradient);

    nubs_test::requireNear(nubs_direct_cost, minco_direct_cost,
                           2.0e-6, 2.0e-8, "direct energy NUBS vs MINCO");
    nubs_test::requireNear(nubs_full_cost, minco_full_cost,
                           2.0e-6, 2.0e-8, "full energy NUBS vs MINCO");
    nubs_test::requireMatrixNear(nubs_point_gradient, minco_point_gradient,
                                 5.0e-4, 2.0e-7,
                                 "point gradient NUBS vs MINCO");
    nubs_test::requireMatrixNear(nubs_time_gradient, minco_time_gradient,
                                 5.0e-3, 5.0e-6,
                                 "time gradient NUBS vs MINCO");

    const auto nubs_direct_timing = nubs_test::measureRepeated(
        [&]()
        {
            double cost = 0.0;
            Eigen::MatrixXd coefficient_gradient;
            Eigen::VectorXd direct_time_gradient;
            nubs_trajectory.getEnergyPartialGradLocalAD(
                cost, coefficient_gradient, direct_time_gradient);
            keepAlive(cost + coefficient_gradient.squaredNorm() +
                      direct_time_gradient.squaredNorm());
        },
        runs);
    const auto minco_direct_timing = nubs_test::measureRepeated(
        [&]()
        {
            const double cost = minco_trajectory.getEnergy();
            Eigen::MatrixXd coefficient_gradient;
            Eigen::VectorXd direct_time_gradient;
            minco_trajectory.getEnergyPartialGradByCoeffs(coefficient_gradient);
            minco_trajectory.getEnergyPartialGradByTimes(direct_time_gradient);
            keepAlive(cost + coefficient_gradient.squaredNorm() +
                      direct_time_gradient.squaredNorm());
        },
        runs);
    const auto nubs_full_timing = nubs_test::measureRepeated(
        [&]()
        {
            double cost = 0.0;
            Eigen::MatrixXd point_gradient;
            Eigen::VectorXd time_gradient;
            nubs_trajectory.getEnergyAndGrad(cost, point_gradient, time_gradient);
            keepAlive(cost + point_gradient.squaredNorm() +
                      time_gradient.squaredNorm());
        },
        runs);
    const auto minco_full_timing = nubs_test::measureRepeated(
        [&]()
        {
            double cost = 0.0;
            Eigen::MatrixXd coefficient_gradient;
            Eigen::Matrix<double, Eigen::Dynamic, Dim> point_gradient;
            Eigen::VectorXd time_gradient;
            Eigen::VectorXd partial_time_gradient;
            minco_trajectory.getEnergyAndGrad(
                cost, coefficient_gradient, point_gradient, time_gradient,
                partial_time_gradient);
            keepAlive(cost + point_gradient.squaredNorm() +
                      time_gradient.squaredNorm());
        },
        runs);

    BenchmarkRow row;
    row.order = S;
    row.piece_num = piece_num;
    row.runs = runs;
    row.point_variables = std::max(0, piece_num - 1) * Dim;
    row.time_variables = piece_num;
    row.nubs_control_variables = nubs_controls.size();
    row.minco_coefficient_variables = minco_direct_coefficients.size();
    row.nubs_direct_gradient_us = nubs_direct_timing.avg_us;
    row.minco_direct_gradient_us = minco_direct_timing.avg_us;
    row.minco_vs_nubs_direct_speedup =
        nubs_test::speedupRatio(minco_direct_timing, nubs_direct_timing);
    row.nubs_full_propagation_us = nubs_full_timing.avg_us;
    row.minco_full_propagation_us = minco_full_timing.avg_us;
    row.minco_vs_nubs_full_speedup =
        nubs_test::speedupRatio(minco_full_timing, nubs_full_timing);
    row.cost_abs_error = std::abs(nubs_full_cost - minco_full_cost);
    row.point_gradient_max_abs_error =
        nubs_test::maxAbsCoeffDiff(nubs_point_gradient, minco_point_gradient);
    row.time_gradient_max_abs_error =
        nubs_test::maxAbsCoeffDiff(nubs_time_gradient, minco_time_gradient);
    printRow(row);
    return row;
}

} // namespace

int main(int argc, char **argv)
{
    const int runs = argc > 1 ? std::max(1, std::atoi(argv[1])) : 500;
    const std::string csv_path = argc > 2
                                     ? argv[2]
                                     : "gradient_propagation_benchmark.csv";
    try
    {
        const std::vector<int> piece_counts = {2, 4, 8, 16, 32, 64};
        std::vector<BenchmarkRow> rows;
        rows.reserve(3 * piece_counts.size());
        std::cout << "[Scalar-reverse NUBS vs MINCO gradient benchmark]" << std::endl;
        std::cout << "direct = native internal gradient construction; full = propagated gradient in shared waypoint/time variables" << std::endl;
        for (const int piece_num : piece_counts)
        {
            rows.push_back(runOne<2>(piece_num, runs));
            rows.push_back(runOne<3>(piece_num, runs));
            rows.push_back(runOne<4>(piece_num, runs));
        }
        writeCsv(csv_path, rows);
        std::cout << "wrote CSV: " << csv_path << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "bench_gradient_propagation failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
