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

template <typename T>
inline void keepAlive(const T &value)
{
    volatile const void *sink = static_cast<const void *>(&value);
    (void)sink;
}

struct BenchmarkRow
{
    int order = 0;
    int piece_num = 0;
    int runs = 0;
    double nubs_avg_us = 0.0;
    double minco_avg_us = 0.0;
    double minco_vs_nubs_speedup = 0.0;
    double max_position_error = 0.0;
    double energy_abs_error = 0.0;
    double energy_rel_error = 0.0;
};

void printRow(const BenchmarkRow &row)
{
    std::cout << std::fixed << std::setprecision(3)
              << "S=" << row.order
              << ", M=" << row.piece_num
              << ", runs=" << row.runs
              << ", nubs_avg_us=" << row.nubs_avg_us
              << ", minco_avg_us=" << row.minco_avg_us
              << ", minco_vs_nubs_speedup="
              << row.minco_vs_nubs_speedup << "x"
              << std::scientific << std::setprecision(3)
              << ", max_pos_err=" << row.max_position_error
              << ", energy_abs_err=" << row.energy_abs_error
              << std::endl;
}

void writeCsv(const std::string &path,
              const std::vector<BenchmarkRow> &rows)
{
    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error("failed to open CSV output: " + path);
    }

    out << "order,piece_count,runs,"
        << "nubs_avg_us,minco_avg_us,minco_vs_nubs_speedup,"
        << "max_position_error,energy_abs_error,energy_rel_error\n";
    out << std::fixed << std::setprecision(9);
    for (const auto &row : rows)
    {
        out << row.order << ','
            << row.piece_num << ','
            << row.runs << ','
            << row.nubs_avg_us << ','
            << row.minco_avg_us << ','
            << row.minco_vs_nubs_speedup << ','
            << std::scientific << std::setprecision(17)
            << row.max_position_error << ','
            << row.energy_abs_error << ','
            << row.energy_rel_error
            << std::fixed << std::setprecision(9) << '\n';
    }
}

template <int S>
void checkEquivalence(const nubs_test::ProblemData<3> &data,
                      double &max_position_error,
                      double &energy_abs_error,
                      double &energy_rel_error)
{
    nubs::NUBSTrajectoryT<3, S> nubs_traj;
    Eigen::MatrixXd control_points;
    nubs_traj.generate(data.inner_points, data.head_state, data.tail_state,
                       data.durations, control_points);

    nubs_test::SplineMincoAdapter<3> minco_traj;
    minco_traj.generate(data.inner_points, data.head_state, data.tail_state,
                        data.durations, S);

    constexpr int samples = 96;
    const double total_duration = data.durations.sum();
    max_position_error = 0.0;
    for (int i = 0; i <= samples; ++i)
    {
        double t = total_duration * static_cast<double>(i) /
                   static_cast<double>(samples);
        if (i == samples)
        {
            t = nubs_test::previousTime(total_duration);
        }
        const double err =
            (nubs_traj.evaluate(t, 0) - minco_traj.evaluate(t, 0)).norm();
        max_position_error = std::max(max_position_error, err);
    }

    const double nubs_energy = nubs_traj.getEnergy();
    const double minco_energy = minco_traj.getEnergy();
    energy_abs_error = std::abs(nubs_energy - minco_energy);
    energy_rel_error = nubs_test::relativeError(nubs_energy, minco_energy);

    nubs_test::require(max_position_error < 1.0e-7,
                       "NUBS vs MINCO non-uniform position mismatch");
    nubs_test::require(energy_rel_error < 1.0e-9,
                       "NUBS vs MINCO non-uniform energy mismatch");
}

template <int S>
BenchmarkRow runOneBenchmark(const int piece_num,
                             const int runs)
{
    const auto data = nubs_test::makeRandomProblem<3>(
        S, piece_num, 13000 + 100 * S + piece_num);

    const auto nubs_time = nubs_test::measureRepeated(
        [&]()
        {
            nubs::NUBSTrajectoryT<3, S> trajectory;
            Eigen::MatrixXd control_points;
            trajectory.generate(data.inner_points, data.head_state,
                                data.tail_state, data.durations,
                                control_points);
            keepAlive(trajectory);
            keepAlive(control_points);
        },
        runs);

    const auto minco_time = nubs_test::measureRepeated(
        [&]()
        {
            nubs_test::SplineMincoAdapter<3> trajectory;
            trajectory.generate(data.inner_points, data.head_state,
                                data.tail_state, data.durations, S);
            keepAlive(trajectory);
        },
        runs);

    BenchmarkRow row;
    row.order = S;
    row.piece_num = piece_num;
    row.runs = runs;
    row.nubs_avg_us = nubs_time.avg_us;
    row.minco_avg_us = minco_time.avg_us;
    row.minco_vs_nubs_speedup = row.minco_avg_us / row.nubs_avg_us;
    checkEquivalence<S>(data, row.max_position_error,
                        row.energy_abs_error, row.energy_rel_error);
    printRow(row);
    return row;
}

std::vector<BenchmarkRow> runBenchmark(const int runs)
{
    const std::vector<int> piece_counts = {2, 4, 8, 16, 32, 64};
    std::cout << "[non-uniform construction benchmark]" << std::endl;
    std::cout << "methods: NUBS, MINCO" << std::endl;
    std::cout << "piece counts: 2, 4, 8, 16, 32, 64" << std::endl;
    std::cout << "runs per piece count: " << runs << std::endl;

    std::vector<BenchmarkRow> rows;
    rows.reserve(2 * piece_counts.size());
    for (const int piece_num : piece_counts)
    {
        rows.push_back(runOneBenchmark<3>(piece_num, runs));
        rows.push_back(runOneBenchmark<4>(piece_num, runs));
    }
    return rows;
}

} // namespace

int main(int argc, char **argv)
{
    int runs = 1000;
    std::string csv_path;
    if (argc > 1)
    {
        runs = std::max(1, std::atoi(argv[1]));
    }
    if (argc > 2)
    {
        csv_path = argv[2];
    }

    try
    {
        const auto rows = runBenchmark(runs);
        if (!csv_path.empty())
        {
            writeCsv(csv_path, rows);
            std::cout << "wrote CSV: " << csv_path << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "bench_nonuniform_construction failed: "
                  << e.what() << std::endl;
        return 1;
    }
    return 0;
}
