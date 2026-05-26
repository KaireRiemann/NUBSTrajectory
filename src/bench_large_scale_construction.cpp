#include "NUBSTrajectory.hpp"
#include "large_scale_traj_opt/traj_min_jerk.hpp"
#include "large_scale_traj_opt/traj_min_snap.hpp"
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

template <int S>
struct LargeScaleUniformAdapter;

template <>
struct LargeScaleUniformAdapter<3>
{
    using Optimizer = min_jerk::JerkOpt;
    using HeadState = Eigen::Matrix3d;

    static const char *name()
    {
        return "JerkOpt";
    }

    static void reset(Optimizer &opt,
                      const Eigen::Matrix<double, 3, 3> &head_state,
                      const Eigen::Matrix<double, 3, 3> &tail_state,
                      const int piece_num)
    {
        opt.reset(HeadState(head_state), HeadState(tail_state), piece_num);
    }

    static void generate(Optimizer &opt,
                         const Eigen::Matrix<double, Eigen::Dynamic, 3> &inner_points,
                         const Eigen::VectorXd &durations)
    {
        opt.generate(inner_points.transpose(), durations);
    }
};

template <>
struct LargeScaleUniformAdapter<4>
{
    using Optimizer = min_snap::SnapOpt;
    using HeadState = Eigen::Matrix<double, 3, 4>;

    static const char *name()
    {
        return "SnapOpt";
    }

    static void reset(Optimizer &opt,
                      const Eigen::Matrix<double, 3, 4> &head_state,
                      const Eigen::Matrix<double, 3, 4> &tail_state,
                      const int piece_num)
    {
        opt.reset(HeadState(head_state), HeadState(tail_state), piece_num);
    }

    static void generate(Optimizer &opt,
                         const Eigen::Matrix<double, Eigen::Dynamic, 3> &inner_points,
                         const Eigen::VectorXd &durations)
    {
        opt.generate(inner_points.transpose(), durations);
    }
};

template <typename T>
inline void keepAlive(const T &value)
{
    volatile const void *sink = static_cast<const void *>(&value);
    (void)sink;
}

struct ConstructionBenchmarkRow
{
    int order = 0;
    int piece_num = 0;
    int runs = 0;
    double ubs_avg_us = 0.0;
    double nubs_avg_us = 0.0;
    double minco_avg_us = 0.0;
    double large_scale_avg_us = 0.0;
    std::string large_scale_name;
};

void printRow(const ConstructionBenchmarkRow &row)
{
    std::cout << std::fixed << std::setprecision(3)
              << "S=" << row.order
              << ", M=" << row.piece_num
              << ", runs=" << row.runs
              << ", ubs_avg_us=" << row.ubs_avg_us
              << ", nubs_avg_us=" << row.nubs_avg_us
              << ", minco_avg_us=" << row.minco_avg_us
              << ", " << row.large_scale_name
              << "_avg_us=" << row.large_scale_avg_us
              << ", ubs_vs_nubs_speedup="
              << row.nubs_avg_us / row.ubs_avg_us << "x"
              << ", ubs_vs_large_speedup="
              << row.large_scale_avg_us / row.ubs_avg_us << "x"
              << ", ubs_vs_minco_speedup="
              << row.minco_avg_us / row.ubs_avg_us << "x"
              << std::endl;
}

void writeCsv(const std::string &path,
              const std::vector<ConstructionBenchmarkRow> &rows)
{
    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error("failed to open CSV output: " + path);
    }

    out << "order,piece_count,runs,"
        << "ubs_avg_us,nubs_avg_us,minco_avg_us,large_scale_avg_us,"
        << "large_scale_name,"
        << "ubs_vs_nubs_speedup,ubs_vs_large_speedup,ubs_vs_minco_speedup\n";
    out << std::fixed << std::setprecision(9);
    for (const auto &row : rows)
    {
        out << row.order << ','
            << row.piece_num << ','
            << row.runs << ','
            << row.ubs_avg_us << ','
            << row.nubs_avg_us << ','
            << row.minco_avg_us << ','
            << row.large_scale_avg_us << ','
            << row.large_scale_name << ','
            << row.nubs_avg_us / row.ubs_avg_us << ','
            << row.large_scale_avg_us / row.ubs_avg_us << ','
            << row.minco_avg_us / row.ubs_avg_us << '\n';
    }
}

template <int S>
ConstructionBenchmarkRow runOneConstructionBenchmark(const int piece_num,
                                                     const int runs)
{
    using Adapter = LargeScaleUniformAdapter<S>;
    const auto data = nubs_test::makeRandomProblem<3>(
        S, piece_num, 11000 + 100 * S + piece_num);
    const double total_duration = data.durations.sum();
    const Eigen::VectorXd uniform_durations =
        nubs::UBSTrajectoryT<3, S>::uniformDurations(piece_num,
                                                     total_duration);

    const auto ubs = nubs_test::measureRepeated(
        [&]()
        {
            nubs::UBSTrajectoryT<3, S> trajectory;
            Eigen::MatrixXd control_points;
            trajectory.generateUniform(data.inner_points, data.head_state,
                                       data.tail_state, total_duration,
                                       control_points);
            keepAlive(trajectory);
            keepAlive(control_points);
        },
        runs);

    const auto full_nubs = nubs_test::measureRepeated(
        [&]()
        {
            nubs::NUBSTrajectoryT<3, S> trajectory;
            Eigen::MatrixXd control_points;
            trajectory.generate(data.inner_points, data.head_state,
                                data.tail_state, uniform_durations,
                                control_points);
            keepAlive(trajectory);
            keepAlive(control_points);
        },
        runs);

    const auto minco = nubs_test::measureRepeated(
        [&]()
        {
            nubs_test::SplineMincoAdapter<3> trajectory;
            trajectory.generate(data.inner_points, data.head_state,
                                data.tail_state, uniform_durations, S);
            keepAlive(trajectory);
        },
        runs);

    const auto large_scale = nubs_test::measureRepeated(
        [&]()
        {
            typename Adapter::Optimizer optimizer;
            Adapter::reset(optimizer, data.head_state, data.tail_state,
                           piece_num);
            Adapter::generate(optimizer, data.inner_points,
                              uniform_durations);
            keepAlive(optimizer);
        },
        runs);

    ConstructionBenchmarkRow row;
    row.order = S;
    row.piece_num = piece_num;
    row.runs = runs;
    row.ubs_avg_us = ubs.avg_us;
    row.nubs_avg_us = full_nubs.avg_us;
    row.minco_avg_us = minco.avg_us;
    row.large_scale_avg_us = large_scale.avg_us;
    row.large_scale_name = Adapter::name();
    printRow(row);
    return row;
}

std::vector<ConstructionBenchmarkRow> runConstructionBenchmark(const int runs)
{
    const std::vector<int> medium_counts = {
        2, 4, 8, 16, 32, 64, 128, 256, 512};
    const std::vector<int> large_counts = {
        1000, 2000, 3000, 4000, 5000};

    std::cout << "[large-scale construction benchmark]" << std::endl;
    std::cout << "piece counts: 2..512 powers of two, plus 1000..5000"
              << std::endl;
    std::cout << "runs per piece count: " << runs << std::endl;

    std::vector<ConstructionBenchmarkRow> rows;
    rows.reserve(2 * (medium_counts.size() + large_counts.size()));
    for (const int piece_num : medium_counts)
    {
        rows.push_back(runOneConstructionBenchmark<3>(piece_num, runs));
        rows.push_back(runOneConstructionBenchmark<4>(piece_num, runs));
    }
    for (const int piece_num : large_counts)
    {
        rows.push_back(runOneConstructionBenchmark<3>(piece_num, runs));
        rows.push_back(runOneConstructionBenchmark<4>(piece_num, runs));
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
        const auto rows = runConstructionBenchmark(runs);
        if (!csv_path.empty())
        {
            writeCsv(csv_path, rows);
            std::cout << "wrote CSV: " << csv_path << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "bench_large_scale_construction failed: "
                  << e.what() << std::endl;
        return 1;
    }
    return 0;
}
