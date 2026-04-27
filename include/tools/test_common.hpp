#ifndef NUBS_TEST_COMMON_HPP
#define NUBS_TEST_COMMON_HPP

#include "NUBSTrajectory.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace nubs_test
{

struct TimingStats
{
    double total_us = 0.0;
    double avg_us = 0.0;
    double min_us = 0.0;
    double max_us = 0.0;
    int runs = 0;
};

struct ErrorStats
{
    double max_norm = 0.0;
    double avg_norm = 0.0;
    int samples = 0;
};

inline void require(const bool condition, const std::string &message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

inline void requireNear(const double actual,
                        const double expected,
                        const double abs_tol,
                        const double rel_tol,
                        const std::string &message)
{
    const double diff = std::abs(actual - expected);
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    if (diff > abs_tol + rel_tol * scale)
    {
        std::ostringstream oss;
        oss << message << ": actual=" << actual << ", expected=" << expected
            << ", diff=" << diff << ", tolerance=" << abs_tol + rel_tol * scale;
        throw std::runtime_error(oss.str());
    }
}

template <typename DerivedA, typename DerivedB>
inline void requireMatrixNear(const Eigen::MatrixBase<DerivedA> &actual,
                              const Eigen::MatrixBase<DerivedB> &expected,
                              const double abs_tol,
                              const double rel_tol,
                              const std::string &message)
{
    require(actual.rows() == expected.rows() && actual.cols() == expected.cols(),
            message + ": matrix shape mismatch");

    for (Eigen::Index r = 0; r < actual.rows(); ++r)
    {
        for (Eigen::Index c = 0; c < actual.cols(); ++c)
        {
            std::ostringstream entry_msg;
            entry_msg << message << "(" << r << "," << c << ")";
            requireNear(actual(r, c), expected(r, c), abs_tol, rel_tol, entry_msg.str());
        }
    }
}

inline const char *orderName(const int sys_order)
{
    switch (sys_order)
    {
    case 2:
        return "Cubic";
    case 3:
        return "Quintic";
    case 4:
        return "Septic";
    default:
        return "NUBS";
    }
}

inline TimingStats measureRepeated(const std::function<void()> &func,
                                   const int runs)
{
    TimingStats result;
    result.runs = runs;
    if (runs <= 0)
    {
        return result;
    }

    std::vector<double> times;
    times.reserve(runs);
    const auto total_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < runs; ++i)
    {
        const auto start = std::chrono::high_resolution_clock::now();
        func();
        const auto end = std::chrono::high_resolution_clock::now();
        times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() /
                        1000.0);
    }
    const auto total_end = std::chrono::high_resolution_clock::now();

    result.total_us =
        std::chrono::duration_cast<std::chrono::nanoseconds>(total_end - total_start).count() /
        1000.0;
    result.avg_us = std::accumulate(times.begin(), times.end(), 0.0) /
                    static_cast<double>(runs);
    result.min_us = *std::min_element(times.begin(), times.end());
    result.max_us = *std::max_element(times.begin(), times.end());
    return result;
}

inline TimingStats measureQueries(const std::function<void(int)> &func,
                                  const int queries)
{
    TimingStats result;
    result.runs = queries;
    if (queries <= 0)
    {
        return result;
    }

    std::vector<double> times;
    times.reserve(queries);
    const auto total_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < queries; ++i)
    {
        const auto start = std::chrono::high_resolution_clock::now();
        func(i);
        const auto end = std::chrono::high_resolution_clock::now();
        times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() /
                        1000.0);
    }
    const auto total_end = std::chrono::high_resolution_clock::now();

    result.total_us =
        std::chrono::duration_cast<std::chrono::nanoseconds>(total_end - total_start).count() /
        1000.0;
    result.avg_us = std::accumulate(times.begin(), times.end(), 0.0) /
                    static_cast<double>(queries);
    result.min_us = *std::min_element(times.begin(), times.end());
    result.max_us = *std::max_element(times.begin(), times.end());
    return result;
}

inline void printTiming(const std::string &label, const TimingStats &timing)
{
    std::cout << std::fixed << std::setprecision(3)
              << label << ": runs=" << timing.runs
              << ", total=" << timing.total_us << " us"
              << ", avg=" << timing.avg_us << " us"
              << ", min=" << timing.min_us << " us"
              << ", max=" << timing.max_us << " us" << std::endl;
}

inline double relativeError(const double actual, const double expected)
{
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    return std::abs(actual - expected) / scale;
}

template <typename DerivedA, typename DerivedB>
inline double maxAbsCoeffDiff(const Eigen::MatrixBase<DerivedA> &actual,
                              const Eigen::MatrixBase<DerivedB> &expected)
{
    if (actual.rows() == 0 || actual.cols() == 0)
    {
        return 0.0;
    }
    return (actual - expected).cwiseAbs().maxCoeff();
}

template <typename Derived>
inline void accumulateSampleError(ErrorStats &stats,
                                  const Eigen::MatrixBase<Derived> &diff)
{
    const double norm = diff.norm();
    stats.max_norm = std::max(stats.max_norm, norm);
    stats.avg_norm += norm;
    ++stats.samples;
}

inline void finalizeErrorStats(ErrorStats &stats)
{
    if (stats.samples > 0)
    {
        stats.avg_norm /= static_cast<double>(stats.samples);
    }
}

inline void printErrorStats(const std::string &label, const ErrorStats &stats)
{
    std::cout << std::scientific << std::setprecision(4)
              << label << ": samples=" << stats.samples
              << ", max_norm=" << stats.max_norm
              << ", avg_norm=" << stats.avg_norm << std::endl;
}

inline Eigen::VectorXd makeDurations(const int piece_num)
{
    Eigen::VectorXd T(piece_num);
    for (int i = 0; i < piece_num; ++i)
    {
        T(i) = 0.72 + 0.11 * static_cast<double>((i * 3 + 1) % 5) +
               0.035 * static_cast<double>(i);
    }
    return T;
}

inline double cumulativeTime(const Eigen::VectorXd &T, const int prefix_count)
{
    double t = 0.0;
    for (int i = 0; i < prefix_count; ++i)
    {
        t += T(i);
    }
    return t;
}

template <int Dim>
struct ProblemData
{
    int sys_order = 0;
    Eigen::VectorXd durations;
    Eigen::Matrix<double, Eigen::Dynamic, Dim> waypoints;
    Eigen::Matrix<double, Eigen::Dynamic, Dim> inner_points;
    Eigen::Matrix<double, Dim, Eigen::Dynamic> head_state;
    Eigen::Matrix<double, Dim, Eigen::Dynamic> tail_state;
};

template <int Dim>
inline ProblemData<Dim> makeProblem(const int sys_order, const int piece_num)
{
    ProblemData<Dim> data;
    data.sys_order = sys_order;
    data.durations = makeDurations(piece_num);

    data.waypoints.resize(piece_num + 1, Dim);
    for (int i = 0; i <= piece_num; ++i)
    {
        for (int d = 0; d < Dim; ++d)
        {
            const double x = static_cast<double>(i + 1);
            const double y = static_cast<double>(d + 1);
            data.waypoints(i, d) =
                0.42 * static_cast<double>(i) +
                0.31 * std::sin(0.37 * x * y) +
                0.19 * std::cos(0.23 * (x + 2.0 * y));
        }
    }

    data.inner_points.resize(std::max(0, piece_num - 1), Dim);
    for (int i = 1; i < piece_num; ++i)
    {
        data.inner_points.row(i - 1) = data.waypoints.row(i);
    }

    data.head_state.resize(Dim, sys_order);
    data.tail_state.resize(Dim, sys_order);
    data.head_state.col(0) = data.waypoints.row(0).transpose();
    data.tail_state.col(0) = data.waypoints.row(piece_num).transpose();
    for (int der = 1; der < sys_order; ++der)
    {
        for (int d = 0; d < Dim; ++d)
        {
            const double scale = 0.10 / static_cast<double>(der);
            data.head_state(d, der) =
                scale * std::sin(0.29 * static_cast<double>((d + 1) * (der + 2)));
            data.tail_state(d, der) =
                scale * std::cos(0.31 * static_cast<double>((d + 2) * (der + 1)));
        }
    }
    return data;
}

template <int Dim>
inline double energyFor(const int sys_order,
                        const Eigen::Matrix<double, Eigen::Dynamic, Dim> &inner_points,
                        const Eigen::Matrix<double, Dim, Eigen::Dynamic> &head_state,
                        const Eigen::Matrix<double, Dim, Eigen::Dynamic> &tail_state,
                        const Eigen::VectorXd &durations)
{
    nubs::NUBSTrajectory<Dim> trajectory(sys_order);
    Eigen::MatrixXd control_points;
    trajectory.generate(inner_points, head_state, tail_state, durations, control_points);
    return trajectory.getEnergy();
}

template <int Dim, int S>
inline double energyForFixed(
    const Eigen::Matrix<double, Eigen::Dynamic, Dim> &inner_points,
    const Eigen::Matrix<double, Dim, Eigen::Dynamic> &head_state,
    const Eigen::Matrix<double, Dim, Eigen::Dynamic> &tail_state,
    const Eigen::VectorXd &durations)
{
    nubs::NUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd control_points;
    trajectory.generate(inner_points, head_state, tail_state, durations, control_points);
    return trajectory.getEnergy();
}

inline double previousTime(const double t)
{
    return std::nextafter(t, -std::numeric_limits<double>::infinity());
}

} // namespace nubs_test

#endif
