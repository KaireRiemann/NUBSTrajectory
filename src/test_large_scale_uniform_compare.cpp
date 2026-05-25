#include "NUBSTrajectory.hpp"
#include "large_scale_traj_opt/traj_min_jerk.hpp"
#include "large_scale_traj_opt/traj_min_snap.hpp"
#include "tools/minco_adapter.hpp"
#include "tools/optimization_test_cases.hpp"

#include <Eigen/Dense>
#include <iostream>
#include <stdexcept>

namespace
{

template <int S>
struct LargeScaleUniformAdapter;

template <>
struct LargeScaleUniformAdapter<3>
{
    using Optimizer = min_jerk::JerkOpt;
    using Trajectory = min_jerk::Trajectory;
    using HeadState = Eigen::Matrix3d;

    static const char *name()
    {
        return "large_scale_traj_opt JerkOpt";
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
    using Trajectory = min_snap::Trajectory;
    using HeadState = Eigen::Matrix<double, 3, 4>;

    static const char *name()
    {
        return "large_scale_traj_opt SnapOpt";
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

template <int S>
Eigen::Vector3d evaluateLargeScale(
    typename LargeScaleUniformAdapter<S>::Trajectory &trajectory,
    const double t,
    const int derivative)
{
    switch (derivative)
    {
    case 0:
        return trajectory.getPos(t);
    case 1:
        return trajectory.getVel(t);
    case 2:
        return trajectory.getAcc(t);
    default:
        throw std::runtime_error("large-scale reference only exposes derivatives up to acceleration.");
    }
}

template <int S>
typename LargeScaleUniformAdapter<S>::Trajectory buildLargeScaleTrajectory(
    const nubs_test::ProblemData<3> &data,
    const Eigen::VectorXd &durations)
{
    using Adapter = LargeScaleUniformAdapter<S>;
    typename Adapter::Optimizer optimizer;
    Adapter::reset(optimizer, data.head_state, data.tail_state,
                   static_cast<int>(durations.size()));
    Adapter::generate(optimizer, data.inner_points, durations);

    typename Adapter::Trajectory trajectory;
    optimizer.getTraj(trajectory);
    return trajectory;
}

template <int S>
void compareUniformConsistency(const int piece_num,
                               const unsigned int seed)
{
    using Adapter = LargeScaleUniformAdapter<S>;
    const auto data = nubs_test::makeRandomProblem<3>(S, piece_num, seed);
    const double total_duration = data.durations.sum();
    const Eigen::VectorXd uniform_durations =
        nubs::UniformNUBSTrajectoryT<3, S>::uniformDurations(piece_num,
                                                             total_duration);

    nubs::UniformNUBSTrajectoryT<3, S> nubs_traj;
    Eigen::MatrixXd control_points;
    nubs_traj.generateUniform(data.inner_points, data.head_state,
                              data.tail_state, total_duration,
                              control_points);
    nubs_traj.prepareEvaluationCache();
    nubs_traj.prepareEvaluationCache();

    nubs::NUBSTrajectoryT<3, S> full_nubs_traj;
    Eigen::MatrixXd full_control_points;
    full_nubs_traj.generate(data.inner_points, data.head_state,
                            data.tail_state, uniform_durations,
                            full_control_points);

    nubs_test::SplineMincoAdapter<3> minco_traj;
    minco_traj.generate(data.inner_points, data.head_state,
                        data.tail_state, uniform_durations, S);

    typename Adapter::Optimizer reference_optimizer;
    Adapter::reset(reference_optimizer, data.head_state, data.tail_state,
                   piece_num);
    Adapter::generate(reference_optimizer, data.inner_points,
                      uniform_durations);
    typename Adapter::Trajectory reference_traj;
    reference_optimizer.getTraj(reference_traj);

    nubs_test::requireNear(nubs_traj.getEnergy(),
                           reference_optimizer.getObjective(),
                           2.0e-7, 2.0e-9,
                           "NUBS uniform vs large-scale objective");
    nubs_test::requireNear(full_nubs_traj.getEnergy(),
                           reference_optimizer.getObjective(),
                           2.0e-7, 2.0e-9,
                           "full NUBS vs large-scale objective");
    nubs_test::requireNear(minco_traj.getEnergy(),
                           reference_optimizer.getObjective(),
                           2.0e-7, 2.0e-9,
                           "MINCO vs large-scale objective");

    nubs_test::ErrorStats pos_error;
    nubs_test::ErrorStats vel_error;
    nubs_test::ErrorStats acc_error;
    nubs_test::ErrorStats full_pos_error;
    nubs_test::ErrorStats minco_pos_error;
    constexpr int samples = 80;
    for (int i = 0; i <= samples; ++i)
    {
        double t = total_duration * static_cast<double>(i) /
                   static_cast<double>(samples);
        if (i == samples)
        {
            t = nubs_test::previousTime(total_duration);
        }

        const Eigen::Vector3d nubs_pos = nubs_traj.evaluate(t, 0);
        const Eigen::Vector3d nubs_vel = nubs_traj.evaluate(t, 1);
        const Eigen::Vector3d nubs_acc = nubs_traj.evaluate(t, 2);
        const Eigen::Vector3d full_pos = full_nubs_traj.evaluate(t, 0);
        const Eigen::Vector3d minco_pos = minco_traj.evaluate(t, 0);
        const Eigen::Vector3d ref_pos = evaluateLargeScale<S>(reference_traj, t, 0);
        nubs_test::accumulateSampleError(
            pos_error, nubs_pos - ref_pos);
        nubs_test::accumulateSampleError(
            vel_error, nubs_vel - evaluateLargeScale<S>(reference_traj, t, 1));
        nubs_test::accumulateSampleError(
            acc_error, nubs_acc - evaluateLargeScale<S>(reference_traj, t, 2));
        nubs_test::accumulateSampleError(full_pos_error, full_pos - ref_pos);
        nubs_test::accumulateSampleError(minco_pos_error, minco_pos - ref_pos);
    }
    nubs_test::finalizeErrorStats(pos_error);
    nubs_test::finalizeErrorStats(vel_error);
    nubs_test::finalizeErrorStats(acc_error);
    nubs_test::finalizeErrorStats(full_pos_error);
    nubs_test::finalizeErrorStats(minco_pos_error);

    nubs_test::require(pos_error.max_norm < 2.0e-7,
                       "NUBS uniform vs large-scale position mismatch");
    nubs_test::require(vel_error.max_norm < 2.0e-6,
                       "NUBS uniform vs large-scale velocity mismatch");
    nubs_test::require(acc_error.max_norm < 2.0e-5,
                       "NUBS uniform vs large-scale acceleration mismatch");
    nubs_test::require(full_pos_error.max_norm < 2.0e-7,
                       "full NUBS vs large-scale position mismatch");
    nubs_test::require(minco_pos_error.max_norm < 2.0e-7,
                       "MINCO vs large-scale position mismatch");

    std::cout << "  S=" << S << ", M=" << piece_num
              << ", ref=" << Adapter::name() << std::endl;
    nubs_test::printErrorStats("    pos error", pos_error);
    nubs_test::printErrorStats("    vel error", vel_error);
    nubs_test::printErrorStats("    acc error", acc_error);
    nubs_test::printErrorStats("    full NUBS pos error", full_pos_error);
    nubs_test::printErrorStats("    MINCO pos error", minco_pos_error);
}

template <int S>
void benchmarkUniformConstructionAndEvaluation(const int piece_num,
                                               const int construction_runs,
                                               const int query_runs)
{
    using Adapter = LargeScaleUniformAdapter<S>;
    const auto data = nubs_test::makeRandomProblem<3>(
        S, piece_num, 9000 + 100 * S + piece_num);
    const double total_duration = data.durations.sum();
    const Eigen::VectorXd uniform_durations =
        nubs::UniformNUBSTrajectoryT<3, S>::uniformDurations(piece_num,
                                                             total_duration);

    const auto nubs_construct = nubs_test::measureRepeated(
        [&]()
        {
            nubs::UniformNUBSTrajectoryT<3, S> trajectory;
            Eigen::MatrixXd control_points;
            trajectory.generateUniform(data.inner_points, data.head_state,
                                       data.tail_state, total_duration,
                                       control_points);
            volatile double sink = control_points.squaredNorm();
            (void)sink;
        },
        construction_runs);

    const auto full_nubs_construct = nubs_test::measureRepeated(
        [&]()
        {
            nubs::NUBSTrajectoryT<3, S> trajectory;
            Eigen::MatrixXd control_points;
            trajectory.generate(data.inner_points, data.head_state,
                                data.tail_state, uniform_durations,
                                control_points);
            volatile double sink = control_points.squaredNorm();
            (void)sink;
        },
        construction_runs);

    const auto minco_construct = nubs_test::measureRepeated(
        [&]()
        {
            nubs_test::SplineMincoAdapter<3> trajectory;
            trajectory.generate(data.inner_points, data.head_state,
                                data.tail_state, uniform_durations, S);
            volatile double sink = trajectory.getEnergy();
            (void)sink;
        },
        construction_runs);

    const auto large_construct = nubs_test::measureRepeated(
        [&]()
        {
            typename Adapter::Optimizer optimizer;
            Adapter::reset(optimizer, data.head_state, data.tail_state,
                           piece_num);
            Adapter::generate(optimizer, data.inner_points,
                              uniform_durations);
            volatile double sink = optimizer.getObjective();
            (void)sink;
        },
        construction_runs);

    nubs::UniformNUBSTrajectoryT<3, S> nubs_traj;
    Eigen::MatrixXd control_points;
    nubs_traj.generateUniform(data.inner_points, data.head_state,
                              data.tail_state, total_duration,
                              control_points);
    nubs::NUBSTrajectoryT<3, S> full_nubs_traj;
    Eigen::MatrixXd full_control_points;
    full_nubs_traj.generate(data.inner_points, data.head_state,
                            data.tail_state, uniform_durations,
                            full_control_points);
    nubs_test::SplineMincoAdapter<3> minco_traj;
    minco_traj.generate(data.inner_points, data.head_state,
                        data.tail_state, uniform_durations, S);
    typename Adapter::Trajectory large_traj =
        buildLargeScaleTrajectory<S>(data, uniform_durations);

    const auto sample_time = [&](const int i)
    {
        const double alpha =
            static_cast<double>((i * 37) % query_runs) /
            static_cast<double>(query_runs);
        return std::min(nubs_test::previousTime(total_duration),
                        total_duration * alpha);
    };

    const auto nubs_eval_pos = nubs_test::measureQueries(
        [&](const int i)
        {
            const double t = sample_time(i);
            const Eigen::Vector3d value = nubs_traj.evaluate(t, 0);
            volatile double sink = value.squaredNorm();
            (void)sink;
        },
        query_runs);

    const auto large_eval_pos = nubs_test::measureQueries(
        [&](const int i)
        {
            const double t = sample_time(i);
            const Eigen::Vector3d value = large_traj.getPos(t);
            volatile double sink = value.squaredNorm();
            (void)sink;
        },
        query_runs);

    const auto full_nubs_eval_pos = nubs_test::measureQueries(
        [&](const int i)
        {
            const double t = sample_time(i);
            const Eigen::Vector3d value = full_nubs_traj.evaluate(t, 0);
            volatile double sink = value.squaredNorm();
            (void)sink;
        },
        query_runs);

    const auto minco_eval_pos = nubs_test::measureQueries(
        [&](const int i)
        {
            const double t = sample_time(i);
            const Eigen::Vector3d value = minco_traj.evaluate(t, 0);
            volatile double sink = value.squaredNorm();
            (void)sink;
        },
        query_runs);

    const auto nubs_eval_pva = nubs_test::measureQueries(
        [&](const int i)
        {
            const double t = sample_time(i);
            Eigen::Vector3d pos;
            Eigen::Vector3d vel;
            Eigen::Vector3d acc;
            nubs_traj.evaluatePVA(t, pos, vel, acc);
            volatile double sink =
                pos.squaredNorm() + vel.squaredNorm() + acc.squaredNorm();
            (void)sink;
        },
        query_runs);

    const auto full_nubs_eval_pva = nubs_test::measureQueries(
        [&](const int i)
        {
            const double t = sample_time(i);
            const Eigen::Vector3d pos = full_nubs_traj.evaluate(t, 0);
            const Eigen::Vector3d vel = full_nubs_traj.evaluate(t, 1);
            const Eigen::Vector3d acc = full_nubs_traj.evaluate(t, 2);
            volatile double sink =
                pos.squaredNorm() + vel.squaredNorm() + acc.squaredNorm();
            (void)sink;
        },
        query_runs);

    const auto minco_eval_pva = nubs_test::measureQueries(
        [&](const int i)
        {
            const double t = sample_time(i);
            const Eigen::Vector3d pos = minco_traj.evaluate(t, 0);
            const Eigen::Vector3d vel = minco_traj.evaluate(t, 1);
            const Eigen::Vector3d acc = minco_traj.evaluate(t, 2);
            volatile double sink =
                pos.squaredNorm() + vel.squaredNorm() + acc.squaredNorm();
            (void)sink;
        },
        query_runs);

    const auto large_eval_pva = nubs_test::measureQueries(
        [&](const int i)
        {
            const double t = sample_time(i);
            const Eigen::Vector3d pos = large_traj.getPos(t);
            const Eigen::Vector3d vel = large_traj.getVel(t);
            const Eigen::Vector3d acc = large_traj.getAcc(t);
            volatile double sink =
                pos.squaredNorm() + vel.squaredNorm() + acc.squaredNorm();
            (void)sink;
        },
        query_runs);

    std::cout << "  [speed] S=" << S << ", M=" << piece_num
              << ", ref=" << Adapter::name() << std::endl;
    nubs_test::printTiming("    uniform NUBS construct", nubs_construct);
    nubs_test::printTiming("    full NUBS construct", full_nubs_construct);
    nubs_test::printTiming("    MINCO construct", minco_construct);
    nubs_test::printTiming("    large construct", large_construct);
    nubs_test::printTiming("    uniform NUBS eval pos", nubs_eval_pos);
    nubs_test::printTiming("    full NUBS eval pos", full_nubs_eval_pos);
    nubs_test::printTiming("    MINCO eval pos", minco_eval_pos);
    nubs_test::printTiming("    large eval pos", large_eval_pos);
    nubs_test::printTiming("    uniform NUBS eval pva", nubs_eval_pva);
    nubs_test::printTiming("    full NUBS eval pva", full_nubs_eval_pva);
    nubs_test::printTiming("    MINCO eval pva", minco_eval_pva);
    nubs_test::printTiming("    large eval pva", large_eval_pva);
}

void runLargeScaleUniformCompare()
{
    std::cout << "[large-scale optimizer vs uniform NUBS]" << std::endl;
    compareUniformConsistency<3>(6, 7303);
    compareUniformConsistency<4>(6, 7404);

    constexpr int construction_runs = 500;
    constexpr int query_runs = 8000;
    for (const int piece_num : {4, 8, 16, 32})
    {
        benchmarkUniformConstructionAndEvaluation<3>(
            piece_num, construction_runs, query_runs);
        benchmarkUniformConstructionAndEvaluation<4>(
            piece_num, construction_runs, query_runs);
    }
    std::cout << "  PASS" << std::endl;
}

} // namespace

int main()
{
    try
    {
        runLargeScaleUniformCompare();
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_large_scale_uniform_compare failed: "
                  << e.what() << std::endl;
        return 1;
    }
    return 0;
}
