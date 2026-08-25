#ifndef NUBS_OPTIMIZATION_TEST_CASES_HPP
#define NUBS_OPTIMIZATION_TEST_CASES_HPP

#include "tools/test_common.hpp"

#include <random>

namespace nubs_test
{

template <typename Derived>
inline void requireFiniteMatrix(const Eigen::MatrixBase<Derived> &mat,
                                const std::string &message)
{
    for (Eigen::Index r = 0; r < mat.rows(); ++r)
    {
        for (Eigen::Index c = 0; c < mat.cols(); ++c)
        {
            if (!std::isfinite(mat(r, c)))
            {
                throw std::runtime_error(message);
            }
        }
    }
}

template <int Dim>
inline ProblemData<Dim> makeRandomProblem(const int sys_order,
                                          const int piece_num,
                                          const unsigned int seed,
                                          const bool small_durations = false)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> point_dist(-1.5, 1.5);
    std::uniform_real_distribution<double> state_dist(-0.2, 0.2);
    std::uniform_real_distribution<double> duration_dist(0.35, 1.45);

    ProblemData<Dim> data;
    data.sys_order = sys_order;
    data.durations.resize(piece_num);
    for (int i = 0; i < piece_num; ++i)
    {
        if (small_durations)
        {
            data.durations(i) =
                (i % 2 == 0 ? 1.0e-4 : 1.0e-3) *
                (1.0 + 0.15 * static_cast<double>(i));
        }
        else
        {
            data.durations(i) = duration_dist(gen);
        }
    }

    data.waypoints.resize(piece_num + 1, Dim);
    for (int r = 0; r <= piece_num; ++r)
    {
        for (int c = 0; c < Dim; ++c)
        {
            data.waypoints(r, c) =
                0.35 * static_cast<double>(r) + point_dist(gen);
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
            data.head_state(d, der) =
                state_dist(gen) / static_cast<double>(der);
            data.tail_state(d, der) =
                state_dist(gen) / static_cast<double>(der);
        }
    }
    return data;
}

inline double externalDurationStep(const double Ti)
{
    double h = 1.0e-6 * std::max(1.0, std::abs(Ti));
    h = std::min(h, 0.25 * Ti);
    return h > 0.0 ? h : 1.0e-8;
}

template <int Dim, int S>
inline void checkBoundaryOne(const unsigned int seed)
{
    constexpr int M = 5;
    const auto data = makeRandomProblem<Dim>(S, M, seed);

    nubs::NUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd control_points;
    trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                        data.durations, control_points);

    double max_boundary = 0.0;
    double max_waypoint = 0.0;
    for (int der = 0; der < S; ++der)
    {
        const auto head_value = trajectory.evaluate(0.0, der);
        const auto tail_value =
            trajectory.evaluate(previousTime(data.durations.sum()), der);
        max_boundary =
            std::max(max_boundary,
                     (head_value - data.head_state.col(der)).norm());
        max_boundary =
            std::max(max_boundary,
                     (tail_value - data.tail_state.col(der)).norm());
        const double tol = der <= 1 ? 1.0e-7 : 5.0e-6;
        requireMatrixNear(head_value, data.head_state.col(der),
                          tol, tol, "random head boundary");
        requireMatrixNear(tail_value, data.tail_state.col(der),
                          5.0e-6, 5.0e-6, "random tail boundary");
    }

    for (int i = 1; i < M; ++i)
    {
        const auto value = trajectory.evaluate(cumulativeTime(data.durations, i), 0);
        const Eigen::Matrix<double, Dim, 1> expected =
            data.inner_points.row(i - 1).transpose();
        max_waypoint = std::max(max_waypoint, (value - expected).norm());
        requireMatrixNear(value, expected, 1.0e-7, 1.0e-7,
                          "random waypoint");
    }

    std::cout << "  S=" << S << ", Dim=" << Dim
              << ": boundary=" << std::scientific << std::setprecision(3)
              << max_boundary << ", waypoint=" << max_waypoint << std::endl;
}

inline void runBoundaryGridCase()
{
    std::cout << "[boundary grid]" << std::endl;
    checkBoundaryOne<1, 2>(1201);
    checkBoundaryOne<2, 2>(1202);
    checkBoundaryOne<3, 2>(1203);
    checkBoundaryOne<1, 3>(1301);
    checkBoundaryOne<2, 3>(1302);
    checkBoundaryOne<3, 3>(1303);
    checkBoundaryOne<1, 4>(1401);
    checkBoundaryOne<2, 4>(1402);
    checkBoundaryOne<3, 4>(1403);
    std::cout << "  PASS" << std::endl;
}

template <int Dim, int S>
inline void checkGenericSpecializedOne(const unsigned int seed)
{
    constexpr int M = 5;
    const auto data = makeRandomProblem<Dim>(S, M, seed);

    nubs::NUBSTrajectory<Dim, 2 * S - 1> generic(S);
    nubs::NUBSTrajectoryT<Dim, S> fixed;
    Eigen::MatrixXd generic_control;
    Eigen::MatrixXd fixed_control;
    generic.generate(data.inner_points, data.head_state, data.tail_state,
                     data.durations, generic_control);
    fixed.generate(data.inner_points, data.head_state, data.tail_state,
                   data.durations, fixed_control);

    requireMatrixNear(fixed_control, generic_control,
                      2.0e-10, 2.0e-10, "generic vs fixed control");
    requireNear(fixed.getEnergy(), generic.getEnergy(),
                2.0e-9, 2.0e-9, "generic vs fixed energy");

    const double total_t = data.durations.sum();
    for (int sample = 0; sample <= 12; ++sample)
    {
        double t = total_t * static_cast<double>(sample) / 12.0;
        if (sample == 12)
        {
            t = previousTime(total_t);
        }
        for (int der = 0; der <= S; ++der)
        {
            requireMatrixNear(fixed.evaluate(t, der), generic.evaluate(t, der),
                              2.0e-8, 2.0e-8,
                              "generic vs fixed sampled derivative");
        }
    }

    double generic_cost = 0.0;
    double fixed_cost = 0.0;
    Eigen::MatrixXd generic_gp;
    Eigen::MatrixXd fixed_gp;
    Eigen::VectorXd generic_gt;
    Eigen::VectorXd fixed_gt;
    generic.getEnergyAndGrad(generic_cost, generic_gp, generic_gt);
    fixed.getEnergyAndGrad(fixed_cost, fixed_gp, fixed_gt);
    requireNear(fixed_cost, generic_cost, 2.0e-9, 2.0e-9,
                "generic vs fixed Local-AD cost");
    requireMatrixNear(fixed_gp, generic_gp, 5.0e-7, 5.0e-7,
                      "generic vs fixed point gradient");
    requireMatrixNear(fixed_gt, generic_gt, 2.0e-5, 2.0e-5,
                      "generic vs fixed time gradient");

    double generic_fixed_ratio_cost = 0.0;
    double fixed_fixed_ratio_cost = 0.0;
    Eigen::MatrixXd generic_fixed_ratio_gp;
    Eigen::MatrixXd fixed_fixed_ratio_gp;
    double generic_gtotal = 0.0;
    double fixed_gtotal = 0.0;
    generic.getEnergyAndFixedRatioGrad(generic_fixed_ratio_cost,
                                       generic_fixed_ratio_gp,
                                       generic_gtotal);
    fixed.getEnergyAndFixedRatioGrad(fixed_fixed_ratio_cost,
                                     fixed_fixed_ratio_gp,
                                     fixed_gtotal);
    requireNear(fixed_fixed_ratio_cost, generic_fixed_ratio_cost,
                2.0e-9, 2.0e-9, "generic vs fixed fixed-ratio cost");
    requireMatrixNear(fixed_fixed_ratio_gp, generic_fixed_ratio_gp,
                      5.0e-7, 5.0e-7,
                      "generic vs fixed fixed-ratio point gradient");
    requireNear(fixed_gtotal, generic_gtotal, 2.0e-5, 2.0e-5,
                "generic vs fixed total-duration gradient");
}

inline void runGenericSpecializedEquivalenceCase()
{
    std::cout << "[generic vs specialized]" << std::endl;
    checkGenericSpecializedOne<1, 2>(2201);
    checkGenericSpecializedOne<2, 2>(2202);
    checkGenericSpecializedOne<3, 2>(2203);
    checkGenericSpecializedOne<1, 3>(2301);
    checkGenericSpecializedOne<2, 3>(2302);
    checkGenericSpecializedOne<3, 3>(2303);
    checkGenericSpecializedOne<1, 4>(2401);
    checkGenericSpecializedOne<2, 4>(2402);
    checkGenericSpecializedOne<3, 4>(2403);
    std::cout << "  PASS" << std::endl;
}

template <int Dim, int S>
inline void checkExternalGradientOne(const unsigned int seed)
{
    constexpr int M = 4;
    const auto data = makeRandomProblem<Dim>(S, M, seed);

    nubs::NUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd control_points;
    trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                        data.durations, control_points);

    double cost = 0.0;
    Eigen::MatrixXd grad_points;
    Eigen::VectorXd grad_times;
    trajectory.getEnergyAndGrad(cost, grad_points, grad_times);

    Eigen::Matrix<double, Eigen::Dynamic, Dim> numeric_gp =
        Eigen::Matrix<double, Eigen::Dynamic, Dim>::Zero(M - 1, Dim);
    Eigen::VectorXd numeric_gt = Eigen::VectorXd::Zero(M);

    const double point_eps = 1.0e-6;
    for (int r = 0; r < data.inner_points.rows(); ++r)
    {
        for (int c = 0; c < Dim; ++c)
        {
            auto plus_points = data.inner_points;
            auto minus_points = data.inner_points;
            plus_points(r, c) += point_eps;
            minus_points(r, c) -= point_eps;
            const double plus_cost =
                energyForFixed<Dim, S>(plus_points, data.head_state,
                                       data.tail_state, data.durations);
            const double minus_cost =
                energyForFixed<Dim, S>(minus_points, data.head_state,
                                       data.tail_state, data.durations);
            numeric_gp(r, c) =
                (plus_cost - minus_cost) / (2.0 * point_eps);
        }
    }

    for (int i = 0; i < M; ++i)
    {
        const double h = externalDurationStep(data.durations(i));
        Eigen::VectorXd plus_times = data.durations;
        Eigen::VectorXd minus_times = data.durations;
        plus_times(i) += h;
        minus_times(i) = std::max(1.0e-8, minus_times(i) - h);
        const double plus_delta = plus_times(i) - data.durations(i);
        const double minus_delta = minus_times(i) - data.durations(i);
        const double plus_cost =
            energyForFixed<Dim, S>(data.inner_points, data.head_state,
                                   data.tail_state, plus_times);
        const double minus_cost =
            energyForFixed<Dim, S>(data.inner_points, data.head_state,
                                   data.tail_state, minus_times);
        numeric_gt(i) =
            (plus_cost - minus_cost) / (plus_delta - minus_delta);
    }

    requireMatrixNear(grad_points, numeric_gp, 5.0e-4, 2.0e-4,
                      "external point gradient");
    requireMatrixNear(grad_times, numeric_gt, 5.0e-3, 5.0e-4,
                      "external time gradient");
}

inline void runExternalGradientCase()
{
    std::cout << "[external finite-diff gradient]" << std::endl;
    checkExternalGradientOne<1, 2>(3201);
    checkExternalGradientOne<2, 2>(3202);
    checkExternalGradientOne<3, 2>(3203);
    checkExternalGradientOne<1, 3>(3301);
    checkExternalGradientOne<2, 3>(3302);
    checkExternalGradientOne<3, 3>(3303);
    checkExternalGradientOne<1, 4>(3401);
    checkExternalGradientOne<2, 4>(3402);
    checkExternalGradientOne<3, 4>(3403);
    std::cout << "  PASS" << std::endl;
}

inline void externalTotalDurationDeltas(const Eigen::VectorXd &ratios,
                                        const double total_duration,
                                        double &plus_delta,
                                        double &minus_delta)
{
    double h = 1.0e-6 * std::max(1.0, std::abs(total_duration));
    h = std::min(h, 0.25 * total_duration);
    plus_delta = h > 0.0 ? h : 1.0e-8;

    const double ratio_sum = ratios.sum();
    double min_allowed_total = 1.0e-8;
    for (int i = 0; i < ratios.size(); ++i)
    {
        const double normalized_ratio = ratios(i) / ratio_sum;
        min_allowed_total =
            std::max(min_allowed_total, 1.0e-8 / normalized_ratio);
    }

    const double allowed_minus =
        std::max(0.0, total_duration - min_allowed_total);
    const double h_minus = std::min(plus_delta, allowed_minus);
    minus_delta = h_minus > 0.0 ? -h_minus : 0.0;
}

template <int Dim, int S>
inline double fixedRatioEnergy(
    const Eigen::Matrix<double, Eigen::Dynamic, Dim> &inner_points,
    const Eigen::Matrix<double, Dim, Eigen::Dynamic> &head_state,
    const Eigen::Matrix<double, Dim, Eigen::Dynamic> &tail_state,
    const Eigen::VectorXd &ratios,
    const double total_duration)
{
    const Eigen::VectorXd durations =
        nubs::NUBSTrajectoryT<Dim, S>::durationsFromRatios(ratios,
                                                           total_duration);
    return energyForFixed<Dim, S>(inner_points, head_state, tail_state,
                                  durations);
}

template <int Dim, int S>
inline void checkFixedRatioTotalDurationOne(const unsigned int seed)
{
    constexpr int M = 5;
    const auto data = makeRandomProblem<Dim>(S, M, seed);
    const double total_duration = data.durations.sum();
    const Eigen::VectorXd ratios = 2.5 * data.durations;

    nubs::NUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd control_points;
    trajectory.generateWithTotalDuration(data.inner_points, data.head_state,
                                         data.tail_state, ratios,
                                         total_duration, control_points);

    double fixed_ratio_cost = 0.0;
    Eigen::MatrixXd fixed_ratio_gp;
    double fixed_ratio_gtotal = 0.0;
    trajectory.getEnergyAndFixedRatioGrad(fixed_ratio_cost, fixed_ratio_gp,
                                          fixed_ratio_gtotal);

    double segment_cost = 0.0;
    Eigen::MatrixXd segment_gp;
    Eigen::VectorXd segment_gt;
    trajectory.getEnergyAndGrad(segment_cost, segment_gp, segment_gt);

    const Eigen::VectorXd normalized_ratios =
        trajectory.getDurations() / trajectory.getTotalDuration();
    const double contracted_time_grad = normalized_ratios.dot(segment_gt);

    requireNear(fixed_ratio_cost, segment_cost, 2.0e-9, 2.0e-9,
                "fixed-ratio cost");
    requireMatrixNear(fixed_ratio_gp, segment_gp, 5.0e-7, 5.0e-7,
                      "fixed-ratio point gradient");
    requireNear(fixed_ratio_gtotal, contracted_time_grad,
                2.0e-4, 2.0e-5, "fixed-ratio contracted time gradient");

    double plus_delta = 0.0;
    double minus_delta = 0.0;
    externalTotalDurationDeltas(ratios, total_duration,
                                plus_delta, minus_delta);
    const double plus_cost =
        fixedRatioEnergy<Dim, S>(data.inner_points, data.head_state,
                                 data.tail_state, ratios,
                                 total_duration + plus_delta);
    const double minus_cost =
        fixedRatioEnergy<Dim, S>(data.inner_points, data.head_state,
                                 data.tail_state, ratios,
                                 total_duration + minus_delta);
    const double numeric_gtotal =
        (plus_cost - minus_cost) / (plus_delta - minus_delta);

    requireNear(fixed_ratio_gtotal, numeric_gtotal,
                3.0e-3, 5.0e-4, "fixed-ratio total-duration gradient");
}

inline void runFixedRatioTotalDurationCase()
{
    std::cout << "[fixed-ratio total-duration gradient]" << std::endl;
    checkFixedRatioTotalDurationOne<1, 2>(5201);
    checkFixedRatioTotalDurationOne<2, 2>(5202);
    checkFixedRatioTotalDurationOne<3, 2>(5203);
    checkFixedRatioTotalDurationOne<1, 3>(5301);
    checkFixedRatioTotalDurationOne<2, 3>(5302);
    checkFixedRatioTotalDurationOne<3, 3>(5303);
    checkFixedRatioTotalDurationOne<1, 4>(5401);
    checkFixedRatioTotalDurationOne<2, 4>(5402);
    checkFixedRatioTotalDurationOne<3, 4>(5403);
    std::cout << "  PASS" << std::endl;
}

template <int Dim, int S>
inline double uniformTimeEnergy(
    const Eigen::Matrix<double, Eigen::Dynamic, Dim> &inner_points,
    const Eigen::Matrix<double, Dim, Eigen::Dynamic> &head_state,
    const Eigen::Matrix<double, Dim, Eigen::Dynamic> &tail_state,
    const double total_duration)
{
    nubs::UniformNUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd control_points;
    trajectory.generateUniform(inner_points, head_state, tail_state,
                               total_duration, control_points);
    return trajectory.getEnergy();
}

template <int Dim, int S>
inline void checkUniformTimeOne(const unsigned int seed)
{
    constexpr int M = 5;
    const auto data = makeRandomProblem<Dim>(S, M, seed);
    const double total_duration = data.durations.sum();
    const Eigen::VectorXd uniform_durations =
        nubs::UniformNUBSTrajectoryT<Dim, S>::uniformDurations(M,
                                                               total_duration);

    nubs::UniformNUBSTrajectoryT<Dim, S> uniform;
    Eigen::MatrixXd uniform_control;
    uniform.generateUniform(data.inner_points, data.head_state,
                            data.tail_state, total_duration,
                            uniform_control);

    nubs::NUBSTrajectoryT<Dim, S> reference;
    Eigen::MatrixXd reference_control;
    reference.generate(data.inner_points, data.head_state, data.tail_state,
                       uniform_durations, reference_control);

    nubs::UniformNUBSTrajectory<Dim, 2 * S - 1> generic(S);
    Eigen::MatrixXd generic_control;
    generic.generateUniform(data.inner_points, data.head_state,
                            data.tail_state, total_duration,
                            generic_control);

    requireMatrixNear(uniform.getDurations(), uniform_durations,
                      1.0e-12, 1.0e-12, "uniform durations");
    requireMatrixNear(uniform_control, reference_control,
                      2.0e-10, 2.0e-10, "uniform fixed control");
    requireMatrixNear(generic_control, uniform_control,
                      2.0e-10, 2.0e-10, "uniform generic control");
    requireNear(uniform.getEnergy(), reference.getEnergy(),
                2.0e-9, 2.0e-9, "uniform fixed energy");
    requireNear(generic.getEnergy(), uniform.getEnergy(),
                2.0e-9, 2.0e-9, "uniform generic energy");

    const double dt = total_duration / static_cast<double>(M);
    for (int i = 1; i < M; ++i)
    {
        const auto value = uniform.evaluate(dt * static_cast<double>(i), 0);
        const Eigen::Matrix<double, Dim, 1> expected =
            data.inner_points.row(i - 1).transpose();
        requireMatrixNear(value, expected, 1.0e-7, 1.0e-7,
                          "uniform waypoint");
    }

    for (int i = 0; i <= 20; ++i)
    {
        double t = total_duration * static_cast<double>(i) / 20.0;
        if (i == 20)
        {
            t = previousTime(total_duration);
        }

        Eigen::Matrix<double, Dim, 1> pos;
        Eigen::Matrix<double, Dim, 1> vel;
        Eigen::Matrix<double, Dim, 1> acc;
        Eigen::Matrix<double, Dim, 1> jerk;
        Eigen::Matrix<double, Dim, 1> snap;
        uniform.evaluatePVAJS(t, pos, vel, acc, jerk, snap);

        requireMatrixNear(pos, uniform.evaluate(t, 0),
                          1.0e-8, 1.0e-8, "uniform PVAJS position");
        requireMatrixNear(vel, uniform.evaluate(t, 1),
                          1.0e-8, 1.0e-8, "uniform PVAJS velocity");
        requireMatrixNear(acc, uniform.evaluate(t, 2),
                          1.0e-8, 1.0e-8, "uniform PVAJS acceleration");
        requireMatrixNear(jerk, uniform.evaluate(t, 3),
                          1.0e-8, 1.0e-8, "uniform PVAJS jerk");
        requireMatrixNear(snap, uniform.evaluate(t, 4),
                          1.0e-8, 1.0e-8, "uniform PVAJS snap");

        requireMatrixNear(uniform.getPos(t), pos,
                          1.0e-8, 1.0e-8, "uniform getPos");
        requireMatrixNear(uniform.getVel(t), vel,
                          1.0e-8, 1.0e-8, "uniform getVel");
        requireMatrixNear(uniform.getAcc(t), acc,
                          1.0e-8, 1.0e-8, "uniform getAcc");
        requireMatrixNear(uniform.getJer(t), jerk,
                          1.0e-8, 1.0e-8, "uniform getJer");
        requireMatrixNear(uniform.getJerk(t), jerk,
                          1.0e-8, 1.0e-8, "uniform getJerk");
        requireMatrixNear(uniform.getSnap(t), snap,
                          1.0e-8, 1.0e-8, "uniform getSnap");
    }

    double uniform_cost = 0.0;
    Eigen::MatrixXd uniform_gp;
    double uniform_gtotal = 0.0;
    uniform.getEnergyAndUniformTimeGrad(uniform_cost, uniform_gp,
                                        uniform_gtotal);

    double segment_cost = 0.0;
    Eigen::MatrixXd segment_gp;
    Eigen::VectorXd segment_gt;
    uniform.getEnergyAndGrad(segment_cost, segment_gp, segment_gt);

    requireNear(uniform_cost, segment_cost, 2.0e-9, 2.0e-9,
                "uniform gradient cost");
    requireMatrixNear(uniform_gp, segment_gp, 5.0e-7, 5.0e-7,
                      "uniform point gradient");
    requireNear(uniform_gtotal, segment_gt.sum() / static_cast<double>(M),
                2.0e-4, 2.0e-5, "uniform contracted time gradient");

    double generic_cost = 0.0;
    Eigen::MatrixXd generic_gp;
    double generic_gtotal = 0.0;
    generic.getEnergyAndUniformTimeGrad(generic_cost, generic_gp,
                                        generic_gtotal);
    requireNear(generic_cost, uniform_cost, 2.0e-9, 2.0e-9,
                "uniform generic gradient cost");
    requireMatrixNear(generic_gp, uniform_gp, 5.0e-7, 5.0e-7,
                      "uniform generic point gradient");
    requireNear(generic_gtotal, uniform_gtotal, 2.0e-5, 2.0e-5,
                "uniform generic total-duration gradient");

    Eigen::Matrix<double, Eigen::Dynamic, Dim> numeric_gp =
        Eigen::Matrix<double, Eigen::Dynamic, Dim>::Zero(M - 1, Dim);
    const double point_eps = 1.0e-6;
    for (int r = 0; r < data.inner_points.rows(); ++r)
    {
        for (int c = 0; c < Dim; ++c)
        {
            auto plus_points = data.inner_points;
            auto minus_points = data.inner_points;
            plus_points(r, c) += point_eps;
            minus_points(r, c) -= point_eps;
            const double plus_cost =
                uniformTimeEnergy<Dim, S>(plus_points, data.head_state,
                                          data.tail_state, total_duration);
            const double minus_cost =
                uniformTimeEnergy<Dim, S>(minus_points, data.head_state,
                                          data.tail_state, total_duration);
            numeric_gp(r, c) =
                (plus_cost - minus_cost) / (2.0 * point_eps);
        }
    }

    Eigen::VectorXd ones = Eigen::VectorXd::Ones(M);
    double plus_delta = 0.0;
    double minus_delta = 0.0;
    externalTotalDurationDeltas(ones, total_duration,
                                plus_delta, minus_delta);
    const double plus_cost =
        uniformTimeEnergy<Dim, S>(data.inner_points, data.head_state,
                                  data.tail_state,
                                  total_duration + plus_delta);
    const double minus_cost =
        uniformTimeEnergy<Dim, S>(data.inner_points, data.head_state,
                                  data.tail_state,
                                  total_duration + minus_delta);
    const double numeric_gtotal =
        (plus_cost - minus_cost) / (plus_delta - minus_delta);

    requireMatrixNear(uniform_gp, numeric_gp, 5.0e-4, 2.0e-4,
                      "uniform external point gradient");
    requireNear(uniform_gtotal, numeric_gtotal,
                3.0e-3, 5.0e-4, "uniform external total-duration gradient");
}

inline void runUniformTimeCase()
{
    std::cout << "[uniform-time NUBS]" << std::endl;
    checkUniformTimeOne<1, 2>(6201);
    checkUniformTimeOne<2, 2>(6202);
    checkUniformTimeOne<3, 2>(6203);
    checkUniformTimeOne<1, 3>(6301);
    checkUniformTimeOne<2, 3>(6302);
    checkUniformTimeOne<3, 3>(6303);
    checkUniformTimeOne<1, 4>(6401);
    checkUniformTimeOne<2, 4>(6402);
    checkUniformTimeOne<3, 4>(6403);
    std::cout << "  PASS" << std::endl;
}

template <int Dim, int S>
inline void checkLocalVsFullOne(const unsigned int seed)
{
    constexpr int M = 6;
    const auto data = makeRandomProblem<Dim>(S, M, seed);

    nubs::NUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd control_points;
    trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                        data.durations, control_points);

    double local_cost = 0.0;
    double full_cost = 0.0;
    Eigen::MatrixXd local_gp;
    Eigen::MatrixXd full_gp;
    Eigen::VectorXd local_gt;
    Eigen::VectorXd full_gt;
    trajectory.getEnergyAndFiniteDiffGrad(local_cost, local_gp, local_gt);
    trajectory.getEnergyAndFiniteDiffGradFull(full_cost, full_gp, full_gt);

    requireNear(local_cost, full_cost, 1.0e-10, 1.0e-10,
                "local vs full cost");
    requireMatrixNear(local_gp, full_gp, 1.0e-10, 1.0e-10,
                      "local vs full point gradient");
    requireMatrixNear(local_gt, full_gt, 2.0e-6, 2.0e-8,
                      "local vs full time gradient");
}

inline void runLocalVsFullFiniteDiffCase()
{
    std::cout << "[local vs full finite-diff]" << std::endl;
    checkLocalVsFullOne<3, 2>(4201);
    checkLocalVsFullOne<3, 3>(4301);
    checkLocalVsFullOne<3, 4>(4401);
    std::cout << "  PASS" << std::endl;
}

template <int S>
inline void checkSmallDurationOne()
{
    constexpr int Dim = 3;
    constexpr int M = 5;
    auto data = makeRandomProblem<Dim>(S, M, 5000 + S, true);

    nubs::NUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd control_points;
    trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                        data.durations, control_points);

    double cost = 0.0;
    Eigen::MatrixXd grad_points;
    Eigen::VectorXd grad_times;
    trajectory.getEnergyAndFiniteDiffGrad(cost, grad_points, grad_times);
    require(std::isfinite(cost), "small duration cost is not finite");
    requireFiniteMatrix(grad_points, "small duration point gradient is not finite");
    requireFiniteMatrix(grad_times, "small duration time gradient is not finite");
}

inline void runSmallDurationRobustnessCase()
{
    std::cout << "[small duration robustness]" << std::endl;
    checkSmallDurationOne<2>();
    checkSmallDurationOne<3>();
    checkSmallDurationOne<4>();
    std::cout << "  PASS" << std::endl;
}

inline void runFiniteDiffBenchmarkS3D3()
{
    constexpr int Dim = 3;
    constexpr int S = 3;
    constexpr int runs = 1000;
    std::cout << "[finite-diff benchmark S=3 Dim=3]" << std::endl;

    for (const int M : {4, 8, 16, 32})
    {
        const auto data = makeRandomProblem<Dim>(S, M, 6300 + M);
        nubs::QuinticNUBS<Dim> trajectory;
        Eigen::MatrixXd control_points;
        trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                            data.durations, control_points);

        const auto generate_time = measureRepeated(
            [&]()
            {
                nubs::QuinticNUBS<Dim> tmp;
                Eigen::MatrixXd out;
                tmp.generate(data.inner_points, data.head_state, data.tail_state,
                             data.durations, out);
            },
            runs);

        const auto coeff_time = measureRepeated(
            [&]()
            {
                double cost = 0.0;
                Eigen::MatrixXd gdC;
                trajectory.getEnergyPartialGradByCoeffs(cost, gdC);
                volatile double sink = cost + gdC.squaredNorm();
                (void)sink;
            },
            runs);

        const auto local_time = measureRepeated(
            [&]()
            {
                double cost = 0.0;
                Eigen::MatrixXd gp;
                Eigen::VectorXd gt;
                trajectory.getEnergyAndFiniteDiffGrad(cost, gp, gt);
                volatile double sink = cost + gp.squaredNorm() + gt.squaredNorm();
                (void)sink;
            },
            runs);

        const auto full_time = measureRepeated(
            [&]()
            {
                double cost = 0.0;
                Eigen::MatrixXd gp;
                Eigen::VectorXd gt;
                trajectory.getEnergyAndFiniteDiffGradFull(cost, gp, gt);
                volatile double sink = cost + gp.squaredNorm() + gt.squaredNorm();
                (void)sink;
            },
            runs);

        std::cout << std::fixed << std::setprecision(3)
                  << "  M=" << M
                  << ": generate=" << generate_time.avg_us << " us"
                  << ", coeff=" << coeff_time.avg_us << " us"
                  << ", local_fd=" << local_time.avg_us << " us"
                  << ", full_fd=" << full_time.avg_us << " us"
                  << ", speedup=" << full_time.avg_us / local_time.avg_us
                  << "x" << std::endl;

        if (M >= 8)
        {
            require(local_time.avg_us < full_time.avg_us,
                    "optimized local finite-diff path should beat full path");
        }
    }
    std::cout << "  PASS" << std::endl;
}

} // namespace nubs_test

#endif
