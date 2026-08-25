#ifndef NUBS_TEST_CASES_HPP
#define NUBS_TEST_CASES_HPP

#include "tools/minco_adapter.hpp"
#include "tools/test_common.hpp"

namespace nubs_test
{

template <int Dim, int S>
inline void runBasicConstructionCase(const int piece_num = 5)
{
    const auto data = makeProblem<Dim>(S, piece_num);

    nubs::NUBSTrajectoryT<Dim, S> fixed_trajectory;
    nubs::NUBSTrajectory<Dim> generic_trajectory(S);
    Eigen::MatrixXd fixed_control_points;
    Eigen::MatrixXd generic_control_points;

    fixed_trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                              data.durations, fixed_control_points);
    generic_trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                                data.durations, generic_control_points);

    const int degree = 2 * S - 1;
    const int expected_ctrl_pts = piece_num + 2 * S - 1;
    require(fixed_trajectory.getS() == S, "fixed sys_order mismatch");
    require(fixed_trajectory.getP() == degree, "fixed degree mismatch");
    require(fixed_control_points.rows() == expected_ctrl_pts &&
                fixed_control_points.cols() == Dim,
            "fixed control point shape mismatch");
    requireMatrixNear(fixed_control_points, generic_control_points,
                      1.0e-12, 1.0e-12, "fixed vs generic control points");

    const Eigen::VectorXd &knots = fixed_trajectory.getKnots();
    double max_knot_error = 0.0;
    double max_boundary_error = 0.0;
    double max_waypoint_error = 0.0;

    for (int i = 0; i <= degree; ++i)
    {
        max_knot_error = std::max(max_knot_error, std::abs(knots(i)));
        max_knot_error =
            std::max(max_knot_error,
                     std::abs(knots(knots.size() - 1 - i) - data.durations.sum()));
    }
    for (int i = 1; i <= piece_num; ++i)
    {
        max_knot_error =
            std::max(max_knot_error,
                     std::abs(knots(degree + i) - cumulativeTime(data.durations, i)));
    }

    for (int der = 0; der < S; ++der)
    {
        const Eigen::Matrix<double, Dim, 1> head_value =
            fixed_trajectory.evaluate(0.0, der);
        const Eigen::Matrix<double, Dim, 1> tail_value =
            fixed_trajectory.evaluate(previousTime(data.durations.sum()), der);
        max_boundary_error =
            std::max(max_boundary_error,
                     (head_value - data.head_state.col(der)).norm());
        max_boundary_error =
            std::max(max_boundary_error,
                     (tail_value - data.tail_state.col(der)).norm());
        requireMatrixNear(head_value,
                          data.head_state.col(der),
                          2.0e-8, 2.0e-8, "head boundary");
        requireMatrixNear(tail_value,
                          data.tail_state.col(der),
                          2.0e-5, 2.0e-6, "tail boundary");
    }

    for (int i = 1; i < piece_num; ++i)
    {
        const double t = cumulativeTime(data.durations, i);
        const Eigen::Matrix<double, Dim, 1> value =
            fixed_trajectory.evaluate(t, 0);
        const Eigen::Matrix<double, Dim, 1> expected =
            data.inner_points.row(i - 1).transpose();
        const Eigen::Matrix<double, Dim, 1> diff = value - expected;
        max_waypoint_error = std::max(max_waypoint_error, diff.norm());
        requireMatrixNear(value,
                          expected,
                          2.0e-8, 2.0e-8, "inner waypoint");
    }

    constexpr int timing_runs = 40;
    const auto fixed_build = measureRepeated(
        [&]()
        {
            nubs::NUBSTrajectoryT<Dim, S> traj;
            Eigen::MatrixXd out;
            traj.generate(data.inner_points, data.head_state, data.tail_state,
                          data.durations, out);
        },
        timing_runs);
    const auto generic_build = measureRepeated(
        [&]()
        {
            nubs::NUBSTrajectory<Dim> traj(S);
            Eigen::MatrixXd out;
            traj.generate(data.inner_points, data.head_state, data.tail_state,
                          data.durations, out);
        },
        timing_runs);

    std::cout << "[" << orderName(S) << " basic] Dim=" << Dim
              << ", M=" << piece_num
              << ", ctrl=" << fixed_control_points.rows()
              << ", knots=" << knots.size() << std::endl;
    std::cout << std::scientific << std::setprecision(4)
              << "  errors: knot=" << max_knot_error
              << ", boundary=" << max_boundary_error
              << ", waypoint=" << max_waypoint_error
              << ", fixed_vs_generic_ctrl="
              << maxAbsCoeffDiff(fixed_control_points, generic_control_points)
              << std::endl;
    std::cout << std::fixed << std::setprecision(3)
              << "  build avg: fixed=" << fixed_build.avg_us
              << " us, generic=" << generic_build.avg_us << " us" << std::endl;
    std::cout << "  PASS" << std::endl;
}

template <int S>
inline void runMincoComparison3DCase(const int piece_num = 6,
                                     const int sample_num = 120)
{
    constexpr int Dim = 3;
    const auto data = makeProblem<Dim>(S, piece_num);

    constexpr int timing_runs = 40;
    const auto nubs_build = measureRepeated(
        [&]()
        {
            nubs::NUBSTrajectoryT<Dim, S> traj;
            Eigen::MatrixXd out;
            traj.generate(data.inner_points, data.head_state, data.tail_state,
                          data.durations, out);
        },
        timing_runs);
    const auto minco_build = measureRepeated(
        [&]()
        {
            SplineMincoAdapter<Dim> ref;
            ref.generate(data.inner_points, data.head_state, data.tail_state,
                         data.durations, S);
        },
        timing_runs);

    nubs::NUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd control_points;
    trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                        data.durations, control_points);

    SplineMincoAdapter<Dim> reference;
    reference.generate(data.inner_points, data.head_state, data.tail_state,
                       data.durations, S);

    std::array<ErrorStats, S + 1> derivative_errors;
    std::vector<double> sample_times(sample_num + 1);
    const double total_t = data.durations.sum();
    for (int sample = 0; sample <= sample_num; ++sample)
    {
        double t = total_t * static_cast<double>(sample) /
                   static_cast<double>(sample_num);
        if (sample == sample_num)
        {
            t = previousTime(total_t);
        }
        sample_times[sample] = t;

        for (int der = 0; der <= S; ++der)
        {
            const auto nubs_value = trajectory.evaluate(t, der);
            const auto minco_value = reference.evaluate(t, der);
            accumulateSampleError(derivative_errors[der],
                                  nubs_value - minco_value);
            requireMatrixNear(nubs_value, minco_value,
                              5.0e-7, 5.0e-7, "NUBS vs MINCO");
        }
    }
    for (auto &stats : derivative_errors)
    {
        finalizeErrorStats(stats);
    }

    const double nubs_energy = trajectory.getEnergy();
    const double minco_energy = reference.getEnergy();
    requireNear(nubs_energy, minco_energy,
                2.0e-7, 2.0e-7, "NUBS vs MINCO energy");

    const auto nubs_query = measureQueries(
        [&](const int i)
        {
            volatile double sink = trajectory.evaluate(sample_times[i], 0).squaredNorm();
            (void)sink;
        },
        static_cast<int>(sample_times.size()));
    const auto minco_query = measureQueries(
        [&](const int i)
        {
            volatile double sink = reference.evaluate(sample_times[i], 0).squaredNorm();
            (void)sink;
        },
        static_cast<int>(sample_times.size()));

    std::cout << "[" << orderName(S) << " vs MINCO 3D] M=" << piece_num
              << ", samples=" << sample_times.size()
              << ", ctrl=" << control_points.rows() << std::endl;
    std::cout << std::scientific << std::setprecision(4) << "  derivative max:";
    for (int der = 0; der <= S; ++der)
    {
        std::cout << " d" << der << "=" << derivative_errors[der].max_norm;
    }
    std::cout << std::endl;
    std::cout << std::scientific << std::setprecision(8)
              << "  energy: nubs=" << nubs_energy
              << ", minco=" << minco_energy
              << ", abs_err=" << std::abs(nubs_energy - minco_energy)
              << ", rel_err=" << relativeError(nubs_energy, minco_energy)
              << std::endl;
    std::cout << std::fixed << std::setprecision(3)
              << "  build avg: nubs=" << nubs_build.avg_us
              << " us, minco=" << minco_build.avg_us << " us" << std::endl;
    std::cout << std::fixed << std::setprecision(3)
              << "  query avg: nubs=" << nubs_query.avg_us
              << " us, minco=" << minco_query.avg_us << " us" << std::endl;
    std::cout << "  PASS" << std::endl;
}

inline double factorialValue(const int n)
{
    double value = 1.0;
    for (int i = 2; i <= n; ++i)
    {
        value *= static_cast<double>(i);
    }
    return value;
}

inline double fallingFactorial(const int n, const int k)
{
    double value = 1.0;
    for (int i = 0; i < k; ++i)
    {
        value *= static_cast<double>(n - i);
    }
    return value;
}

template <int Dim, int S>
inline Eigen::Matrix<double, Eigen::Dynamic, Dim> extractLocalPolynomialCoeffs(
    const nubs::NUBSTrajectoryT<Dim, S> &trajectory,
    const Eigen::VectorXd &durations)
{
    constexpr int Degree = 2 * S - 1;
    constexpr int CoeffPerPiece = Degree + 1;
    const int piece_num = static_cast<int>(durations.size());
    Eigen::Matrix<double, Eigen::Dynamic, Dim> coeffs(piece_num * CoeffPerPiece,
                                                      Dim);

    for (int piece = 0; piece < piece_num; ++piece)
    {
        const double t0 = cumulativeTime(durations, piece);
        for (int k = 0; k <= Degree; ++k)
        {
            coeffs.row(piece * CoeffPerPiece + k) =
                (trajectory.evaluate(t0, k) / factorialValue(k)).transpose();
        }
    }
    return coeffs;
}

template <int Dim, int S>
inline void getLocalPolynomialEnergyGrad(
    const Eigen::Matrix<double, Eigen::Dynamic, Dim> &coeffs,
    const Eigen::VectorXd &durations,
    double &cost,
    Eigen::Matrix<double, Eigen::Dynamic, Dim> &gradByCoeffs,
    Eigen::VectorXd &partialGradByTimes)
{
    constexpr int Degree = 2 * S - 1;
    constexpr int CoeffPerPiece = Degree + 1;
    using Gauss = nubs::FixedGaussRule<S>;

    cost = 0.0;
    gradByCoeffs.setZero(coeffs.rows(), Dim);
    partialGradByTimes.setZero(durations.size());

    for (int piece = 0; piece < durations.size(); ++piece)
    {
        const double T = durations(piece);
        for (int q = 0; q < Gauss::Num; ++q)
        {
            const double tau = 0.5 * T * (Gauss::nodes[q] + 1.0);
            const double w = 0.5 * T * Gauss::weights[q];
            Eigen::Matrix<double, Dim, 1> value =
                Eigen::Matrix<double, Dim, 1>::Zero();
            double basis[CoeffPerPiece] = {0.0};

            for (int k = S; k <= Degree; ++k)
            {
                basis[k] = fallingFactorial(k, S) *
                           std::pow(tau, static_cast<double>(k - S));
                value += basis[k] *
                         coeffs.row(piece * CoeffPerPiece + k).transpose();
            }

            cost += w * value.squaredNorm();
            for (int k = S; k <= Degree; ++k)
            {
                gradByCoeffs.row(piece * CoeffPerPiece + k) +=
                    2.0 * w * basis[k] * value.transpose();
            }
        }

        Eigen::Matrix<double, Dim, 1> terminal_value =
            Eigen::Matrix<double, Dim, 1>::Zero();
        for (int k = S; k <= Degree; ++k)
        {
            const double basis =
                fallingFactorial(k, S) *
                std::pow(T, static_cast<double>(k - S));
            terminal_value += basis *
                              coeffs.row(piece * CoeffPerPiece + k).transpose();
        }
        partialGradByTimes(piece) = terminal_value.squaredNorm();
    }
}

template <int S>
inline void runMincoEnergyGradient3DCase(const int piece_num = 6)
{
    constexpr int Dim = 3;
    constexpr int Degree = 2 * S - 1;
    constexpr int CoeffPerPiece = Degree + 1;
    const auto data = makeProblem<Dim>(S, piece_num);

    nubs::NUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd control_points;
    trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                        data.durations, control_points);

    SplineMincoAdapter<Dim> reference;
    reference.generate(data.inner_points, data.head_state, data.tail_state,
                       data.durations, S);

    const auto nubs_poly_coeffs =
        extractLocalPolynomialCoeffs<Dim, S>(trajectory, data.durations);
    const Eigen::MatrixXd minco_coeffs = reference.getCoeffs();

    double nubs_poly_cost = 0.0;
    Eigen::Matrix<double, Eigen::Dynamic, Dim> nubs_poly_gdC;
    Eigen::VectorXd nubs_partial_gdT;
    getLocalPolynomialEnergyGrad<Dim, S>(nubs_poly_coeffs, data.durations,
                                         nubs_poly_cost, nubs_poly_gdC,
                                         nubs_partial_gdT);

    double nubs_cost = 0.0;
    Eigen::MatrixXd nubs_grad_points;
    Eigen::VectorXd nubs_grad_times;
    trajectory.getEnergyAndGrad(nubs_cost, nubs_grad_points,
                                nubs_grad_times);

    double minco_cost = 0.0;
    Eigen::MatrixXd minco_gdC;
    Eigen::Matrix<double, Eigen::Dynamic, Dim> minco_grad_points;
    Eigen::VectorXd minco_grad_times;
    Eigen::VectorXd minco_partial_gdT;
    reference.getEnergyAndGrad(minco_cost, minco_gdC, minco_grad_points,
                               minco_grad_times, minco_partial_gdT);

    requireMatrixNear(nubs_poly_coeffs, minco_coeffs,
                      2.0e-7, 2.0e-8, "NUBS local coeffs vs MINCO coeffs");
    requireNear(nubs_cost, minco_cost,
                2.0e-7, 2.0e-8, "NUBS vs MINCO gradient cost");
    requireNear(nubs_poly_cost, minco_cost,
                2.0e-7, 2.0e-8, "NUBS poly cost vs MINCO cost");
    requireMatrixNear(nubs_poly_gdC, minco_gdC,
                      2.0e-5, 2.0e-7, "NUBS poly coeff grad vs MINCO");
    requireMatrixNear(nubs_partial_gdT, minco_partial_gdT,
                      2.0e-5, 2.0e-7, "NUBS direct time grad vs MINCO");
    requireMatrixNear(nubs_grad_points, minco_grad_points,
                      2.0e-5, 2.0e-7, "NUBS propagated point grad vs MINCO");
    requireMatrixNear(nubs_grad_times, minco_grad_times,
                      5.0e-3, 5.0e-5, "NUBS propagated time grad vs MINCO");

    constexpr int timing_runs = 40;
    const auto nubs_grad_time = measureRepeated(
        [&]()
        {
            double cost = 0.0;
            Eigen::MatrixXd gp;
            Eigen::VectorXd gt;
            trajectory.getEnergyAndGrad(cost, gp, gt);
        },
        timing_runs);
    const auto minco_grad_time = measureRepeated(
        [&]()
        {
            double cost = 0.0;
            Eigen::MatrixXd gdC;
            Eigen::Matrix<double, Eigen::Dynamic, Dim> gp;
            Eigen::VectorXd gt;
            Eigen::VectorXd gtd;
            reference.getEnergyAndGrad(cost, gdC, gp, gt, gtd);
        },
        timing_runs);

    std::cout << "[" << orderName(S) << " MINCO gradient 3D] M=" << piece_num
              << ", coeff_vars=" << piece_num * CoeffPerPiece * Dim
              << ", point_vars=" << nubs_grad_points.size()
              << ", time_vars=" << nubs_grad_times.size() << std::endl;
    std::cout << std::scientific << std::setprecision(4)
              << "  coeff: value_max="
              << maxAbsCoeffDiff(nubs_poly_coeffs, minco_coeffs)
              << ", grad_max=" << maxAbsCoeffDiff(nubs_poly_gdC, minco_gdC)
              << ", direct_time_max="
              << maxAbsCoeffDiff(nubs_partial_gdT, minco_partial_gdT)
              << std::endl;
    std::cout << std::scientific << std::setprecision(4)
              << "  propagated: point_max="
              << maxAbsCoeffDiff(nubs_grad_points, minco_grad_points)
              << ", time_max="
              << maxAbsCoeffDiff(nubs_grad_times, minco_grad_times)
              << ", cost_abs=" << std::abs(nubs_cost - minco_cost)
              << std::endl;
    std::cout << std::fixed << std::setprecision(3)
              << "  grad avg: nubs=" << nubs_grad_time.avg_us
              << " us, minco=" << minco_grad_time.avg_us << " us" << std::endl;
    std::cout << "  PASS" << std::endl;
}

template <int Dim, int S>
inline void runCenteredGradientCase(const int piece_num = 5)
{
    const auto data = makeProblem<Dim>(S, piece_num);

    nubs::NUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd control_points;
    trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                        data.durations, control_points);

    constexpr int timing_runs = 20;
    const auto centered_grad_time = measureRepeated(
        [&]()
        {
            double local_cost = 0.0;
            Eigen::MatrixXd local_gp;
            Eigen::VectorXd local_gt;
            trajectory.getEnergyAndFiniteDiffGrad(local_cost, local_gp, local_gt);
        },
        timing_runs);

    double centered_cost = 0.0;
    Eigen::MatrixXd centered_grad_points;
    Eigen::VectorXd centered_grad_times;
    trajectory.getEnergyAndFiniteDiffGrad(centered_cost,
                                          centered_grad_points,
                                          centered_grad_times);

    Eigen::Matrix<double, Eigen::Dynamic, Dim> numeric_grad_points =
        Eigen::Matrix<double, Eigen::Dynamic, Dim>::Zero(piece_num - 1, Dim);
    Eigen::VectorXd numeric_grad_times = Eigen::VectorXd::Zero(piece_num);

    const double point_eps = 1.0e-6;
    const auto numeric_start = std::chrono::high_resolution_clock::now();
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
            numeric_grad_points(r, c) =
                (plus_cost - minus_cost) / (2.0 * point_eps);
        }
    }

    for (int i = 0; i < piece_num; ++i)
    {
        const double time_eps = 1.0e-6 * std::max(1.0, data.durations(i));
        Eigen::VectorXd plus_times = data.durations;
        Eigen::VectorXd minus_times = data.durations;
        plus_times(i) += time_eps;
        minus_times(i) -= time_eps;
        const double plus_cost =
            energyForFixed<Dim, S>(data.inner_points, data.head_state,
                                   data.tail_state, plus_times);
        const double minus_cost =
            energyForFixed<Dim, S>(data.inner_points, data.head_state,
                                   data.tail_state, minus_times);
        numeric_grad_times(i) =
            (plus_cost - minus_cost) / (2.0 * time_eps);
    }
    const auto numeric_end = std::chrono::high_resolution_clock::now();
    const double numeric_total_us =
        std::chrono::duration_cast<std::chrono::nanoseconds>(numeric_end - numeric_start).count() /
        1000.0;

    requireMatrixNear(centered_grad_points, numeric_grad_points,
                      2.0e-4, 2.0e-4, "centered point gradients");
    requireMatrixNear(centered_grad_times, numeric_grad_times,
                      3.0e-3, 5.0e-4, "centered time gradients");

    std::cout << "[" << orderName(S) << " centered gradient] Dim=" << Dim
              << ", M=" << piece_num
              << ", point_vars=" << centered_grad_points.size()
              << ", time_vars=" << centered_grad_times.size() << std::endl;
    std::cout << std::scientific << std::setprecision(8)
              << "  cost=" << centered_cost
              << ", point_max_abs_err="
              << maxAbsCoeffDiff(centered_grad_points, numeric_grad_points)
              << ", time_max_abs_err="
              << maxAbsCoeffDiff(centered_grad_times, numeric_grad_times)
              << std::endl;
    std::cout << std::fixed << std::setprecision(3)
              << "  reduced centered avg=" << centered_grad_time.avg_us
              << " us, full numeric total=" << numeric_total_us << " us" << std::endl;
    std::cout << "  PASS" << std::endl;
}

} // namespace nubs_test

#endif
