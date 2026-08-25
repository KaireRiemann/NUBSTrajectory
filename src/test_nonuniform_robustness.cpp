#include "tools/optimization_test_cases.hpp"

#include <exception>
#include <iostream>

namespace
{

template <int S>
void checkRatio(const double ratio, const unsigned int seed)
{
    constexpr int Dim = 3;
    constexpr int M = 6;
    auto data = nubs_test::makeRandomProblem<Dim>(S, M, seed);
    for (int i = 0; i < M; ++i)
    {
        // Alternating short and long segments exercises both sides of every
        // local timing stencil while keeping the absolute scale well posed.
        data.durations(i) = 0.02 * (i % 2 == 0 ? 1.0 : ratio);
    }

    nubs::NUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd controls;
    trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                        data.durations, controls);
    nubs_test::requireFiniteMatrix(controls, "nonuniform controls are not finite");
    nubs_test::require(trajectory.A.factorizationStats().near_zero_pivot_count == 0,
                       "unexpected near-zero pivot");
    nubs_test::require(std::isfinite(trajectory.A.factorizationStats().minimum_abs_pivot) &&
                           trajectory.A.factorizationStats().minimum_abs_pivot > 0.0,
                       "invalid pivot statistic");
    nubs_test::require(std::isfinite(trajectory.getLastLinearSolveRelativeResidual()) &&
                           trajectory.getLastLinearSolveRelativeResidual() < 1.0e-8,
                       "linear-solve residual is too large");

    double cost = 0.0;
    Eigen::MatrixXd grad_points;
    Eigen::VectorXd grad_times;
    trajectory.getEnergyAndGrad(cost, grad_points, grad_times);
    nubs_test::require(std::isfinite(cost), "nonuniform energy is not finite");
    nubs_test::requireFiniteMatrix(grad_points,
                                   "nonuniform point gradient is not finite");
    nubs_test::requireFiniteMatrix(grad_times,
                                   "nonuniform time gradient is not finite");

    // Construction interpolation is the relevant end-to-end residual for this
    // standalone repository: all boundary derivatives and interior waypoints
    // must survive timing-ratio stress.
    for (int waypoint = 1; waypoint < M; ++waypoint)
    {
        nubs_test::requireMatrixNear(
            trajectory.evaluate(nubs_test::cumulativeTime(data.durations, waypoint)),
            data.inner_points.row(waypoint - 1).transpose(),
            1.0e-5, 1.0e-6, "nonuniform waypoint residual");
    }
    std::cout << "  S=" << S << ", ratio=" << ratio
              << ", min_pivot=" << trajectory.A.factorizationStats().minimum_abs_pivot
              << ", residual=" << trajectory.getLastLinearSolveRelativeResidual()
              << ", energy=" << cost << std::endl;
}

template <int S>
void runOrder()
{
    checkRatio<S>(1.0, 5100 + S);
    checkRatio<S>(10.0, 5200 + S);
    checkRatio<S>(100.0, 5300 + S);
    checkRatio<S>(1000.0, 5400 + S);
}

} // namespace

int main()
{
    try
    {
        std::cout << "[non-uniform timing robustness]" << std::endl;
        runOrder<2>();
        runOrder<3>();
        runOrder<4>();
        std::cout << "  PASS" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_nonuniform_robustness failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
