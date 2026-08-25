#include "tools/optimization_test_cases.hpp"

#include <exception>
#include <iostream>

namespace
{

template <int Dim, int S>
void checkLocalAD(const unsigned int seed)
{
    constexpr int M = 6;
    const auto data = nubs_test::makeRandomProblem<Dim>(S, M, seed);

    nubs::NUBSTrajectoryT<Dim, S> trajectory;
    Eigen::MatrixXd control_points;
    trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                        data.durations, control_points);

    double local_cost = 0.0;
    Eigen::MatrixXd local_points;
    Eigen::VectorXd local_times;
    trajectory.getEnergyAndGrad(local_cost, local_points, local_times);

    double dense_cost = 0.0;
    Eigen::MatrixXd dense_points;
    Eigen::VectorXd dense_times;
    trajectory.getEnergyAndAnalyticGrad(dense_cost, dense_points, dense_times);

    double finite_diff_cost = 0.0;
    Eigen::MatrixXd finite_diff_points;
    Eigen::VectorXd finite_diff_times;
    trajectory.getEnergyAndFiniteDiffGradFull(
        finite_diff_cost, finite_diff_points, finite_diff_times);

    nubs_test::requireNear(local_cost, dense_cost, 1.0e-9, 1.0e-10,
                           "Local-AD vs dense cost");
    nubs_test::requireMatrixNear(local_points, dense_points,
                                 2.0e-8, 2.0e-9,
                                 "Local-AD vs dense point gradient");
    nubs_test::requireMatrixNear(local_times, dense_times,
                                 5.0e-7, 5.0e-8,
                                 "Local-AD vs dense time gradient");
    nubs_test::requireNear(local_cost, finite_diff_cost, 1.0e-9, 1.0e-10,
                           "Local-AD vs finite-diff cost");
    nubs_test::requireMatrixNear(local_points, finite_diff_points,
                                 2.0e-8, 2.0e-9,
                                 "Local-AD vs finite-diff point gradient");
    nubs_test::requireMatrixNear(local_times, finite_diff_times,
                                 2.0e-3, 2.0e-5,
                                 "Local-AD vs finite-diff time gradient");

    // The runtime-order API is also part of the public contract.
    nubs::NUBSTrajectory<Dim, 2 * S - 1> runtime_trajectory(S);
    runtime_trajectory.generate(data.inner_points, data.head_state,
                                data.tail_state, data.durations,
                                control_points);
    double runtime_cost = 0.0;
    Eigen::MatrixXd runtime_points;
    Eigen::VectorXd runtime_times;
    runtime_trajectory.getEnergyAndGrad(runtime_cost, runtime_points,
                                         runtime_times);
    nubs_test::requireNear(runtime_cost, local_cost, 1.0e-9, 1.0e-10,
                           "runtime Local-AD cost");
    nubs_test::requireMatrixNear(runtime_points, local_points,
                                 2.0e-8, 2.0e-9,
                                 "runtime Local-AD point gradient");
    nubs_test::requireMatrixNear(runtime_times, local_times,
                                 5.0e-7, 5.0e-8,
                                 "runtime Local-AD time gradient");

    std::cout << "  S=" << S << ", Dim=" << Dim
              << ", local-vs-dense time max="
              << nubs_test::maxAbsCoeffDiff(local_times, dense_times)
              << std::endl;
}

} // namespace

int main()
{
    try
    {
        std::cout << "[Local-AD timing gradient]" << std::endl;
        checkLocalAD<1, 2>(4101);
        checkLocalAD<3, 2>(4103);
        checkLocalAD<1, 3>(4201);
        checkLocalAD<3, 3>(4203);
        checkLocalAD<1, 4>(4301);
        checkLocalAD<3, 4>(4303);
        std::cout << "  PASS" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_local_ad_gradient failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
