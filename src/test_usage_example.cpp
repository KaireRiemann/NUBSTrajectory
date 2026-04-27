#include "NUBSTrajectory.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>

int main()
{
    try
    {
        constexpr int Dim = 3;
        constexpr int S = 3; // minimum jerk, degree p = 5

        nubs::NUBSTrajectoryT<Dim, S> traj;

        Eigen::MatrixXd headState(Dim, S);
        Eigen::MatrixXd tailState(Dim, S);

        // columns: position, velocity, acceleration
        headState.col(0) = Eigen::Vector3d(0.0, 0.0, 0.0);
        headState.col(1) = Eigen::Vector3d(0.0, 0.0, 0.0);
        headState.col(2) = Eigen::Vector3d(0.0, 0.0, 0.0);

        tailState.col(0) = Eigen::Vector3d(5.0, 3.0, 1.0);
        tailState.col(1) = Eigen::Vector3d(0.0, 0.0, 0.0);
        tailState.col(2) = Eigen::Vector3d(0.0, 0.0, 0.0);

        Eigen::MatrixXd P_inner(2, Dim);
        P_inner.row(0) = Eigen::RowVector3d(1.5, 1.0, 0.5);
        P_inner.row(1) = Eigen::RowVector3d(3.5, 2.5, 0.8);

        Eigen::VectorXd T(3);
        T << 1.0, 1.2, 1.0;

        Eigen::MatrixXd control_points;
        traj.generate(P_inner, headState, tailState, T, control_points);

        const double t = 1.5;
        const Eigen::Vector3d pos = traj.evaluate(t, 0);
        const Eigen::Vector3d vel = traj.evaluate(t, 1);
        const Eigen::Vector3d acc = traj.evaluate(t, 2);
        const Eigen::Vector3d jerk = traj.evaluate(t, 3);
        const double energy = traj.getEnergy();

        if (!pos.allFinite() || !vel.allFinite() || !acc.allFinite() ||
            !jerk.allFinite() || !std::isfinite(energy))
        {
            throw std::runtime_error("usage example produced non-finite output");
        }

        std::cout << "[usage example] ctrl=" << control_points.rows()
                  << ", t=" << t
                  << ", pos=" << pos.transpose()
                  << ", vel_norm=" << vel.norm()
                  << ", acc_norm=" << acc.norm()
                  << ", jerk_norm=" << jerk.norm()
                  << ", energy=" << energy << std::endl;
        std::cout << "  PASS" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_usage_example failed: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
