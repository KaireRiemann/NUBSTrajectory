#ifndef NUBS_TEST_MINCO_ADAPTER_HPP
#define NUBS_TEST_MINCO_ADAPTER_HPP

#include "gcopter/minco.hpp"
#include "tools/test_common.hpp"

#include <vector>

namespace nubs_test
{

template <int Dim>
class SplineMincoAdapter
{
public:
    void generate(const Eigen::Matrix<double, Eigen::Dynamic, Dim> &inner_points,
                  const Eigen::Matrix<double, Dim, Eigen::Dynamic> &head_state,
                  const Eigen::Matrix<double, Dim, Eigen::Dynamic> &tail_state,
                  const Eigen::VectorXd &durations,
                  const int sys_order)
    {
        sys_order_ = sys_order;
        piece_num_ = static_cast<int>(durations.size());
        instance_num_ = (Dim + 2) / 3;
        durations_ = durations;

        if (sys_order_ == 2)
        {
            s2_solvers_.clear();
            s2_solvers_.resize(instance_num_);
            s2_trajs_.resize(instance_num_);
            for (int instance = 0; instance < instance_num_; ++instance)
            {
                Eigen::Matrix<double, 3, 2> head =
                    Eigen::Matrix<double, 3, 2>::Zero();
                Eigen::Matrix<double, 3, 2> tail =
                    Eigen::Matrix<double, 3, 2>::Zero();
                fillBoundary<2>(head_state, tail_state, instance, head, tail);

                s2_solvers_[instance].setConditions(head, tail, piece_num_);
                s2_solvers_[instance].setParameters(makeInnerPoints(inner_points, instance),
                                                    durations_);
                s2_solvers_[instance].getTrajectory(s2_trajs_[instance]);
            }
            return;
        }

        if (sys_order_ == 3)
        {
            s3_solvers_.clear();
            s3_solvers_.resize(instance_num_);
            s3_trajs_.resize(instance_num_);
            for (int instance = 0; instance < instance_num_; ++instance)
            {
                Eigen::Matrix3d head = Eigen::Matrix3d::Zero();
                Eigen::Matrix3d tail = Eigen::Matrix3d::Zero();
                fillBoundary<3>(head_state, tail_state, instance, head, tail);

                s3_solvers_[instance].setConditions(head, tail, piece_num_);
                s3_solvers_[instance].setParameters(makeInnerPoints(inner_points, instance),
                                                    durations_);
                s3_solvers_[instance].getTrajectory(s3_trajs_[instance]);
            }
            return;
        }

        if (sys_order_ == 4)
        {
            s4_solvers_.clear();
            s4_solvers_.resize(instance_num_);
            s4_trajs_.resize(instance_num_);
            for (int instance = 0; instance < instance_num_; ++instance)
            {
                Eigen::Matrix<double, 3, 4> head =
                    Eigen::Matrix<double, 3, 4>::Zero();
                Eigen::Matrix<double, 3, 4> tail =
                    Eigen::Matrix<double, 3, 4>::Zero();
                fillBoundary<4>(head_state, tail_state, instance, head, tail);

                s4_solvers_[instance].setConditions(head, tail, piece_num_);
                s4_solvers_[instance].setParameters(makeInnerPoints(inner_points, instance),
                                                    durations_);
                s4_solvers_[instance].getTrajectory(s4_trajs_[instance]);
            }
            return;
        }

        throw std::runtime_error("SplineMincoAdapter supports only sys_order 2, 3, and 4");
    }

    Eigen::Matrix<double, Dim, 1> evaluate(const double t, const int derivative) const
    {
        Eigen::Matrix<double, Dim, 1> value =
            Eigen::Matrix<double, Dim, 1>::Zero();
        for (int instance = 0; instance < instance_num_; ++instance)
        {
            const int start_dim = 3 * instance;
            const int active_dim = std::min(3, Dim - start_dim);
            const Eigen::Vector3d block_value = evaluate3(instance, t, derivative);
            for (int d = 0; d < active_dim; ++d)
            {
                value(start_dim + d) = block_value(d);
            }
        }
        return value;
    }

    double getEnergy() const
    {
        double energy = 0.0;
        double block_energy = 0.0;
        if (sys_order_ == 2)
        {
            for (const auto &solver : s2_solvers_)
            {
                solver.getEnergy(block_energy);
                energy += block_energy;
            }
        }
        else if (sys_order_ == 3)
        {
            for (const auto &solver : s3_solvers_)
            {
                solver.getEnergy(block_energy);
                energy += block_energy;
            }
        }
        else if (sys_order_ == 4)
        {
            for (const auto &solver : s4_solvers_)
            {
                solver.getEnergy(block_energy);
                energy += block_energy;
            }
        }
        return energy;
    }

    Eigen::MatrixXd getCoeffs() const
    {
        Eigen::MatrixXd coeffs =
            Eigen::MatrixXd::Zero(2 * sys_order_ * piece_num_, Dim);
        if (sys_order_ == 2)
        {
            copyCoeffBlocks(s2_solvers_, coeffs);
        }
        else if (sys_order_ == 3)
        {
            copyCoeffBlocks(s3_solvers_, coeffs);
        }
        else if (sys_order_ == 4)
        {
            copyCoeffBlocks(s4_solvers_, coeffs);
        }
        return coeffs;
    }

    void getEnergyPartialGradByCoeffs(Eigen::MatrixXd &gdC) const
    {
        gdC.setZero(2 * sys_order_ * piece_num_, Dim);
        if (sys_order_ == 2)
        {
            copyCoeffGradBlocks(s2_solvers_, gdC);
        }
        else if (sys_order_ == 3)
        {
            copyCoeffGradBlocks(s3_solvers_, gdC);
        }
        else if (sys_order_ == 4)
        {
            copyCoeffGradBlocks(s4_solvers_, gdC);
        }
    }

    void getEnergyPartialGradByTimes(Eigen::VectorXd &gdT) const
    {
        gdT = Eigen::VectorXd::Zero(piece_num_);
        if (sys_order_ == 2)
        {
            accumulateTimeGradBlocks(s2_solvers_, gdT);
        }
        else if (sys_order_ == 3)
        {
            accumulateTimeGradBlocks(s3_solvers_, gdT);
        }
        else if (sys_order_ == 4)
        {
            accumulateTimeGradBlocks(s4_solvers_, gdT);
        }
    }

    void getEnergyAndGrad(double &energy,
                          Eigen::MatrixXd &gradByCoeffs,
                          Eigen::Matrix<double, Eigen::Dynamic, Dim> &gradByPoints,
                          Eigen::VectorXd &gradByTimes,
                          Eigen::VectorXd &partialGradByTimes)
    {
        energy = 0.0;
        gradByCoeffs.setZero(2 * sys_order_ * piece_num_, Dim);
        gradByPoints.setZero(std::max(0, piece_num_ - 1), Dim);
        gradByTimes = Eigen::VectorXd::Zero(piece_num_);
        partialGradByTimes = Eigen::VectorXd::Zero(piece_num_);

        if (sys_order_ == 2)
        {
            accumulateEnergyGradBlocks(s2_solvers_, energy, gradByCoeffs,
                                       gradByPoints, gradByTimes,
                                       partialGradByTimes);
        }
        else if (sys_order_ == 3)
        {
            accumulateEnergyGradBlocks(s3_solvers_, energy, gradByCoeffs,
                                       gradByPoints, gradByTimes,
                                       partialGradByTimes);
        }
        else if (sys_order_ == 4)
        {
            accumulateEnergyGradBlocks(s4_solvers_, energy, gradByCoeffs,
                                       gradByPoints, gradByTimes,
                                       partialGradByTimes);
        }
    }

private:
    int sys_order_ = 0;
    int piece_num_ = 0;
    int instance_num_ = 0;
    Eigen::VectorXd durations_;
    std::vector<minco::MINCO_S2NU> s2_solvers_;
    std::vector<minco::MINCO_S3NU> s3_solvers_;
    std::vector<minco::MINCO_S4NU> s4_solvers_;
    std::vector<Trajectory<3>> s2_trajs_;
    std::vector<Trajectory<5>> s3_trajs_;
    std::vector<Trajectory<7>> s4_trajs_;

    template <int S>
    void fillBoundary(const Eigen::Matrix<double, Dim, Eigen::Dynamic> &head_state,
                      const Eigen::Matrix<double, Dim, Eigen::Dynamic> &tail_state,
                      const int instance,
                      Eigen::Matrix<double, 3, S> &head,
                      Eigen::Matrix<double, 3, S> &tail) const
    {
        const int start_dim = 3 * instance;
        const int active_dim = std::min(3, Dim - start_dim);
        head.topRows(active_dim) = head_state.middleRows(start_dim, active_dim);
        tail.topRows(active_dim) = tail_state.middleRows(start_dim, active_dim);
    }

    Eigen::Matrix3Xd makeInnerPoints(
        const Eigen::Matrix<double, Eigen::Dynamic, Dim> &inner_points,
        const int instance) const
    {
        Eigen::Matrix3Xd out = Eigen::Matrix3Xd::Zero(3, inner_points.rows());
        const int start_dim = 3 * instance;
        const int active_dim = std::min(3, Dim - start_dim);
        for (Eigen::Index i = 0; i < inner_points.rows(); ++i)
        {
            out.block(0, i, active_dim, 1) =
                inner_points.row(i).segment(start_dim, active_dim).transpose();
        }
        return out;
    }

    template <typename SolverVec>
    void copyCoeffBlocks(const SolverVec &solvers,
                         Eigen::MatrixXd &coeffs) const
    {
        for (int instance = 0; instance < instance_num_; ++instance)
        {
            const int start_dim = 3 * instance;
            const int active_dim = std::min(3, Dim - start_dim);
            coeffs.block(0, start_dim, coeffs.rows(), active_dim) =
                solvers[instance].getCoeffs().leftCols(active_dim);
        }
    }

    template <typename SolverVec>
    void copyCoeffGradBlocks(const SolverVec &solvers,
                             Eigen::MatrixXd &gdC) const
    {
        Eigen::MatrixX3d block_grad;
        for (int instance = 0; instance < instance_num_; ++instance)
        {
            const int start_dim = 3 * instance;
            const int active_dim = std::min(3, Dim - start_dim);
            solvers[instance].getEnergyPartialGradByCoeffs(block_grad);
            gdC.block(0, start_dim, gdC.rows(), active_dim) =
                block_grad.leftCols(active_dim);
        }
    }

    template <typename SolverVec>
    void accumulateTimeGradBlocks(const SolverVec &solvers,
                                  Eigen::VectorXd &gdT) const
    {
        Eigen::VectorXd block_grad;
        for (const auto &solver : solvers)
        {
            solver.getEnergyPartialGradByTimes(block_grad);
            gdT += block_grad;
        }
    }

    template <typename SolverVec>
    void accumulateEnergyGradBlocks(
        SolverVec &solvers,
        double &energy,
        Eigen::MatrixXd &gradByCoeffs,
        Eigen::Matrix<double, Eigen::Dynamic, Dim> &gradByPoints,
        Eigen::VectorXd &gradByTimes,
        Eigen::VectorXd &partialGradByTimes)
    {
        for (int instance = 0; instance < instance_num_; ++instance)
        {
            const int start_dim = 3 * instance;
            const int active_dim = std::min(3, Dim - start_dim);

            double block_energy = 0.0;
            Eigen::MatrixX3d block_coeff_grad;
            Eigen::VectorXd block_time_partial;
            Eigen::Matrix3Xd block_point_grad;
            Eigen::VectorXd block_time_grad;

            auto &solver = solvers[instance];
            solver.getEnergy(block_energy);
            solver.getEnergyPartialGradByCoeffs(block_coeff_grad);
            solver.getEnergyPartialGradByTimes(block_time_partial);
            solver.propogateGrad(block_coeff_grad, block_time_partial,
                                 block_point_grad, block_time_grad);

            energy += block_energy;
            gradByCoeffs.block(0, start_dim, gradByCoeffs.rows(), active_dim) =
                block_coeff_grad.leftCols(active_dim);
            partialGradByTimes += block_time_partial;
            gradByTimes += block_time_grad;

            for (int i = 0; i < piece_num_ - 1; ++i)
            {
                for (int d = 0; d < active_dim; ++d)
                {
                    gradByPoints(i, start_dim + d) =
                        block_point_grad(d, i);
                }
            }
        }
    }

    Eigen::Vector3d evaluate3(const int instance,
                              const double t,
                              const int derivative) const
    {
        if (sys_order_ == 2)
        {
            return evaluateTrajectory(s2_trajs_[instance], t, derivative);
        }
        if (sys_order_ == 3)
        {
            return evaluateTrajectory(s3_trajs_[instance], t, derivative);
        }
        return evaluateTrajectory(s4_trajs_[instance], t, derivative);
    }

    template <int Degree>
    static Eigen::Vector3d evaluateTrajectory(const Trajectory<Degree> &trajectory,
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
        case 3:
            return trajectory.getJer(t);
        case 4:
            return trajectory.getSnap(t);
        default:
            return Eigen::Vector3d::Zero();
        }
    }
};

} // namespace nubs_test

#endif
