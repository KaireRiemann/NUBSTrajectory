#ifndef NUBS_TRAJECTORY_HPP
#define NUBS_TRAJECTORY_HPP

#include <Eigen/Dense>
#include "nubs/banded_system.hpp"
#include "nubs/basis.hpp"
#include "nubs/dual.hpp"
#include "nubs/timing_stencil.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace nubs
{

template <int Dim, int MaxP = 7>
class NUBSTrajectory
{
protected:
    int s;
    int p;
    int order;
    int N_c;

    Eigen::VectorXd durations_;
    Eigen::VectorXd knots;
    mutable Eigen::MatrixXd knotJacobian;
    Eigen::Matrix<double, Eigen::Dynamic, Dim> control_points;
    double last_linear_solve_relative_residual_ = 0.0;

    static constexpr double min_duration = 1.0e-8;
    static constexpr double finite_diff_rel_eps = 1.0e-5;

    static inline double finiteDiffStepForDuration(const double Ti)
    {
        const double abs_t = std::abs(Ti);
        double h = finite_diff_rel_eps * std::max(1.0, abs_t);
        if (Ti > 0.0)
        {
            h = std::min(h, 0.25 * Ti);
        }
        if (!(h > 0.0) || !std::isfinite(h))
        {
            h = finite_diff_rel_eps;
        }
        return h;
    }

    static inline void finiteDiffDeltasForDuration(const double Ti,
                                                   double &plus_delta,
                                                   double &minus_delta)
    {
        const double h = finiteDiffStepForDuration(Ti);
        plus_delta = h;

        const double allowed_minus = std::max(0.0, Ti - min_duration);
        const double h_minus = std::min(h, allowed_minus);
        minus_delta = h_minus > 0.0 ? -h_minus : 0.0;
    }

    inline void validateGenerateInputs(const Eigen::MatrixXd &P_inner,
                                       const Eigen::MatrixXd &headState,
                                       const Eigen::MatrixXd &tailState,
                                       const Eigen::VectorXd &T) const
    {
        const int M = static_cast<int>(T.size());
        if (M <= 0)
        {
            throw std::runtime_error("NUBSTrajectory::generate(): T.size() must be positive.");
        }
        for (int i = 0; i < M; ++i)
        {
            if (!std::isfinite(T(i)) || T(i) <= min_duration)
            {
                throw std::runtime_error("NUBSTrajectory::generate(): all durations must be finite and greater than 1e-8.");
            }
        }
        if (P_inner.rows() != M - 1)
        {
            throw std::runtime_error("NUBSTrajectory::generate(): P_inner.rows() must equal T.size() - 1.");
        }
        if (P_inner.cols() != Dim)
        {
            throw std::runtime_error("NUBSTrajectory::generate(): P_inner.cols() must equal Dim.");
        }
        if (headState.rows() != Dim)
        {
            throw std::runtime_error("NUBSTrajectory::generate(): headState.rows() must equal Dim.");
        }
        if (tailState.rows() != Dim)
        {
            throw std::runtime_error("NUBSTrajectory::generate(): tailState.rows() must equal Dim.");
        }
        if (headState.cols() < s)
        {
            throw std::runtime_error("NUBSTrajectory::generate(): headState.cols() must be at least system order.");
        }
        if (tailState.cols() < s)
        {
            throw std::runtime_error("NUBSTrajectory::generate(): tailState.cols() must be at least system order.");
        }
    }

    struct Jet
    {
        double val = 0.0;
        Eigen::VectorXd dT;

        Jet() = default;
        explicit Jet(const int M) : val(0.0), dT(Eigen::VectorXd::Zero(M)) {}
    };

    struct BasisJetCache
    {
        const Eigen::VectorXd &u;
        const Eigen::MatrixXd &du_dT;
        const int numCtrl;
        const int Mdur;
        const double t;
        mutable std::map<std::tuple<int, int, int>, Jet> memo;

        BasisJetCache(const Eigen::VectorXd &u_,
                      const Eigen::MatrixXd &du_dT_,
                      const int numCtrl_,
                      const int Mdur_,
                      const double t_)
            : u(u_), du_dT(du_dT_), numCtrl(numCtrl_), Mdur(Mdur_), t(t_) {}

        inline Jet zeroJet() const
        {
            return Jet(Mdur);
        }

        inline Eigen::VectorXd knotSens(const int k) const
        {
            return du_dT.row(k).transpose();
        }

        Jet eval(const int i, const int q, const int d) const
        {
            if (d > q || q < 0 || i < 0 || i + q + 1 >= u.size())
            {
                return zeroJet();
            }

            const auto key = std::make_tuple(i, q, d);
            const auto iter = memo.find(key);
            if (iter != memo.end())
            {
                return iter->second;
            }

            const double eps = 1.0e-12;
            Jet out(Mdur);
            if (q == 0)
            {
                if (d == 0)
                {
                    out.val = ((u(i) <= t && t < u(i + 1)) ||
                               (std::abs(t - u(numCtrl)) <= eps && i == numCtrl - 1))
                                  ? 1.0
                                  : 0.0;
                }
                memo[key] = out;
                return out;
            }

            if (d == 0)
            {
                const Jet N1 = eval(i, q - 1, 0);
                const Jet N2 = eval(i + 1, q - 1, 0);
                const double den1 = u(i + q) - u(i);
                const double den2 = u(i + q + 1) - u(i + 1);

                double alpha = 0.0;
                double beta = 0.0;
                Eigen::VectorXd dAlpha = Eigen::VectorXd::Zero(Mdur);
                Eigen::VectorXd dBeta = Eigen::VectorXd::Zero(Mdur);

                if (std::abs(den1) > eps)
                {
                    const double n1 = t - u(i);
                    const Eigen::VectorXd dn1 = -knotSens(i);
                    const Eigen::VectorXd dden1 = knotSens(i + q) - knotSens(i);
                    alpha = n1 / den1;
                    dAlpha = (dn1 * den1 - n1 * dden1) / (den1 * den1);
                }
                if (std::abs(den2) > eps)
                {
                    const double n2 = u(i + q + 1) - t;
                    const Eigen::VectorXd dn2 = knotSens(i + q + 1);
                    const Eigen::VectorXd dden2 = knotSens(i + q + 1) - knotSens(i + 1);
                    beta = n2 / den2;
                    dBeta = (dn2 * den2 - n2 * dden2) / (den2 * den2);
                }

                out.val = alpha * N1.val + beta * N2.val;
                out.dT = dAlpha * N1.val + alpha * N1.dT +
                         dBeta * N2.val + beta * N2.dT;
                memo[key] = out;
                return out;
            }

            const Jet D1 = eval(i, q - 1, d - 1);
            const Jet D2 = eval(i + 1, q - 1, d - 1);
            const double den1 = u(i + q) - u(i);
            const double den2 = u(i + q + 1) - u(i + 1);

            double c1 = 0.0;
            double c2 = 0.0;
            Eigen::VectorXd dc1 = Eigen::VectorXd::Zero(Mdur);
            Eigen::VectorXd dc2 = Eigen::VectorXd::Zero(Mdur);
            if (std::abs(den1) > eps)
            {
                const Eigen::VectorXd dden1 = knotSens(i + q) - knotSens(i);
                c1 = static_cast<double>(q) / den1;
                dc1 = -static_cast<double>(q) * dden1 / (den1 * den1);
            }
            if (std::abs(den2) > eps)
            {
                const Eigen::VectorXd dden2 = knotSens(i + q + 1) - knotSens(i + 1);
                c2 = static_cast<double>(q) / den2;
                dc2 = -static_cast<double>(q) * dden2 / (den2 * den2);
            }

            out.val = c1 * D1.val - c2 * D2.val;
            out.dT = dc1 * D1.val + c1 * D1.dT -
                     dc2 * D2.val - c2 * D2.dT;
            memo[key] = out;
            return out;
        }
    };

    static inline void gaussRule(const int n,
                                 std::vector<double> &nodes,
                                 std::vector<double> &weights)
    {
        switch (n)
        {
        case 1:
            nodes = {0.0};
            weights = {2.0};
            break;
        case 2:
            nodes = {-0.5773502691896257, 0.5773502691896257};
            weights = {1.0, 1.0};
            break;
        case 3:
            nodes = {-0.7745966692414834, 0.0, 0.7745966692414834};
            weights = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
            break;
        case 4:
            nodes = {-0.8611363115940526, -0.3399810435848563,
                     0.3399810435848563, 0.8611363115940526};
            weights = {0.3478548451374539, 0.6521451548625461,
                       0.6521451548625461, 0.3478548451374539};
            break;
        default:
            nodes = {-0.9061798459386640, -0.5384693101056831, 0.0,
                     0.5384693101056831, 0.9061798459386640};
            weights = {0.2369268850561891, 0.4786286704993665,
                       0.5688888888888889, 0.4786286704993665,
                       0.2369268850561891};
            break;
        }
    }

    inline void buildSystemMatrixA(const int M,
                                   const Eigen::VectorXd &u_vec,
                                   BandedSystem &A_out) const
    {
        A_out.create(N_c, p, p);
        int row = 0;
        Eigen::Matrix<double, MaxP + 1, MaxP + 1> ders;

        for (int d = 0; d < s; ++d)
        {
            dersBasisFuns(d, p, u_vec(p), u_vec, ders);
            for (int j = 0; j <= p; ++j)
            {
                A_out(row, j) = ders(d, j);
            }
            ++row;
        }

        for (int i = 1; i < M; ++i)
        {
            const int span = p + i;
            dersBasisFuns(0, span, u_vec(span), u_vec, ders);
            for (int j = 0; j <= p; ++j)
            {
                A_out(row, span - p + j) = ders(0, j);
            }
            ++row;
        }

        const double t_end = u_vec(N_c);
        for (int d = s - 1; d >= 0; --d)
        {
            dersBasisFuns(d, N_c - 1, t_end, u_vec, ders);
            for (int j = 0; j <= p; ++j)
            {
                A_out(row, N_c - 1 - p + j) = ders(d, j);
            }
            ++row;
        }
    }

    inline Eigen::MatrixXd buildKnotJacobian(const Eigen::VectorXd &T,
                                             const int nc) const
    {
        const int M = static_cast<int>(T.size());
        const int num_knots = nc + p + 1;
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(num_knots, M);
        for (int k = 0; k < num_knots; ++k)
        {
            if (k <= p)
            {
                continue;
            }
            if (k <= nc)
            {
                const int last_seg = k - p - 1;
                for (int i = 0; i <= last_seg && i < M; ++i)
                {
                    J(k, i) = 1.0;
                }
            }
            else
            {
                J.row(k).setOnes();
            }
        }
        return J;
    }

    inline void ensureKnotJacobian() const
    {
        const int M = static_cast<int>(durations_.size());
        const int expected_rows = N_c + p + 1;
        if (knotJacobian.rows() != expected_rows ||
            knotJacobian.cols() != M)
        {
            knotJacobian = buildKnotJacobian(durations_, N_c);
        }
    }

    inline void shiftKnotsForDuration(const int duration_idx,
                                      const double delta,
                                      Eigen::VectorXd &u_vec) const
    {
        if (delta == 0.0)
        {
            return;
        }

        const int first_knot = p + 1 + duration_idx;
        for (int k = first_knot; k < u_vec.size(); ++k)
        {
            u_vec(k) += delta;
        }
    }

    inline void finiteDiffDeltasForFixedRatioTotalDuration(double &plus_delta,
                                                           double &minus_delta) const
    {
        const double total_duration = getTotalDuration();
        if (!std::isfinite(total_duration) || total_duration <= min_duration)
        {
            throw std::runtime_error(
                "NUBSTrajectory::getEnergyAndFixedRatioGrad(): total duration must be finite and positive.");
        }

        const double h = finiteDiffStepForDuration(total_duration);
        plus_delta = h;

        double min_allowed_total = min_duration;
        for (int i = 0; i < durations_.size(); ++i)
        {
            const double ratio = durations_(i) / total_duration;
            min_allowed_total =
                std::max(min_allowed_total, min_duration / ratio);
        }

        const double allowed_minus =
            std::max(0.0, total_duration - min_allowed_total);
        const double h_minus = std::min(h, allowed_minus);
        minus_delta = h_minus > 0.0 ? -h_minus : 0.0;
    }

    inline void scaleKnotsForTotalDurationDelta(const double delta,
                                                Eigen::VectorXd &u_vec) const
    {
        const double total_duration = getTotalDuration();
        const double new_total_duration = total_duration + delta;
        if (!std::isfinite(new_total_duration) ||
            new_total_duration <= min_duration)
        {
            throw std::runtime_error(
                "NUBSTrajectory::getEnergyAndFixedRatioGrad(): perturbed total duration is invalid.");
        }
        u_vec = knots * (new_total_duration / total_duration);
    }

    inline bool isKnotShiftedByDuration(const int knot_idx,
                                        const int duration_idx) const
    {
        return knot_idx >= p + 1 + duration_idx;
    }

    inline bool hasUniformShiftFlags(const int duration_idx,
                                     const std::vector<int> &knot_indices) const
    {
        bool initialized = false;
        bool reference = false;
        for (const int idx : knot_indices)
        {
            if (idx < 0 || idx >= knots.size())
            {
                continue;
            }
            const bool shifted = isKnotShiftedByDuration(idx, duration_idx);
            if (!initialized)
            {
                reference = shifted;
                initialized = true;
            }
            else if (reference != shifted)
            {
                return false;
            }
        }
        return true;
    }

    inline bool isEnergySpanAffectedByDuration(const int span,
                                               const int duration_idx) const
    {
        bool initialized = false;
        bool reference = false;
        bool uniform = true;
        auto consume = [&](const int idx)
        {
            if (idx < 0 || idx >= knots.size())
            {
                return;
            }
            const bool shifted = isKnotShiftedByDuration(idx, duration_idx);
            if (!initialized)
            {
                reference = shifted;
                initialized = true;
            }
            else if (reference != shifted)
            {
                uniform = false;
            }
        };

        consume(span);
        consume(span + 1);
        for (int k = span - p + 1; k <= span + p; ++k)
        {
            consume(k);
        }
        return !uniform;
    }

    inline std::pair<int, int> affectedEnergySpanRangeByDuration(
        const int duration_idx) const
    {
        const int valid_begin = p;
        const int valid_end = static_cast<int>(knots.size()) - p - 2;
        int affected_begin = valid_end + 1;
        int affected_end = valid_begin - 1;
        for (int span = valid_begin; span <= valid_end; ++span)
        {
            if (isEnergySpanAffectedByDuration(span, duration_idx))
            {
                affected_begin = std::min(affected_begin, span);
                affected_end = std::max(affected_end, span);
            }
        }
        return {affected_begin, affected_end};
    }

    struct ConstraintRowInfo
    {
        int derivative = 0;
        int span = 0;
        int first_col = 0;
        int eval_knot = 0;
    };

    inline ConstraintRowInfo constraintRowInfo(const int row) const
    {
        const int M = static_cast<int>(durations_.size());
        ConstraintRowInfo info;
        if (row < s)
        {
            info.derivative = row;
            info.span = p;
            info.first_col = 0;
            info.eval_knot = p;
            return info;
        }

        const int tail_start = s + M - 1;
        if (row < tail_start)
        {
            const int waypoint_idx = row - s + 1;
            info.derivative = 0;
            info.span = p + waypoint_idx;
            info.first_col = waypoint_idx;
            info.eval_knot = info.span;
            return info;
        }

        const int local_tail_idx = row - tail_start;
        info.derivative = s - 1 - local_tail_idx;
        info.span = N_c - 1;
        info.first_col = N_c - 1 - p;
        info.eval_knot = N_c;
        return info;
    }

    inline bool isConstraintRowAffectedByDuration(const int row,
                                                  const int duration_idx) const
    {
        const ConstraintRowInfo info = constraintRowInfo(row);
        bool initialized = false;
        bool reference = false;
        bool uniform = true;
        auto consume = [&](const int idx)
        {
            if (idx < 0 || idx >= knots.size())
            {
                return;
            }
            const bool shifted = isKnotShiftedByDuration(idx, duration_idx);
            if (!initialized)
            {
                reference = shifted;
                initialized = true;
            }
            else if (reference != shifted)
            {
                uniform = false;
            }
        };

        consume(info.eval_knot);
        for (int k = info.span - p + 1; k <= info.span + p; ++k)
        {
            consume(k);
        }
        return !uniform;
    }

    inline double dotSystemRowWithControlAndAdj(
        const int row,
        const Eigen::VectorXd &u_vec,
        const Eigen::Matrix<double, Eigen::Dynamic, Dim> &ctrl,
        const Eigen::MatrixXd &adjGrad) const
    {
        const ConstraintRowInfo info = constraintRowInfo(row);
        Eigen::Matrix<double, MaxP + 1, MaxP + 1> ders;
        dersBasisFuns(info.derivative, info.span, u_vec(info.eval_knot),
                      u_vec, ders);

        double result = 0.0;
        for (int j = 0; j <= p; ++j)
        {
            result += ders(info.derivative, j) *
                      adjGrad.row(row).dot(ctrl.row(info.first_col + j));
        }
        return result;
    }

    // A knot accessor for a single timing derivative.  Creating one of these
    // is O(1): no dense knot-by-duration Jacobian is materialised.
    struct LocalTimingKnotView
    {
        const Eigen::VectorXd &base_knots;
        int first_shifted_knot = 0;

        inline ad::Dual operator()(const int knot_index) const
        {
            return ad::Dual(base_knots(knot_index),
                            knot_index >= first_shifted_knot ? 1.0 : 0.0);
        }
    };

    inline double directEnergyDerivativeLocalAD(const int duration_index) const
    {
        const int M = static_cast<int>(durations_.size());
        const auto physical_spans =
            timing::affectedPhysicalSpans(p, duration_index, M);
        if (physical_spans.empty())
        {
            return 0.0;
        }

        const LocalTimingKnotView local_knots{
            knots, p + 1 + duration_index};
        std::vector<double> nodes;
        std::vector<double> weights;
        gaussRule(std::min(s, 5), nodes, weights);

        ad::Dual energy{};
        std::array<std::array<ad::Dual, MaxP + 1>, MaxP + 1> ders{};
        for (int physical_span = physical_spans.first;
             physical_span <= physical_spans.last;
             ++physical_span)
        {
            const int span = p + physical_span;
            const ad::Dual t_start = local_knots(span);
            const ad::Dual t_end = local_knots(span + 1);
            const ad::Dual length = t_end - t_start;
            if (ad::primal(length) < 1.0e-12)
            {
                continue;
            }
            const ad::Dual midpoint = (t_end + t_start) * 0.5;
            for (std::size_t q = 0; q < nodes.size(); ++q)
            {
                const ad::Dual t = midpoint + length * (0.5 * nodes[q]);
                const ad::Dual weight = length * (0.5 * weights[q]);
                basis::dersBasisFuns<ad::Dual, MaxP>(
                    p, s, span, t, local_knots, ders);

                std::array<ad::Dual, Dim> derivative{};
                for (int j = 0; j <= p; ++j)
                {
                    const int control_index = span - p + j;
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        derivative[dimension] +=
                            ders[s][j] * control_points(control_index, dimension);
                    }
                }
                ad::Dual squared_norm{};
                for (const ad::Dual &component : derivative)
                {
                    squared_norm += component * component;
                }
                energy += weight * squared_norm;
            }
        }
        return energy.derivative;
    }

    inline double systemRowAdjointDerivativeLocalAD(
        const int row,
        const int duration_index,
        const Eigen::MatrixXd &adjoint) const
    {
        const ConstraintRowInfo info = constraintRowInfo(row);
        const LocalTimingKnotView local_knots{
            knots, p + 1 + duration_index};
        std::array<std::array<ad::Dual, MaxP + 1>, MaxP + 1> ders{};
        basis::dersBasisFuns<ad::Dual, MaxP>(
            p, info.derivative, info.span, local_knots(info.eval_knot),
            local_knots, ders);

        ad::Dual projected_row{};
        for (int j = 0; j <= p; ++j)
        {
            projected_row += ders[info.derivative][j] *
                             adjoint.row(row).dot(
                                 control_points.row(info.first_col + j));
        }
        return projected_row.derivative;
    }

    inline void propagateEnergyGradLocalAD(
        const Eigen::MatrixXd &gdC,
        const Eigen::VectorXd &gdT_direct,
        Eigen::MatrixXd &gradByPoints,
        Eigen::VectorXd &gradByTimes) const
    {
        const int M = static_cast<int>(durations_.size());
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes.resize(M);

        Eigen::MatrixXd adjoint = gdC;
        A.solveAdj(adjoint);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjoint.row(s + i);
        }

        for (int duration_index = 0; duration_index < M; ++duration_index)
        {
            double constraint_term = 0.0;
            const std::vector<int> affected_rows = timing::affectedConstraintRows(
                p, s, duration_index, M);
            for (const int row : affected_rows)
            {
                constraint_term += systemRowAdjointDerivativeLocalAD(
                    row, duration_index, adjoint);
            }
            gradByTimes(duration_index) =
                gdT_direct(duration_index) - constraint_term;
        }
    }

    inline void evalLocalBasisAndTimeGrad(
        const int d,
        const int span,
        const double t,
        const Eigen::VectorXd &dt_dT,
        Eigen::Matrix<double, MaxP + 1, 1> &vals,
        Eigen::MatrixXd &dvals_dT_total) const
    {
        ensureKnotJacobian();
        const int M = durations_.size();
        vals.setZero();
        dvals_dT_total = Eigen::MatrixXd::Zero(p + 1, M);
        BasisJetCache cache(knots, knotJacobian, N_c, M, t);
        for (int j = 0; j <= p; ++j)
        {
            const int basis_idx = span - p + j;
            const Jet jd = cache.eval(basis_idx, p, d);
            const Jet jd1 = cache.eval(basis_idx, p, d + 1);
            vals(j) = jd.val;
            dvals_dT_total.row(j) = (jd.dT + jd1.val * dt_dT).transpose();
        }
    }

public:
    BandedSystem A;

    explicit NUBSTrajectory(const int sys_order = 3)
        : s(sys_order), p(2 * sys_order - 1), order(2 * sys_order), N_c(0)
    {
        if (p > MaxP)
        {
            throw std::runtime_error("Increase MaxP template argument.");
        }
    }

    inline int getS() const { return s; }
    inline int getP() const { return p; }
    inline int getCtrlPtNum(const int M) const { return M + 2 * s - 1; }
    inline int getPieceNum() const { return durations_.size(); }
    inline double getTotalDuration() const { return durations_.sum(); }
    inline const Eigen::VectorXd &getDurations() const { return durations_; }
    inline const Eigen::VectorXd &getKnots() const { return knots; }
    inline const Eigen::Matrix<double, Eigen::Dynamic, Dim> &getControlPoints() const
    {
        return control_points;
    }
    inline double getLastLinearSolveRelativeResidual() const
    {
        return last_linear_solve_relative_residual_;
    }

    inline static double conditionNumber2(const BandedSystem &system)
    {
        const int n = system.rows();
        if (n == 0)
        {
            return 1.0;
        }
        Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(n, n);
        for (int row = 0; row < n; ++row)
        {
            for (int col = 0; col < n; ++col)
            {
                dense(row, col) = system.originalCoefficient(row, col);
            }
        }
        const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            dense, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const Eigen::VectorXd singular_values = svd.singularValues();
        const double sigma_max = singular_values(0);
        const double sigma_min = singular_values(singular_values.size() - 1);
        if (sigma_min <= std::numeric_limits<double>::min())
        {
            return std::numeric_limits<double>::infinity();
        }
        return sigma_max / sigma_min;
    }

    inline double getFullSystemConditionNumber() const
    {
        return conditionNumber2(A);
    }

    inline int findSpan(const double t,
                        const int num_ctrl_pts,
                        const Eigen::VectorXd &u) const
    {
        if (t >= u(num_ctrl_pts))
        {
            return num_ctrl_pts - 1;
        }
        if (t <= u(p))
        {
            return p;
        }
        int low = p;
        int high = num_ctrl_pts;
        while (low < high)
        {
            const int mid = (low + high) / 2;
            if (t < u(mid))
            {
                high = mid;
            }
            else
            {
                low = mid + 1;
            }
        }
        return low - 1;
    }

    inline void dersBasisFuns(int n,
                              const int span,
                              const double t,
                              const Eigen::VectorXd &u,
                              Eigen::Matrix<double, MaxP + 1, MaxP + 1> &ders) const
    {
        n = std::min(n, p);
        ders.setZero();
        double ndu[MaxP + 1][MaxP + 1] = {{0.0}};
        double left[MaxP + 1] = {0.0};
        double right[MaxP + 1] = {0.0};
        ndu[0][0] = 1.0;

        for (int j = 1; j <= p; ++j)
        {
            left[j] = t - u(span + 1 - j);
            right[j] = u(span + j) - t;
            double saved = 0.0;
            for (int r = 0; r < j; ++r)
            {
                ndu[j][r] = right[r + 1] + left[j - r];
                const double temp =
                    std::abs(ndu[j][r]) < 1.0e-15 ? 0.0 : ndu[r][j - 1] / ndu[j][r];
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }

        for (int j = 0; j <= p; ++j)
        {
            ders(0, j) = ndu[j][p];
        }
        if (n == 0)
        {
            return;
        }

        double a[2][MaxP + 1] = {{0.0}};
        for (int r = 0; r <= p; ++r)
        {
            int s1 = 0;
            int s2 = 1;
            a[0][0] = 1.0;
            for (int k = 1; k <= n; ++k)
            {
                double d = 0.0;
                const int rk = r - k;
                const int pk = p - k;
                if (r >= k)
                {
                    const double den = ndu[pk + 1][rk];
                    a[s2][0] = std::abs(den) < 1.0e-15 ? 0.0 : a[s1][0] / den;
                    d = a[s2][0] * ndu[rk][pk];
                }
                const int j1 = rk >= -1 ? 1 : -rk;
                const int j2 = r - 1 <= pk ? k - 1 : p - r;
                for (int j = j1; j <= j2; ++j)
                {
                    const double den = ndu[pk + 1][rk + j];
                    a[s2][j] = std::abs(den) < 1.0e-15 ? 0.0 : (a[s1][j] - a[s1][j - 1]) / den;
                    d += a[s2][j] * ndu[rk + j][pk];
                }
                if (r <= pk)
                {
                    const double den = ndu[pk + 1][r];
                    a[s2][k] = std::abs(den) < 1.0e-15 ? 0.0 : -a[s1][k - 1] / den;
                    d += a[s2][k] * ndu[r][pk];
                }
                ders(k, r) = d;
                std::swap(s1, s2);
            }
        }

        double fac = p;
        for (int k = 1; k <= n; ++k)
        {
            for (int j = 0; j <= p; ++j)
            {
                ders(k, j) *= fac;
            }
            fac *= (p - k);
        }
    }

    inline Eigen::VectorXd generateKnots(const Eigen::VectorXd &T,
                                         const int nc) const
    {
        const int num_knots = nc + p + 1;
        Eigen::VectorXd u = Eigen::VectorXd::Zero(num_knots);
        double current_t = 0.0;
        for (int i = 0; i < T.size(); ++i)
        {
            current_t += T(i);
            u(p + 1 + i) = current_t;
        }
        for (int i = p + 1 + T.size(); i < num_knots; ++i)
        {
            u(i) = current_t;
        }
        return u;
    }

    static inline Eigen::VectorXd durationsFromRatios(
        const Eigen::VectorXd &duration_ratios,
        const double total_duration)
    {
        if (duration_ratios.size() <= 0)
        {
            throw std::runtime_error(
                "NUBSTrajectory::durationsFromRatios(): duration_ratios.size() must be positive.");
        }
        if (!std::isfinite(total_duration) || total_duration <= min_duration)
        {
            throw std::runtime_error(
                "NUBSTrajectory::durationsFromRatios(): total_duration must be finite and greater than 1e-8.");
        }

        double ratio_sum = 0.0;
        for (int i = 0; i < duration_ratios.size(); ++i)
        {
            if (!std::isfinite(duration_ratios(i)) ||
                duration_ratios(i) <= 0.0)
            {
                throw std::runtime_error(
                    "NUBSTrajectory::durationsFromRatios(): all ratios must be finite and positive.");
            }
            ratio_sum += duration_ratios(i);
        }
        if (!std::isfinite(ratio_sum) || ratio_sum <= 0.0)
        {
            throw std::runtime_error(
                "NUBSTrajectory::durationsFromRatios(): ratio sum must be finite and positive.");
        }

        Eigen::VectorXd durations =
            duration_ratios * (total_duration / ratio_sum);
        for (int i = 0; i < durations.size(); ++i)
        {
            if (!std::isfinite(durations(i)) ||
                durations(i) <= min_duration)
            {
                throw std::runtime_error(
                    "NUBSTrajectory::durationsFromRatios(): allocated duration is not greater than 1e-8.");
            }
        }
        return durations;
    }

    static inline Eigen::VectorXd uniformDurations(
        const int piece_num,
        const double total_duration)
    {
        if (piece_num <= 0)
        {
            throw std::runtime_error(
                "NUBSTrajectory::uniformDurations(): piece_num must be positive.");
        }
        return durationsFromRatios(Eigen::VectorXd::Ones(piece_num),
                                   total_duration);
    }

    inline void generate(const Eigen::MatrixXd &P_inner,
                         const Eigen::MatrixXd &headState,
                         const Eigen::MatrixXd &tailState,
                         const Eigen::VectorXd &T,
                         Eigen::MatrixXd &P_full)
    {
        validateGenerateInputs(P_inner, headState, tailState, T);
        const int M = T.size();
        N_c = getCtrlPtNum(M);
        durations_ = T;
        P_full.resize(N_c, Dim);
        knots = generateKnots(T, N_c);
        knotJacobian.resize(0, 0);
        buildSystemMatrixA(M, knots, A);

        Eigen::Matrix<double, Eigen::Dynamic, Dim> b =
            Eigen::Matrix<double, Eigen::Dynamic, Dim>::Zero(N_c, Dim);
        int row = 0;
        for (int d = 0; d < s; ++d)
        {
            b.row(row++) = headState.col(d).transpose();
        }
        for (int i = 1; i < M; ++i)
        {
            b.row(row++) = P_inner.row(i - 1);
        }
        for (int d = s - 1; d >= 0; --d)
        {
            b.row(row++) = tailState.col(d).transpose();
        }

        const Eigen::Matrix<double, Eigen::Dynamic, Dim> rhs = b;
        A.factorizeLU();
        A.solve(b);
        last_linear_solve_relative_residual_ =
            A.relativeResidualInfinityNorm(b, rhs);
        P_full = b;
        control_points = P_full;
    }

    inline void generateWithTotalDuration(const Eigen::MatrixXd &P_inner,
                                          const Eigen::MatrixXd &headState,
                                          const Eigen::MatrixXd &tailState,
                                          const Eigen::VectorXd &duration_ratios,
                                          const double total_duration,
                                          Eigen::MatrixXd &P_full)
    {
        const Eigen::VectorXd T =
            durationsFromRatios(duration_ratios, total_duration);
        generate(P_inner, headState, tailState, T, P_full);
    }

    inline void generateFixedRatio(const Eigen::MatrixXd &P_inner,
                                   const Eigen::MatrixXd &headState,
                                   const Eigen::MatrixXd &tailState,
                                   const Eigen::VectorXd &duration_ratios,
                                   const double total_duration,
                                   Eigen::MatrixXd &P_full)
    {
        generateWithTotalDuration(P_inner, headState, tailState,
                                  duration_ratios, total_duration, P_full);
    }

    inline void generateUniform(const Eigen::MatrixXd &P_inner,
                                const Eigen::MatrixXd &headState,
                                const Eigen::MatrixXd &tailState,
                                const double total_duration,
                                Eigen::MatrixXd &P_full)
    {
        const int M = static_cast<int>(P_inner.rows()) + 1;
        generate(P_inner, headState, tailState,
                 uniformDurations(M, total_duration), P_full);
    }

    inline Eigen::Matrix<double, Dim, 1> evaluate(double t,
                                                  const int d_ord = 0) const
    {
        if (d_ord > p)
        {
            return Eigen::Matrix<double, Dim, 1>::Zero();
        }
        if (t <= 0.0)
        {
            t = 0.0;
        }
        const double total_duration = getTotalDuration();
        if (t >= total_duration)
        {
            t = std::max(0.0, total_duration - 1.0e-12);
        }

        const int span = findSpan(t, N_c, knots);
        Eigen::Matrix<double, MaxP + 1, MaxP + 1> ders;
        dersBasisFuns(d_ord, span, t, knots, ders);
        Eigen::Matrix<double, Dim, 1> res =
            Eigen::Matrix<double, Dim, 1>::Zero();
        for (int j = 0; j <= p; ++j)
        {
            res += ders(d_ord, j) * control_points.row(span - p + j).transpose();
        }
        return res;
    }

    inline Eigen::Matrix<double, Dim, 1> getDerivative(
        const double t,
        const int derivative) const
    {
        return evaluate(t, derivative);
    }

    inline Eigen::Matrix<double, Dim, 1> getPos(const double t) const
    {
        return evaluate(t, 0);
    }

    inline Eigen::Matrix<double, Dim, 1> getVel(const double t) const
    {
        return evaluate(t, 1);
    }

    inline Eigen::Matrix<double, Dim, 1> getAcc(const double t) const
    {
        return evaluate(t, 2);
    }

    inline Eigen::Matrix<double, Dim, 1> getJer(const double t) const
    {
        return evaluate(t, 3);
    }

    inline Eigen::Matrix<double, Dim, 1> getJerk(const double t) const
    {
        return getJer(t);
    }

    inline Eigen::Matrix<double, Dim, 1> getSnap(const double t) const
    {
        return evaluate(t, 4);
    }

    inline void evaluatePVA(double t,
                            Eigen::Matrix<double, Dim, 1> &pos,
                            Eigen::Matrix<double, Dim, 1> &vel,
                            Eigen::Matrix<double, Dim, 1> &acc) const
    {
        pos = evaluate(t, 0);
        vel = evaluate(t, 1);
        acc = evaluate(t, 2);
    }

    inline void evaluatePVAJ(double t,
                             Eigen::Matrix<double, Dim, 1> &pos,
                             Eigen::Matrix<double, Dim, 1> &vel,
                             Eigen::Matrix<double, Dim, 1> &acc,
                             Eigen::Matrix<double, Dim, 1> &jerk) const
    {
        evaluatePVA(t, pos, vel, acc);
        jerk = evaluate(t, 3);
    }

    inline void evaluatePVAJS(double t,
                              Eigen::Matrix<double, Dim, 1> &pos,
                              Eigen::Matrix<double, Dim, 1> &vel,
                              Eigen::Matrix<double, Dim, 1> &acc,
                              Eigen::Matrix<double, Dim, 1> &jerk,
                              Eigen::Matrix<double, Dim, 1> &snap) const
    {
        evaluatePVAJ(t, pos, vel, acc, jerk);
        snap = evaluate(t, 4);
    }

    inline double getEnergyForKnots(const Eigen::VectorXd &u_vec) const
    {
        return getEnergyForKnotsRange(
            u_vec, p, static_cast<int>(u_vec.size()) - p - 2);
    }

    inline double getEnergyForKnotsRange(const Eigen::VectorXd &u_vec,
                                         const int span_begin,
                                         const int span_end) const
    {
        const int valid_begin = p;
        const int valid_end = static_cast<int>(u_vec.size()) - p - 2;
        const int begin = std::max(span_begin, valid_begin);
        const int end = std::min(span_end, valid_end);
        if (begin > end)
        {
            return 0.0;
        }

        double cost = 0.0;
        std::vector<double> nodes;
        std::vector<double> weights;
        gaussRule(std::min(s, 5), nodes, weights);
        Eigen::Matrix<double, MaxP + 1, MaxP + 1> ders;
        for (int i = begin; i <= end; ++i)
        {
            const double t_start = u_vec(i);
            const double t_end = u_vec(i + 1);
            if (t_end - t_start < 1.0e-12)
            {
                continue;
            }
            const double len = t_end - t_start;
            const double mid = 0.5 * (t_end + t_start);
            for (std::size_t k = 0; k < nodes.size(); ++k)
            {
                const double t = mid + 0.5 * len * nodes[k];
                const double w = weights[k] * 0.5 * len;
                dersBasisFuns(s, i, t, u_vec, ders);
                Eigen::Matrix<double, Dim, 1> val =
                    Eigen::Matrix<double, Dim, 1>::Zero();
                for (int j = 0; j <= p; ++j)
                {
                    val += ders(s, j) * control_points.row(i - p + j).transpose();
                }
                cost += w * val.squaredNorm();
            }
        }
        return cost;
    }

    inline double getEnergy() const
    {
        return getEnergyForKnots(knots);
    }

    inline void getEnergyPartialGradByCoeffs(double &cost,
                                             Eigen::MatrixXd &gdC) const
    {
        cost = 0.0;
        gdC.setZero(N_c, Dim);
        std::vector<double> nodes;
        std::vector<double> weights;
        gaussRule(std::min(s, 5), nodes, weights);
        Eigen::Matrix<double, MaxP + 1, MaxP + 1> ders;
        for (int i = p; i < knots.size() - p - 1; ++i)
        {
            const double t_start = knots(i);
            const double t_end = knots(i + 1);
            if (t_end - t_start < 1.0e-12)
            {
                continue;
            }
            const double len = t_end - t_start;
            const double mid = 0.5 * (t_end + t_start);
            for (std::size_t k = 0; k < nodes.size(); ++k)
            {
                const double t = mid + 0.5 * len * nodes[k];
                const double w = weights[k] * 0.5 * len;
                dersBasisFuns(s, i, t, knots, ders);
                Eigen::Matrix<double, Dim, 1> val =
                    Eigen::Matrix<double, Dim, 1>::Zero();
                for (int j = 0; j <= p; ++j)
                {
                    val += ders(s, j) * control_points.row(i - p + j).transpose();
                }
                cost += w * val.squaredNorm();
                for (int j = 0; j <= p; ++j)
                {
                    gdC.row(i - p + j) += 2.0 * w * ders(s, j) * val.transpose();
                }
            }
        }
    }

    inline void getEnergyPartialGradByTimesAnalytic(Eigen::VectorXd &gdT_direct) const
    {
        // Dense knot sensitivity path kept mainly for validation.
        ensureKnotJacobian();
        const int M = durations_.size();
        gdT_direct.setZero(M);
        std::vector<double> nodes;
        std::vector<double> weights;
        gaussRule(std::min(s, 5), nodes, weights);

        for (int span = p; span < knots.size() - p - 1; ++span)
        {
            const double t_start = knots(span);
            const double t_end = knots(span + 1);
            if (t_end - t_start < 1.0e-12)
            {
                continue;
            }
            const Eigen::VectorXd du0 = knotJacobian.row(span).transpose();
            const Eigen::VectorXd du1 = knotJacobian.row(span + 1).transpose();
            for (std::size_t k = 0; k < nodes.size(); ++k)
            {
                const double xi = nodes[k];
                const double gw = weights[k];
                const double half = 0.5 * (t_end - t_start);
                const double mid = 0.5 * (t_end + t_start);
                const double t = mid + half * xi;
                const double w = gw * half;
                const Eigen::VectorXd dt_dT =
                    0.5 * (1.0 - xi) * du0 + 0.5 * (1.0 + xi) * du1;
                const Eigen::VectorXd dw_dT = 0.5 * gw * (du1 - du0);

                Eigen::Matrix<double, MaxP + 1, 1> vals;
                Eigen::MatrixXd dvals_dT;
                evalLocalBasisAndTimeGrad(s, span, t, dt_dT, vals, dvals_dT);

                Eigen::Matrix<double, Dim, 1> acc =
                    Eigen::Matrix<double, Dim, 1>::Zero();
                Eigen::MatrixXd dacc_dT = Eigen::MatrixXd::Zero(M, Dim);
                for (int j = 0; j <= p; ++j)
                {
                    const int gidx = span - p + j;
                    acc += vals(j) * control_points.row(gidx).transpose();
                    for (int m = 0; m < M; ++m)
                    {
                        dacc_dT.row(m) += dvals_dT(j, m) * control_points.row(gidx);
                    }
                }

                const double sq = acc.squaredNorm();
                gdT_direct += dw_dT * sq;
                for (int m = 0; m < M; ++m)
                {
                    gdT_direct(m) += 2.0 * w * dacc_dT.row(m).dot(acc.transpose());
                }
            }
        }
    }

    inline void getEnergyPartialGradByTimesFiniteDiff(const Eigen::VectorXd &T,
                                                      Eigen::VectorXd &gdT_direct) const
    {
        gdT_direct.resize(T.size());

        Eigen::VectorXd u_plus = knots;
        Eigen::VectorXd u_minus = knots;
        for (int i = 0; i < T.size(); ++i)
        {
            u_plus = knots;
            u_minus = knots;

            double plus_delta = 0.0;
            double minus_delta = 0.0;
            finiteDiffDeltasForDuration(T(i), plus_delta, minus_delta);
            shiftKnotsForDuration(i, plus_delta, u_plus);
            shiftKnotsForDuration(i, minus_delta, u_minus);

            const auto span_range = affectedEnergySpanRangeByDuration(i);
            const double cost_p =
                getEnergyForKnotsRange(u_plus, span_range.first, span_range.second);
            const double cost_m =
                getEnergyForKnotsRange(u_minus, span_range.first, span_range.second);

            gdT_direct(i) = (cost_p - cost_m) / (plus_delta - minus_delta);
        }
    }

    inline void propagateGradAnalytic(const Eigen::MatrixXd &gdC,
                                      const Eigen::VectorXd &gdT_direct,
                                      Eigen::MatrixXd &gradByPoints,
                                      Eigen::VectorXd &gradByTimes) const
    {
        ensureKnotJacobian();
        const int M = durations_.size();
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes = gdT_direct;

        Eigen::MatrixXd adjGrad = gdC;
        A.solveAdj(adjGrad);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjGrad.row(s + i);
        }

        auto accumulate_row_adj = [&](const int row,
                                      const int d_ord,
                                      const int span,
                                      const double t,
                                      const Eigen::VectorXd &dt_dT)
        {
            Eigen::Matrix<double, MaxP + 1, 1> vals_dummy;
            Eigen::MatrixXd dvals_dT;
            evalLocalBasisAndTimeGrad(d_ord, span, t, dt_dT, vals_dummy, dvals_dT);
            Eigen::MatrixXd rowVecByT = Eigen::MatrixXd::Zero(M, Dim);
            for (int j = 0; j <= p; ++j)
            {
                const int gidx = span - p + j;
                for (int m = 0; m < M; ++m)
                {
                    rowVecByT.row(m) += dvals_dT(j, m) * control_points.row(gidx);
                }
            }
            for (int m = 0; m < M; ++m)
            {
                gradByTimes(m) -= adjGrad.row(row).dot(rowVecByT.row(m));
            }
        };

        int row = 0;
        for (int d = 0; d < s; ++d)
        {
            accumulate_row_adj(row++, d, p, knots(p), Eigen::VectorXd::Zero(M));
        }
        for (int i = 1; i < M; ++i)
        {
            accumulate_row_adj(row++, 0, p + i, knots(p + i),
                               knotJacobian.row(p + i).transpose());
        }
        for (int d = s - 1; d >= 0; --d)
        {
            accumulate_row_adj(row++, d, N_c - 1, knots(N_c),
                               knotJacobian.row(N_c).transpose());
        }
    }

    inline void propagateGradFiniteDiff(const Eigen::MatrixXd &gdC,
                                        const Eigen::VectorXd &gdT_direct,
                                        const Eigen::VectorXd &T,
                                        Eigen::MatrixXd &gradByPoints,
                                        Eigen::VectorXd &gradByTimes) const
    {
        const int M = T.size();
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes.resize(M);

        Eigen::MatrixXd adjGrad = gdC;
        A.solveAdj(adjGrad);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjGrad.row(s + i);
        }

        Eigen::VectorXd u_plus = knots;
        Eigen::VectorXd u_minus = knots;
        for (int i = 0; i < M; ++i)
        {
            u_plus = knots;
            u_minus = knots;

            double plus_delta = 0.0;
            double minus_delta = 0.0;
            finiteDiffDeltasForDuration(T(i), plus_delta, minus_delta);
            const double denom = plus_delta - minus_delta;
            shiftKnotsForDuration(i, plus_delta, u_plus);
            shiftKnotsForDuration(i, minus_delta, u_minus);

            double plus_adj = 0.0;
            double minus_adj = 0.0;
            for (int row = 0; row < N_c; ++row)
            {
                if (!isConstraintRowAffectedByDuration(row, i))
                {
                    continue;
                }
                plus_adj += dotSystemRowWithControlAndAdj(
                    row, u_plus, control_points, adjGrad);
                minus_adj += dotSystemRowWithControlAndAdj(
                    row, u_minus, control_points, adjGrad);
            }
            gradByTimes(i) = gdT_direct(i) - (plus_adj - minus_adj) / denom;
        }
    }

    inline void propagateEnergyGradFiniteDiffFull(const Eigen::MatrixXd &gdC,
                                                  const Eigen::VectorXd &T,
                                                  Eigen::MatrixXd &gradByPoints,
                                                  Eigen::VectorXd &gradByTimes) const
    {
        const int M = T.size();
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes.resize(M);

        Eigen::MatrixXd adjGrad = gdC;
        A.solveAdj(adjGrad);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjGrad.row(s + i);
        }

        Eigen::VectorXd u_plus = knots;
        Eigen::VectorXd u_minus = knots;
        BandedSystem A_plus;
        BandedSystem A_minus;
        for (int i = 0; i < M; ++i)
        {
            u_plus = knots;
            u_minus = knots;

            double plus_delta = 0.0;
            double minus_delta = 0.0;
            finiteDiffDeltasForDuration(T(i), plus_delta, minus_delta);
            const double denom = plus_delta - minus_delta;
            shiftKnotsForDuration(i, plus_delta, u_plus);
            shiftKnotsForDuration(i, minus_delta, u_minus);

            const double cost_p = getEnergyForKnots(u_plus);
            const double cost_m = getEnergyForKnots(u_minus);

            buildSystemMatrixA(M, u_plus, A_plus);
            buildSystemMatrixA(M, u_minus, A_minus);
            const double plus_adj = A_plus.dotMultiply(control_points, adjGrad);
            const double minus_adj = A_minus.dotMultiply(control_points, adjGrad);

            gradByTimes(i) =
                (cost_p - cost_m - plus_adj + minus_adj) / denom;
        }
    }

    inline void propagateEnergyGradFiniteDiff(const Eigen::MatrixXd &gdC,
                                              const Eigen::VectorXd &T,
                                              Eigen::MatrixXd &gradByPoints,
                                              Eigen::VectorXd &gradByTimes) const
    {
#ifdef NUBS_USE_FULL_FD_PROPAGATION
        propagateEnergyGradFiniteDiffFull(gdC, T, gradByPoints, gradByTimes);
        return;
#else
        const int M = T.size();
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes.resize(M);

        Eigen::MatrixXd adjGrad = gdC;
        A.solveAdj(adjGrad);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjGrad.row(s + i);
        }

        Eigen::VectorXd u_plus = knots;
        Eigen::VectorXd u_minus = knots;
        for (int i = 0; i < M; ++i)
        {
            u_plus = knots;
            u_minus = knots;

            double plus_delta = 0.0;
            double minus_delta = 0.0;
            finiteDiffDeltasForDuration(T(i), plus_delta, minus_delta);
            const double denom = plus_delta - minus_delta;
            shiftKnotsForDuration(i, plus_delta, u_plus);
            shiftKnotsForDuration(i, minus_delta, u_minus);

            const auto span_range = affectedEnergySpanRangeByDuration(i);
            const double cost_p =
                getEnergyForKnotsRange(u_plus, span_range.first, span_range.second);
            const double cost_m =
                getEnergyForKnotsRange(u_minus, span_range.first, span_range.second);

            double plus_adj = 0.0;
            double minus_adj = 0.0;
            for (int row = 0; row < N_c; ++row)
            {
                if (!isConstraintRowAffectedByDuration(row, i))
                {
                    continue;
                }
                plus_adj += dotSystemRowWithControlAndAdj(
                    row, u_plus, control_points, adjGrad);
                minus_adj += dotSystemRowWithControlAndAdj(
                    row, u_minus, control_points, adjGrad);
            }

#ifdef NUBS_VALIDATE_LOCAL_FD
            const double full_cost_p = getEnergyForKnots(u_plus);
            const double full_cost_m = getEnergyForKnots(u_minus);
            BandedSystem A_plus;
            BandedSystem A_minus;
            buildSystemMatrixA(M, u_plus, A_plus);
            buildSystemMatrixA(M, u_minus, A_minus);
            const double full_plus_adj =
                A_plus.dotMultiply(control_points, adjGrad);
            const double full_minus_adj =
                A_minus.dotMultiply(control_points, adjGrad);
            const double local_diff = cost_p - cost_m;
            const double full_diff = full_cost_p - full_cost_m;
            const double scale = std::max(1.0, std::abs(full_diff));
            assert(std::abs(local_diff - full_diff) <= 1.0e-8 * scale);
            const double adj_scale =
                std::max(1.0, std::abs(full_plus_adj - full_minus_adj));
            assert(std::abs((plus_adj - minus_adj) -
                            (full_plus_adj - full_minus_adj)) <=
                   1.0e-8 * adj_scale);
#endif

            gradByTimes(i) =
                (cost_p - cost_m - plus_adj + minus_adj) / denom;
        }
#endif
    }

    inline void propagateEnergyGradFixedRatioFiniteDiff(
        const Eigen::MatrixXd &gdC,
        Eigen::MatrixXd &gradByPoints,
        double &gradByTotalDuration) const
    {
        const int M = static_cast<int>(durations_.size());
        gradByPoints.resize(std::max(0, M - 1), Dim);

        Eigen::MatrixXd adjGrad = gdC;
        A.solveAdj(adjGrad);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjGrad.row(s + i);
        }

        double plus_delta = 0.0;
        double minus_delta = 0.0;
        finiteDiffDeltasForFixedRatioTotalDuration(plus_delta, minus_delta);
        const double denom = plus_delta - minus_delta;

        Eigen::VectorXd u_plus;
        Eigen::VectorXd u_minus;
        scaleKnotsForTotalDurationDelta(plus_delta, u_plus);
        scaleKnotsForTotalDurationDelta(minus_delta, u_minus);

        const double cost_p = getEnergyForKnots(u_plus);
        const double cost_m = getEnergyForKnots(u_minus);

        BandedSystem A_plus;
        BandedSystem A_minus;
        buildSystemMatrixA(M, u_plus, A_plus);
        buildSystemMatrixA(M, u_minus, A_minus);
        const double plus_adj =
            A_plus.dotMultiply(control_points, adjGrad);
        const double minus_adj =
            A_minus.dotMultiply(control_points, adjGrad);

        gradByTotalDuration =
            (cost_p - cost_m - plus_adj + minus_adj) / denom;
    }

    inline void getEnergyAndAnalyticGrad(double &cost,
                                         Eigen::MatrixXd &gradByPoints,
                                         Eigen::VectorXd &gradByTimes) const
    {
        // This validation path uses dense knot sensitivity through BasisJetCache.
        // Prefer getEnergyAndGrad() for optimization loops.
        Eigen::MatrixXd gdC;
        getEnergyPartialGradByCoeffs(cost, gdC);
        Eigen::VectorXd gdT_direct;
        getEnergyPartialGradByTimesAnalytic(gdT_direct);
        propagateGradAnalytic(gdC, gdT_direct, gradByPoints, gradByTimes);
    }

    inline void getEnergyAndFiniteDiffGrad(double &cost,
                                           Eigen::MatrixXd &gradByPoints,
                                           Eigen::VectorXd &gradByTimes) const
    {
        Eigen::MatrixXd gdC;
        getEnergyPartialGradByCoeffs(cost, gdC);
        propagateEnergyGradFiniteDiff(gdC, durations_, gradByPoints, gradByTimes);
    }

    // Direct energy sensitivities at fixed control points. This is useful for
    // profiling the local-AD construction separately from adjoint propagation.
    inline void getEnergyPartialGradByTimesLocalAD(Eigen::VectorXd &gdT_direct) const
    {
        const int M = static_cast<int>(durations_.size());
        gdT_direct.resize(M);
        for (int duration_index = 0; duration_index < M; ++duration_index)
        {
            gdT_direct(duration_index) =
                directEnergyDerivativeLocalAD(duration_index);
        }
    }

    inline void getEnergyPartialGradLocalAD(double &cost,
                                            Eigen::MatrixXd &gdC,
                                            Eigen::VectorXd &gdT_direct) const
    {
        getEnergyPartialGradByCoeffs(cost, gdC);
        getEnergyPartialGradByTimesLocalAD(gdT_direct);
    }

    // Default production path: global adjoint for the construction solve and
    // one-dimensional forward AD restricted to each duration's exact stencil.
    // Dense analytic and finite-difference methods above are validation paths.
    inline void getEnergyAndLocalADGrad(double &cost,
                                        Eigen::MatrixXd &gradByPoints,
                                        Eigen::VectorXd &gradByTimes) const
    {
        Eigen::MatrixXd gdC;
        Eigen::VectorXd gdT_direct;
        getEnergyPartialGradLocalAD(cost, gdC, gdT_direct);
        propagateEnergyGradLocalAD(gdC, gdT_direct,
                                   gradByPoints, gradByTimes);
    }

    inline void getEnergyAndGrad(double &cost,
                                 Eigen::MatrixXd &gradByPoints,
                                 Eigen::VectorXd &gradByTimes) const
    {
        getEnergyAndLocalADGrad(cost, gradByPoints, gradByTimes);
    }

    inline void getEnergyAndFiniteDiffGradFull(double &cost,
                                               Eigen::MatrixXd &gradByPoints,
                                               Eigen::VectorXd &gradByTimes) const
    {
        Eigen::MatrixXd gdC;
        getEnergyPartialGradByCoeffs(cost, gdC);
        propagateEnergyGradFiniteDiffFull(gdC, durations_, gradByPoints, gradByTimes);
    }

    inline void getEnergyAndFixedRatioGrad(double &cost,
                                           Eigen::MatrixXd &gradByPoints,
                                           double &gradByTotalDuration) const
    {
        Eigen::MatrixXd gdC;
        getEnergyPartialGradByCoeffs(cost, gdC);
        propagateEnergyGradFixedRatioFiniteDiff(gdC, gradByPoints,
                                                gradByTotalDuration);
    }

    inline void getEnergyAndTotalDurationGrad(double &cost,
                                              Eigen::MatrixXd &gradByPoints,
                                              double &gradByTotalDuration) const
    {
        getEnergyAndFixedRatioGrad(cost, gradByPoints, gradByTotalDuration);
    }

    inline void getEnergyAndUniformTimeGrad(double &cost,
                                            Eigen::MatrixXd &gradByPoints,
                                            double &gradByTotalDuration) const
    {
        getEnergyAndFixedRatioGrad(cost, gradByPoints, gradByTotalDuration);
    }
};

template <int S>
struct FixedGaussRule;

template <>
struct FixedGaussRule<2>
{
    static constexpr int Num = 2;
    inline static constexpr double nodes[Num] = {
        -0.5773502691896257, 0.5773502691896257};
    inline static constexpr double weights[Num] = {1.0, 1.0};
};

template <>
struct FixedGaussRule<3>
{
    static constexpr int Num = 3;
    inline static constexpr double nodes[Num] = {
        -0.7745966692414834, 0.0, 0.7745966692414834};
    inline static constexpr double weights[Num] = {
        0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
};

template <>
struct FixedGaussRule<4>
{
    static constexpr int Num = 4;
    inline static constexpr double nodes[Num] = {
        -0.8611363115940526, -0.3399810435848563,
        0.3399810435848563, 0.8611363115940526};
    inline static constexpr double weights[Num] = {
        0.3478548451374539, 0.6521451548625461,
        0.6521451548625461, 0.3478548451374539};
};

template <int Dim, int S>
class NUBSTrajectoryT : public NUBSTrajectory<Dim, 2 * S - 1>
{
    static_assert(S >= 2 && S <= 4,
                  "NUBSTrajectoryT currently supports S = 2, 3, and 4.");

private:
    using Base = NUBSTrajectory<Dim, 2 * S - 1>;
    using Gauss = FixedGaussRule<S>;
    static constexpr int P = 2 * S - 1;

    inline void setSystemRow(BandedSystem &A_out,
                             const int row,
                             const int first_col,
                             const int derivative,
                             const Eigen::Matrix<double, P + 1, P + 1> &ders) const
    {
        for (int j = 0; j <= P; ++j)
        {
            A_out(row, first_col + j) = ders(derivative, j);
        }
    }

    inline void dersBasisFunsFixed(int n,
                                   const int span,
                                   const double t,
                                   const Eigen::VectorXd &u,
                                   Eigen::Matrix<double, P + 1, P + 1> &ders) const
    {
        n = std::min(n, P);
        ders.setZero();
        double ndu[P + 1][P + 1] = {{0.0}};
        double left[P + 1] = {0.0};
        double right[P + 1] = {0.0};
        ndu[0][0] = 1.0;

        for (int j = 1; j <= P; ++j)
        {
            left[j] = t - u(span + 1 - j);
            right[j] = u(span + j) - t;
            double saved = 0.0;
            for (int r = 0; r < j; ++r)
            {
                ndu[j][r] = right[r + 1] + left[j - r];
                const double temp =
                    std::abs(ndu[j][r]) < 1.0e-15 ? 0.0 : ndu[r][j - 1] / ndu[j][r];
                ndu[r][j] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            ndu[j][j] = saved;
        }

        for (int j = 0; j <= P; ++j)
        {
            ders(0, j) = ndu[j][P];
        }
        if (n == 0)
        {
            return;
        }

        double a[2][P + 1] = {{0.0}};
        for (int r = 0; r <= P; ++r)
        {
            int s1 = 0;
            int s2 = 1;
            a[0][0] = 1.0;
            for (int k = 1; k <= n; ++k)
            {
                double d = 0.0;
                const int rk = r - k;
                const int pk = P - k;
                if (r >= k)
                {
                    const double den = ndu[pk + 1][rk];
                    a[s2][0] =
                        std::abs(den) < 1.0e-15 ? 0.0 : a[s1][0] / den;
                    d = a[s2][0] * ndu[rk][pk];
                }
                const int j1 = rk >= -1 ? 1 : -rk;
                const int j2 = r - 1 <= pk ? k - 1 : P - r;
                for (int j = j1; j <= j2; ++j)
                {
                    const double den = ndu[pk + 1][rk + j];
                    a[s2][j] =
                        std::abs(den) < 1.0e-15 ? 0.0 : (a[s1][j] - a[s1][j - 1]) / den;
                    d += a[s2][j] * ndu[rk + j][pk];
                }
                if (r <= pk)
                {
                    const double den = ndu[pk + 1][r];
                    a[s2][k] =
                        std::abs(den) < 1.0e-15 ? 0.0 : -a[s1][k - 1] / den;
                    d += a[s2][k] * ndu[r][pk];
                }
                ders(k, r) = d;
                std::swap(s1, s2);
            }
        }

        double fac = P;
        for (int k = 1; k <= n; ++k)
        {
            for (int j = 0; j <= P; ++j)
            {
                ders(k, j) *= fac;
            }
            fac *= (P - k);
        }
    }

    inline int findSpanFixed(const double t,
                             const int num_ctrl_pts,
                             const Eigen::VectorXd &u) const
    {
        if (t >= u(num_ctrl_pts))
        {
            return num_ctrl_pts - 1;
        }
        if (t <= u(P))
        {
            return P;
        }

        int low = P;
        int high = num_ctrl_pts;
        while (low < high)
        {
            const int mid = (low + high) / 2;
            if (t < u(mid))
            {
                high = mid;
            }
            else
            {
                low = mid + 1;
            }
        }
        return low - 1;
    }

    inline bool isKnotShiftedByDurationFixed(const int knot_idx,
                                             const int duration_idx) const
    {
        return knot_idx >= P + 1 + duration_idx;
    }

    inline bool hasUniformShiftFlagsFixed(const int duration_idx,
                                          const std::vector<int> &knot_indices) const
    {
        bool initialized = false;
        bool reference = false;
        for (const int idx : knot_indices)
        {
            if (idx < 0 || idx >= this->knots.size())
            {
                continue;
            }
            const bool shifted = isKnotShiftedByDurationFixed(idx, duration_idx);
            if (!initialized)
            {
                reference = shifted;
                initialized = true;
            }
            else if (reference != shifted)
            {
                return false;
            }
        }
        return true;
    }

    inline bool isEnergySpanAffectedByDurationFixed(const int span,
                                                    const int duration_idx) const
    {
        bool initialized = false;
        bool reference = false;
        bool uniform = true;
        auto consume = [&](const int idx)
        {
            if (idx < 0 || idx >= this->knots.size())
            {
                return;
            }
            const bool shifted = isKnotShiftedByDurationFixed(idx, duration_idx);
            if (!initialized)
            {
                reference = shifted;
                initialized = true;
            }
            else if (reference != shifted)
            {
                uniform = false;
            }
        };

        consume(span);
        consume(span + 1);
        for (int k = span - P + 1; k <= span + P; ++k)
        {
            consume(k);
        }
        return !uniform;
    }

    inline std::pair<int, int> affectedEnergySpanRangeByDurationFixed(
        const int duration_idx) const
    {
        const int valid_begin = P;
        const int valid_end = static_cast<int>(this->knots.size()) - P - 2;
        int affected_begin = valid_end + 1;
        int affected_end = valid_begin - 1;
        for (int span = valid_begin; span <= valid_end; ++span)
        {
            if (isEnergySpanAffectedByDurationFixed(span, duration_idx))
            {
                affected_begin = std::min(affected_begin, span);
                affected_end = std::max(affected_end, span);
            }
        }
        return {affected_begin, affected_end};
    }

    struct FixedConstraintRowInfo
    {
        int derivative = 0;
        int span = 0;
        int first_col = 0;
        int eval_knot = 0;
    };

    inline FixedConstraintRowInfo constraintRowInfoFixed(const int row) const
    {
        const int M = static_cast<int>(this->durations_.size());
        FixedConstraintRowInfo info;
        if (row < S)
        {
            info.derivative = row;
            info.span = P;
            info.first_col = 0;
            info.eval_knot = P;
            return info;
        }

        const int tail_start = S + M - 1;
        if (row < tail_start)
        {
            const int waypoint_idx = row - S + 1;
            info.derivative = 0;
            info.span = P + waypoint_idx;
            info.first_col = waypoint_idx;
            info.eval_knot = info.span;
            return info;
        }

        const int local_tail_idx = row - tail_start;
        info.derivative = S - 1 - local_tail_idx;
        info.span = this->N_c - 1;
        info.first_col = this->N_c - 1 - P;
        info.eval_knot = this->N_c;
        return info;
    }

    inline bool isConstraintRowAffectedByDurationFixed(const int row,
                                                       const int duration_idx) const
    {
        const FixedConstraintRowInfo info = constraintRowInfoFixed(row);
        bool initialized = false;
        bool reference = false;
        bool uniform = true;
        auto consume = [&](const int idx)
        {
            if (idx < 0 || idx >= this->knots.size())
            {
                return;
            }
            const bool shifted = isKnotShiftedByDurationFixed(idx, duration_idx);
            if (!initialized)
            {
                reference = shifted;
                initialized = true;
            }
            else if (reference != shifted)
            {
                uniform = false;
            }
        };

        consume(info.eval_knot);
        for (int k = info.span - P + 1; k <= info.span + P; ++k)
        {
            consume(k);
        }
        return !uniform;
    }

    inline double dotSystemRowWithControlAndAdjFixed(
        const int row,
        const Eigen::VectorXd &u_vec,
        const Eigen::Matrix<double, Eigen::Dynamic, Dim> &ctrl,
        const Eigen::MatrixXd &adjGrad) const
    {
        const FixedConstraintRowInfo info = constraintRowInfoFixed(row);
        Eigen::Matrix<double, P + 1, P + 1> ders;
        dersBasisFunsFixed(info.derivative, info.span, u_vec(info.eval_knot),
                           u_vec, ders);

        double result = 0.0;
        for (int j = 0; j <= P; ++j)
        {
            result += ders(info.derivative, j) *
                      adjGrad.row(row).dot(ctrl.row(info.first_col + j));
        }
        return result;
    }

    inline void buildSystemMatrixFixed(const int M,
                                       const Eigen::VectorXd &u_vec,
                                       BandedSystem &A_out) const
    {
        A_out.create(this->N_c, P, P);
        int row = 0;
        Eigen::Matrix<double, P + 1, P + 1> ders;

        for (int d = 0; d < S; ++d)
        {
            dersBasisFunsFixed(d, P, u_vec(P), u_vec, ders);
            setSystemRow(A_out, row++, 0, d, ders);
        }

        for (int i = 1; i < M; ++i)
        {
            const int span = P + i;
            dersBasisFunsFixed(0, span, u_vec(span), u_vec, ders);
            setSystemRow(A_out, row++, span - P, 0, ders);
        }

        const double t_end = u_vec(this->N_c);
        for (int d = S - 1; d >= 0; --d)
        {
            dersBasisFunsFixed(d, this->N_c - 1, t_end, u_vec, ders);
            setSystemRow(A_out, row++, this->N_c - 1 - P, d, ders);
        }
    }

    // One physical B-spline span depends on at most 2P-1 neighbouring
    // durations.  Carry all of those directions together so its basis values
    // and quadrature points are evaluated once, rather than once per T_k.
    static constexpr int LocalTimingWidth = 2 * P - 1;
    using LocalTimingJet = ad::LocalJet<LocalTimingWidth>;

    struct FixedLocalTimingKnotView
    {
        const Eigen::VectorXd &base_knots;
        int first_duration = 0;
        int duration_count = 0;

        inline LocalTimingJet operator()(const int knot_index) const
        {
            LocalTimingJet result(base_knots(knot_index));
            for (int lane = 0; lane < duration_count; ++lane)
            {
                const int duration_index = first_duration + lane;
                result.derivative[lane] =
                    knot_index >= P + 1 + duration_index ? 1.0 : 0.0;
            }
            return result;
        }
    };

    // Reverse-mode counterpart of FixedLocalTimingKnotView.  The local knot
    // window is cached as affine functions of its contiguous duration stencil,
    // so generic basis evaluation can use Reverse without rebuilding suffix
    // sums at every knot accessor call.
    struct FixedLocalReverseKnotView
    {
        ad::ReverseTape &tape;
        int first_knot = 0;
        int first_duration = 0;
        int duration_count = 0;
        std::array<ad::Reverse, LocalTimingWidth> duration_variables{};
        std::array<ad::Reverse, LocalTimingWidth + 1> shifted_sums{};
        std::array<ad::Reverse, 2 * P> local_knots{};

        FixedLocalReverseKnotView(ad::ReverseTape &tape_,
                                  const Eigen::VectorXd &base_knots,
                                  const int first_knot_,
                                  const int first_duration_,
                                  const int duration_count_)
            : tape(tape_),
              first_knot(first_knot_),
              first_duration(first_duration_),
              duration_count(duration_count_)
        {
            shifted_sums[0] = ad::Reverse(0.0);
            for (int lane = 0; lane < duration_count; ++lane)
            {
                duration_variables[lane] = tape.variable(0.0);
                shifted_sums[lane + 1] =
                    shifted_sums[lane] + duration_variables[lane];
            }
            for (int offset = 0; offset < 2 * P; ++offset)
            {
                const int knot_index = first_knot + offset;
                const int shifted_count = std::min(
                    duration_count,
                    std::max(0, knot_index - (P + 1 + first_duration) + 1));
                local_knots[offset] =
                    ad::Reverse(base_knots(knot_index)) +
                    shifted_sums[shifted_count];
            }
        }

        inline ad::Reverse operator()(const int knot_index) const
        {
            return local_knots[knot_index - first_knot];
        }

        inline double gradient(const int lane) const
        {
            return tape.gradient(duration_variables[lane]);
        }
    };

    inline timing::IndexRange affectedDurationRangeForKnotWindowFixed(
        const int first_knot,
        const int last_knot) const
    {
        const int M = static_cast<int>(this->durations_.size());
        return {std::max(0, first_knot - P),
                std::min(M - 1, last_knot - P - 1)};
    }

    // Fused local-jet pass for the direct energy terms.  It accumulates the
    // energy, control-point gradient, and all direct duration derivatives in
    // one traversal of the physical spans.
    inline void getEnergyPartialGradLocalADFused(
        double &cost,
        Eigen::MatrixXd &gdC,
        Eigen::VectorXd &gdT_direct) const
    {
        const int M = static_cast<int>(this->durations_.size());
        cost = 0.0;
        gdC.setZero(this->N_c, Dim);
        gdT_direct.setZero(M);

        std::array<std::array<LocalTimingJet, P + 1>, P + 1> ders{};
        for (int physical_span = 0; physical_span < M; ++physical_span)
        {
            const int first_duration = std::max(0, physical_span - P + 1);
            const int last_duration =
                std::min(M - 1, physical_span + P - 1);
            const int duration_count = last_duration - first_duration + 1;
            const FixedLocalTimingKnotView local_knots{
                this->knots, first_duration, duration_count};
            const int span = P + physical_span;
            const LocalTimingJet t_start = local_knots(span);
            const LocalTimingJet t_end = local_knots(span + 1);
            const LocalTimingJet length = t_end - t_start;
            if (ad::primal(length) < 1.0e-12)
            {
                continue;
            }

            const LocalTimingJet midpoint = (t_end + t_start) * 0.5;
            LocalTimingJet span_energy;
            for (int q = 0; q < Gauss::Num; ++q)
            {
                const LocalTimingJet t =
                    midpoint + length * (0.5 * Gauss::nodes[q]);
                const LocalTimingJet weight =
                    length * (0.5 * Gauss::weights[q]);
                basis::dersBasisFunsFixed<LocalTimingJet, P>(
                    S, span, t, local_knots, ders);

                std::array<LocalTimingJet, Dim> derivative{};
                for (int j = 0; j <= P; ++j)
                {
                    const int control_index = span - P + j;
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        derivative[dimension] +=
                            ders[S][j] *
                            this->control_points(control_index, dimension);
                    }
                }

                LocalTimingJet squared_norm;
                for (const LocalTimingJet &component : derivative)
                {
                    squared_norm += component * component;
                }
                span_energy += weight * squared_norm;

                for (int j = 0; j <= P; ++j)
                {
                    const int control_index = span - P + j;
                    const double coefficient =
                        2.0 * ad::primal(weight) *
                        ad::primal(ders[S][j]);
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        gdC(control_index, dimension) +=
                            coefficient * ad::primal(derivative[dimension]);
                    }
                }
            }

            cost += ad::primal(span_energy);
            for (int lane = 0; lane < duration_count; ++lane)
            {
                gdT_direct(first_duration + lane) +=
                    span_energy.derivative[lane];
            }
        }
    }

    static inline double factorialFixed(const int value)
    {
        double result = 1.0;
        for (int i = 2; i <= value; ++i)
        {
            result *= static_cast<double>(i);
        }
        return result;
    }

    inline double getEnergyExact() const
    {
        const int M = static_cast<int>(this->durations_.size());
        double cost = 0.0;
        Eigen::Matrix<double, P + 1, P + 1> ders;

        for (int physical_span = 0; physical_span < M; ++physical_span)
        {
            const int span = P + physical_span;
            const double t_start = this->knots(span);
            const double length = this->knots(span + 1) - t_start;
            if (length < 1.0e-12)
            {
                continue;
            }

            dersBasisFunsFixed(P, span, t_start, this->knots, ders);
            std::array<std::array<double, Dim>, S> derivatives{};
            for (int a = 0; a < S; ++a)
            {
                for (int j = 0; j <= P; ++j)
                {
                    const int control_index = span - P + j;
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        derivatives[a][dimension] +=
                            ders(S + a, j) *
                            this->control_points(control_index, dimension);
                    }
                }
            }

            std::array<double, P + 1> length_powers{};
            length_powers[0] = 1.0;
            for (int exponent = 1; exponent <= P; ++exponent)
            {
                length_powers[exponent] =
                    length_powers[exponent - 1] * length;
            }
            for (int a = 0; a < S; ++a)
            {
                for (int b = 0; b < S; ++b)
                {
                    const int exponent = a + b + 1;
                    const double gram_weight =
                        length_powers[exponent] /
                        (factorialFixed(a) * factorialFixed(b) *
                         static_cast<double>(exponent));
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        cost += gram_weight * derivatives[a][dimension] *
                                derivatives[b][dimension];
                    }
                }
            }
        }
        return cost;
    }

    // Exact span-energy kernel.  On span r, p^(S)(t_r + tau) is a polynomial
    // of degree S-1 because P = 2S-1.  With v_a = p^(S+a)(t_r),
    //
    // E_r = sum_{a,b=0}^{S-1} v_a^T v_b L^(a+b+1)
    //       / (a! b! (a+b+1)).
    //
    // Thus a single fixed-degree basis-derivative recurrence at the left
    // endpoint replaces the S Gauss-point recurrences previously used for a
    // span.  LocalTimingJet preserves the exact local time derivative while
    // the double-valued Gram form supplies the control-point gradient.
    inline void getEnergyPartialGradLocalADExact(
        double &cost,
        Eigen::MatrixXd &gdC,
        Eigen::VectorXd &gdT_direct) const
    {
        const int M = static_cast<int>(this->durations_.size());
        cost = 0.0;
        gdC.setZero(this->N_c, Dim);
        gdT_direct.setZero(M);

        std::array<std::array<LocalTimingJet, P + 1>, P + 1> ders{};
        for (int physical_span = 0; physical_span < M; ++physical_span)
        {
            const int first_duration = std::max(0, physical_span - P + 1);
            const int last_duration =
                std::min(M - 1, physical_span + P - 1);
            const int duration_count = last_duration - first_duration + 1;
            const FixedLocalTimingKnotView local_knots{
                this->knots, first_duration, duration_count};
            const int span = P + physical_span;
            const LocalTimingJet t_start = local_knots(span);
            const LocalTimingJet t_end = local_knots(span + 1);
            const LocalTimingJet length = t_end - t_start;
            if (ad::primal(length) < 1.0e-12)
            {
                continue;
            }

            basis::dersBasisFunsFixed<LocalTimingJet, P>(
                P, span, t_start, local_knots, ders);

            std::array<std::array<LocalTimingJet, Dim>, S> derivatives{};
            for (int a = 0; a < S; ++a)
            {
                for (int j = 0; j <= P; ++j)
                {
                    const int control_index = span - P + j;
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        derivatives[a][dimension] +=
                            ders[S + a][j] *
                            this->control_points(control_index, dimension);
                    }
                }
            }

            std::array<LocalTimingJet, P + 1> length_powers{};
            length_powers[0] = LocalTimingJet(1.0);
            for (int exponent = 1; exponent <= P; ++exponent)
            {
                length_powers[exponent] =
                    length_powers[exponent - 1] * length;
            }

            LocalTimingJet span_energy;
            std::array<std::array<double, Dim>, S> polynomial_gradient{};
            for (int a = 0; a < S; ++a)
            {
                for (int b = 0; b < S; ++b)
                {
                    const int exponent = a + b + 1;
                    const double gram_weight =
                        1.0 / (factorialFixed(a) * factorialFixed(b) *
                               static_cast<double>(exponent));
                    LocalTimingJet dot_product;
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        dot_product += derivatives[a][dimension] *
                                       derivatives[b][dimension];
                        polynomial_gradient[a][dimension] +=
                            2.0 * gram_weight *
                            ad::primal(length_powers[exponent]) *
                            ad::primal(derivatives[b][dimension]);
                    }
                    span_energy += gram_weight * length_powers[exponent] *
                                   dot_product;
                }
            }

            for (int j = 0; j <= P; ++j)
            {
                const int control_index = span - P + j;
                for (int a = 0; a < S; ++a)
                {
                    const double basis_derivative =
                        ad::primal(ders[S + a][j]);
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        gdC(control_index, dimension) +=
                            basis_derivative * polynomial_gradient[a][dimension];
                    }
                }
            }

            cost += ad::primal(span_energy);
            for (int lane = 0; lane < duration_count; ++lane)
            {
                gdT_direct(first_duration + lane) +=
                    span_energy.derivative[lane];
            }
        }
    }

    // Reverse-mode implementation of the exact Gram kernel.  Each span has a
    // scalar energy output but up to 2P-1 local duration inputs; reverse mode
    // therefore computes all timing sensitivities in one backward sweep.
    inline void getEnergyPartialGradLocalADReverse(
        double &cost,
        Eigen::MatrixXd &gdC,
        Eigen::VectorXd &gdT_direct) const
    {
        const int M = static_cast<int>(this->durations_.size());
        cost = 0.0;
        gdC.setZero(this->N_c, Dim);
        gdT_direct.setZero(M);

        ad::ReverseTape tape(4096);
        std::array<std::array<ad::Reverse, P + 1>, P + 1> ders{};
        for (int physical_span = 0; physical_span < M; ++physical_span)
        {
            tape.reset();
            const int first_duration = std::max(0, physical_span - P + 1);
            const int last_duration =
                std::min(M - 1, physical_span + P - 1);
            const int duration_count = last_duration - first_duration + 1;
            const int span = P + physical_span;
            const FixedLocalReverseKnotView local_knots{
                tape, this->knots, span - P + 1, first_duration, duration_count};
            const ad::Reverse t_start = local_knots(span);
            const ad::Reverse t_end = local_knots(span + 1);
            const ad::Reverse length = t_end - t_start;
            if (ad::primal(length) < 1.0e-12)
            {
                continue;
            }

            basis::dersBasisFunsFixed<ad::Reverse, P>(
                P, span, t_start, local_knots, ders);

            std::array<std::array<ad::Reverse, Dim>, S> derivatives{};
            for (int a = 0; a < S; ++a)
            {
                for (int j = 0; j <= P; ++j)
                {
                    const int control_index = span - P + j;
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        derivatives[a][dimension] +=
                            ders[S + a][j] *
                            this->control_points(control_index, dimension);
                    }
                }
            }

            std::array<ad::Reverse, P + 1> length_powers{};
            length_powers[0] = ad::Reverse(1.0);
            for (int exponent = 1; exponent <= P; ++exponent)
            {
                length_powers[exponent] =
                    length_powers[exponent - 1] * length;
            }

            ad::Reverse span_energy;
            std::array<std::array<double, Dim>, S> polynomial_gradient{};
            for (int a = 0; a < S; ++a)
            {
                for (int b = 0; b < S; ++b)
                {
                    const int exponent = a + b + 1;
                    const double gram_weight =
                        1.0 / (factorialFixed(a) * factorialFixed(b) *
                               static_cast<double>(exponent));
                    ad::Reverse dot_product;
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        dot_product += derivatives[a][dimension] *
                                       derivatives[b][dimension];
                        polynomial_gradient[a][dimension] +=
                            2.0 * gram_weight *
                            ad::primal(length_powers[exponent]) *
                            ad::primal(derivatives[b][dimension]);
                    }
                    span_energy += gram_weight * length_powers[exponent] *
                                   dot_product;
                }
            }

            for (int j = 0; j <= P; ++j)
            {
                const int control_index = span - P + j;
                for (int a = 0; a < S; ++a)
                {
                    const double basis_derivative =
                        ad::primal(ders[S + a][j]);
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        gdC(control_index, dimension) +=
                            basis_derivative * polynomial_gradient[a][dimension];
                    }
                }
            }

            cost += ad::primal(span_energy);
            tape.backward(span_energy);
            for (int lane = 0; lane < duration_count; ++lane)
            {
                gdT_direct(first_duration + lane) += local_knots.gradient(lane);
            }
        }
    }

    // Derivative-control formulation of the exact energy.  Repeatedly
    // differentiating a degree P B-spline curve yields a degree S-1 curve:
    //
    // D_i^(r) = (P-r+1) (D_(i+1)^(r-1) - D_i^(r-1))
    //           / (u_(i+P+1) - u_(i+r)).
    //
    // The degree S-1 derivative curve is integrated with its S-point Gauss
    // rule, which is exact for the degree 2S-2 squared norm.  This replaces a
    // degree-P, order-P basis derivative recurrence with small local divided
    // differences plus a low-degree basis evaluation.  The transpose of the
    // divided-difference recurrence supplies the exact control gradient.
    // Forward local-jet reference implementation.  It remains available for
    // cross-checking the scalar reverse path below, but is not used by the
    // production API.
    inline void getEnergyPartialGradDerivativeControlLocalJet(
        double &cost,
        Eigen::MatrixXd &gdC,
        Eigen::VectorXd &gdT_direct) const
    {
        const int M = static_cast<int>(this->durations_.size());
        cost = 0.0;
        gdC.setZero(this->N_c, Dim);
        gdT_direct.setZero(M);

        std::array<std::array<LocalTimingJet, S>, S> low_degree_basis{};
        for (int physical_span = 0; physical_span < M; ++physical_span)
        {
            const int first_duration = std::max(0, physical_span - P + 1);
            const int last_duration =
                std::min(M - 1, physical_span + P - 1);
            const int duration_count = last_duration - first_duration + 1;
            const int span = P + physical_span;
            const int first_control = span - P;
            const FixedLocalTimingKnotView local_knots{
                this->knots, first_duration, duration_count};
            const LocalTimingJet t_start = local_knots(span);
            const LocalTimingJet t_end = local_knots(span + 1);
            const LocalTimingJet length = t_end - t_start;
            if (ad::primal(length) < 1.0e-12)
            {
                continue;
            }

            std::array<std::array<LocalTimingJet, Dim>, P + 1>
                derivative_controls{};
            for (int local_index = 0; local_index <= P; ++local_index)
            {
                for (int dimension = 0; dimension < Dim; ++dimension)
                {
                    derivative_controls[local_index][dimension] =
                        LocalTimingJet(
                            this->control_points(first_control + local_index,
                                                 dimension));
                }
            }

            std::array<std::array<double, P + 1>, S> difference_scales{};
            for (int derivative = 1; derivative <= S; ++derivative)
            {
                const int output_count = P + 1 - derivative;
                for (int local_index = 0;
                     local_index < output_count;
                     ++local_index)
                {
                    const LocalTimingJet denominator =
                        local_knots(first_control + local_index + P + 1) -
                        local_knots(first_control + local_index + derivative);
                    const LocalTimingJet scale =
                        static_cast<double>(P - derivative + 1) / denominator;
                    difference_scales[derivative - 1][local_index] =
                        ad::primal(scale);
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        derivative_controls[local_index][dimension] =
                            scale *
                            (derivative_controls[local_index + 1][dimension] -
                             derivative_controls[local_index][dimension]);
                    }
                }
            }

            const auto derivative_knots = [&](const int knot_index)
            {
                return local_knots(knot_index + S);
            };
            const int derivative_span = span - S;
            LocalTimingJet span_energy;
            std::array<std::array<double, Dim>, P + 1> derivative_gradient{};
            for (int q = 0; q < Gauss::Num; ++q)
            {
                const LocalTimingJet t =
                    (t_end + t_start) * 0.5 +
                    length * (0.5 * Gauss::nodes[q]);
                const LocalTimingJet weight =
                    length * (0.5 * Gauss::weights[q]);
                basis::dersBasisFunsFixed<LocalTimingJet, S - 1>(
                    0, derivative_span, t, derivative_knots, low_degree_basis);

                std::array<LocalTimingJet, Dim> value{};
                for (int j = 0; j < S; ++j)
                {
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        value[dimension] += low_degree_basis[0][j] *
                                            derivative_controls[j][dimension];
                    }
                }
                LocalTimingJet squared_norm;
                for (const LocalTimingJet &component : value)
                {
                    squared_norm += component * component;
                }
                span_energy += weight * squared_norm;

                for (int j = 0; j < S; ++j)
                {
                    const double coefficient =
                        2.0 * ad::primal(weight) *
                        ad::primal(low_degree_basis[0][j]);
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        derivative_gradient[j][dimension] +=
                            coefficient * ad::primal(value[dimension]);
                    }
                }
            }

            // H^T application: propagate the gradient of D^(S) back to the
            // original P+1 control points through the divided differences.
            for (int derivative = S; derivative >= 1; --derivative)
            {
                std::array<std::array<double, Dim>, P + 1> previous_gradient{};
                const int output_count = P + 1 - derivative;
                for (int local_index = 0;
                     local_index < output_count;
                     ++local_index)
                {
                    const double scale =
                        difference_scales[derivative - 1][local_index];
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        const double gradient =
                            derivative_gradient[local_index][dimension];
                        previous_gradient[local_index][dimension] -=
                            scale * gradient;
                        previous_gradient[local_index + 1][dimension] +=
                            scale * gradient;
                    }
                }
                derivative_gradient = previous_gradient;
            }
            for (int local_index = 0; local_index <= P; ++local_index)
            {
                for (int dimension = 0; dimension < Dim; ++dimension)
                {
                    gdC(first_control + local_index, dimension) +=
                        derivative_gradient[local_index][dimension];
                }
            }

            cost += ad::primal(span_energy);
            for (int lane = 0; lane < duration_count; ++lane)
            {
                gdT_direct(first_duration + lane) +=
                    span_energy.derivative[lane];
            }
        }
    }

    // Reverse Algorithm A2.2 for a zero-order basis row.  The caller supplies
    // the adjoints of its Degree+1 basis values; this routine returns the
    // derivative with respect to the (possibly shifted) knot window and to
    // the evaluation time.  It is deliberately scalar: all time directions
    // are recovered in one reverse sweep instead of being carried through
    // every arithmetic operation by LocalTimingJet.
    template <int Degree>
    inline void accumulateZeroOrderBasisKnotGradient(
        const int span,
        const int knot_offset,
        const double evaluation_time,
        const std::array<double, Degree + 1> &final_basis_adjoint,
        double &evaluation_time_adjoint,
        std::array<double, 2 * Degree> &knot_adjoint) const
    {
        std::array<double, Degree + 1> values{};
        std::array<double, Degree + 1> left{};
        std::array<double, Degree + 1> right{};
        std::array<std::array<double, Degree + 1>, Degree + 1> inputs{};
        std::array<std::array<double, Degree + 1>, Degree + 1> denominators{};
        std::array<std::array<double, Degree + 1>, Degree + 1> temporaries{};
        values[0] = 1.0;

        for (int degree = 1; degree <= Degree; ++degree)
        {
            left[degree] = evaluation_time -
                           this->knots(knot_offset + span + 1 - degree);
            right[degree] = this->knots(knot_offset + span + degree) -
                            evaluation_time;
            double saved = 0.0;
            for (int local_index = 0; local_index < degree; ++local_index)
            {
                inputs[degree][local_index] = values[local_index];
                denominators[degree][local_index] =
                    right[local_index + 1] + left[degree - local_index];
                const double denominator = denominators[degree][local_index];
                const double temporary = std::abs(denominator) < 1.0e-15
                                             ? 0.0
                                             : values[local_index] / denominator;
                temporaries[degree][local_index] = temporary;
                values[local_index] = saved + right[local_index + 1] * temporary;
                saved = left[degree - local_index] * temporary;
            }
            values[degree] = saved;
        }

        std::array<double, Degree + 1> basis_adjoint = final_basis_adjoint;
        std::array<double, Degree + 1> left_adjoint{};
        std::array<double, Degree + 1> right_adjoint{};
        for (int degree = Degree; degree >= 1; --degree)
        {
            std::array<double, Degree + 1> previous_adjoint{};
            double saved_adjoint = basis_adjoint[degree];
            for (int local_index = degree - 1;
                 local_index >= 0;
                 --local_index)
            {
                const double output_adjoint = basis_adjoint[local_index];
                const double temporary = temporaries[degree][local_index];
                const double input = inputs[degree][local_index];
                const double denominator = denominators[degree][local_index];
                const int left_index = degree - local_index;
                const int right_index = local_index + 1;

                right_adjoint[right_index] += output_adjoint * temporary;
                left_adjoint[left_index] += saved_adjoint * temporary;
                const double temporary_adjoint =
                    output_adjoint * right[right_index] +
                    saved_adjoint * left[left_index];
                if (std::abs(denominator) >= 1.0e-15)
                {
                    previous_adjoint[local_index] +=
                        temporary_adjoint / denominator;
                    const double denominator_adjoint =
                        -temporary_adjoint * input /
                        (denominator * denominator);
                    right_adjoint[right_index] += denominator_adjoint;
                    left_adjoint[left_index] += denominator_adjoint;
                }
                saved_adjoint = output_adjoint;
            }
            basis_adjoint = previous_adjoint;
        }

        knot_adjoint.fill(0.0);
        evaluation_time_adjoint = 0.0;
        const int first_knot = knot_offset + span - Degree + 1;
        for (int local_index = 1; local_index <= Degree; ++local_index)
        {
            evaluation_time_adjoint += left_adjoint[local_index] -
                                       right_adjoint[local_index];
            knot_adjoint[knot_offset + span + 1 - local_index - first_knot] -=
                left_adjoint[local_index];
            knot_adjoint[knot_offset + span + local_index - first_knot] +=
                right_adjoint[local_index];
        }
    }

    // Exact derivative-control energy with scalar reverse propagation in the
    // duration variables.  In contrast with the local-jet reference path,
    // it differentiates the divided-difference recurrence and the low-degree
    // basis recurrence explicitly.  Thus its work is proportional to the
    // small local stencil, rather than to (stencil width x every primitive
    // arithmetic operation).
    inline void getEnergyPartialGradDerivativeControl(
        double &cost,
        Eigen::MatrixXd &gdC,
        Eigen::VectorXd &gdT_direct) const
    {
        const int M = static_cast<int>(this->durations_.size());
        cost = 0.0;
        gdC.setZero(this->N_c, Dim);
        gdT_direct.setZero(M);

        for (int physical_span = 0; physical_span < M; ++physical_span)
        {
            const int span = P + physical_span;
            const int first_control = span - P;
            const int first_duration = std::max(0, physical_span - P + 1);
            const int last_duration =
                std::min(M - 1, physical_span + P - 1);
            const timing::IndexRange duration_range{
                first_duration, last_duration};
            const double t_start = this->knots(span);
            const double t_end = this->knots(span + 1);
            const double length = t_end - t_start;
            if (length < 1.0e-12)
            {
                continue;
            }

            using LocalVector = std::array<double, Dim>;
            std::array<std::array<LocalVector, P + 1>, S + 1>
                derivative_history{};
            for (int local_index = 0; local_index <= P; ++local_index)
            {
                for (int dimension = 0; dimension < Dim; ++dimension)
                {
                    derivative_history[0][local_index][dimension] =
                        this->control_points(first_control + local_index,
                                             dimension);
                }
            }

            std::array<std::array<double, P + 1>, S> difference_scales{};
            std::array<std::array<double, P + 1>, S> difference_denominators{};
            for (int derivative = 1; derivative <= S; ++derivative)
            {
                const int output_count = P + 1 - derivative;
                for (int local_index = 0;
                     local_index < output_count;
                     ++local_index)
                {
                    const double denominator =
                        this->knots(first_control + local_index + P + 1) -
                        this->knots(first_control + local_index + derivative);
                    const double scale =
                        static_cast<double>(P - derivative + 1) / denominator;
                    difference_scales[derivative - 1][local_index] = scale;
                    difference_denominators[derivative - 1][local_index] =
                        denominator;
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        derivative_history[derivative][local_index][dimension] =
                            scale *
                            (derivative_history[derivative - 1][local_index + 1]
                                                      [dimension] -
                             derivative_history[derivative - 1][local_index]
                                                      [dimension]);
                    }
                }
            }

            const int derivative_span = span - S;
            std::array<std::array<double, S>, S> low_degree_basis{};
            std::array<LocalVector, P + 1> derivative_gradient{};
            std::array<double, 2 * P + 1> knot_adjoint{};
            const int first_knot = span - P + 1;
            const auto add_knot_adjoint = [&](const int knot_index,
                                              const double gradient)
            {
                knot_adjoint[knot_index - first_knot] += gradient;
            };

            for (int q = 0; q < Gauss::Num; ++q)
            {
                const double midpoint = 0.5 * (t_start + t_end);
                const double t = midpoint +
                                 0.5 * length * Gauss::nodes[q];
                const double weight = 0.5 * length * Gauss::weights[q];
                const auto derivative_knots = [&](const int knot_index)
                {
                    return this->knots(knot_index + S);
                };
                basis::dersBasisFunsFixed<double, S - 1>(
                    0, derivative_span, t, derivative_knots, low_degree_basis);

                LocalVector value{};
                for (int local_index = 0; local_index < S; ++local_index)
                {
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        value[dimension] +=
                            low_degree_basis[0][local_index] *
                            derivative_history[S][local_index][dimension];
                    }
                }
                double squared_norm = 0.0;
                for (const double component : value)
                {
                    squared_norm += component * component;
                }
                cost += weight * squared_norm;

                std::array<double, S> basis_adjoint{};
                for (int local_index = 0; local_index < S; ++local_index)
                {
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        const double value_adjoint =
                            2.0 * weight * value[dimension];
                        derivative_gradient[local_index][dimension] +=
                            low_degree_basis[0][local_index] * value_adjoint;
                        basis_adjoint[local_index] +=
                            derivative_history[S][local_index][dimension] *
                            value_adjoint;
                    }
                }

                std::array<double, 2 * (S - 1)> basis_knot_adjoint{};
                double t_adjoint = 0.0;
                accumulateZeroOrderBasisKnotGradient<S - 1>(
                    derivative_span, S, t, basis_adjoint, t_adjoint,
                    basis_knot_adjoint);
                const int basis_first_knot =
                    S + derivative_span - (S - 1) + 1;
                for (int local_index = 0;
                     local_index < 2 * (S - 1);
                     ++local_index)
                {
                    add_knot_adjoint(basis_first_knot + local_index,
                                     basis_knot_adjoint[local_index]);
                }

                // t = ((1-xi)/2) t_start + ((1+xi)/2) t_end, while
                // weight = (t_end-t_start) * gauss_weight / 2.
                add_knot_adjoint(
                    span, t_adjoint * 0.5 * (1.0 - Gauss::nodes[q]) -
                              0.5 * Gauss::weights[q] * squared_norm);
                add_knot_adjoint(
                    span + 1, t_adjoint * 0.5 * (1.0 + Gauss::nodes[q]) +
                                      0.5 * Gauss::weights[q] * squared_norm);
            }

            // Apply H^T and simultaneously differentiate every local
            // divided-difference scale with respect to its two knots.
            for (int derivative = S; derivative >= 1; --derivative)
            {
                std::array<LocalVector, P + 1> previous_gradient{};
                const int output_count = P + 1 - derivative;
                for (int local_index = 0;
                     local_index < output_count;
                     ++local_index)
                {
                    const double scale =
                        difference_scales[derivative - 1][local_index];
                    const double denominator =
                        difference_denominators[derivative - 1][local_index];
                    double scale_adjoint = 0.0;
                    for (int dimension = 0; dimension < Dim; ++dimension)
                    {
                        const double gradient =
                            derivative_gradient[local_index][dimension];
                        scale_adjoint += gradient *
                            (derivative_history[derivative - 1][local_index + 1]
                                                      [dimension] -
                             derivative_history[derivative - 1][local_index]
                                                      [dimension]);
                        previous_gradient[local_index][dimension] -=
                            scale * gradient;
                        previous_gradient[local_index + 1][dimension] +=
                            scale * gradient;
                    }
                    const double denominator_adjoint =
                        -scale_adjoint * scale / denominator;
                    add_knot_adjoint(
                        first_control + local_index + P + 1,
                        denominator_adjoint);
                    add_knot_adjoint(first_control + local_index + derivative,
                                     -denominator_adjoint);
                }
                derivative_gradient = previous_gradient;
            }
            for (int local_index = 0; local_index <= P; ++local_index)
            {
                for (int dimension = 0; dimension < Dim; ++dimension)
                {
                    gdC(first_control + local_index, dimension) +=
                        derivative_gradient[local_index][dimension];
                }
            }

            for (int duration = duration_range.first;
                 duration <= duration_range.last;
                 ++duration)
            {
                const int first_shifted_knot = P + 1 + duration;
                double duration_adjoint = 0.0;
                for (int local_index = 0; local_index < 2 * P + 1;
                     ++local_index)
                {
                    if (first_knot + local_index >= first_shifted_knot)
                    {
                        duration_adjoint += knot_adjoint[local_index];
                    }
                }
                gdT_direct(duration) += duration_adjoint;
            }
        }
    }

    // Exact reverse of Algorithm A2.2 for a zero-order basis row evaluated at
    // an interior waypoint knot.  Waypoint rows are the O(M) part of A(T).
    // Keeping this reverse scalar avoids carrying a LocalTimingJet through
    // every recurrence operation while still returning the whole local knot
    // stencil derivative of lambda_row^T A_row(T) C.
    inline void accumulateWaypointRowTimeGradient(
        const FixedConstraintRowInfo &info,
        const Eigen::MatrixXd &adjoint,
        const int row,
        Eigen::VectorXd &gradByTimes) const
    {
        const int span = info.span;
        const int first_knot = span - P + 1;
        const timing::IndexRange duration_range =
            affectedDurationRangeForKnotWindowFixed(first_knot, span + P);

        std::array<double, P + 1> basis_adjoint{};
        for (int local_index = 0; local_index <= P; ++local_index)
        {
            basis_adjoint[local_index] = adjoint.row(row).dot(
                this->control_points.row(info.first_col + local_index));
        }
        std::array<double, 2 * P> knot_adjoint{};
        double evaluation_time_adjoint = 0.0;
        accumulateZeroOrderBasisKnotGradient<P>(
            span, 0, this->knots(span), basis_adjoint,
            evaluation_time_adjoint, knot_adjoint);
        knot_adjoint[span - first_knot] += evaluation_time_adjoint;

        for (int duration = duration_range.first;
             duration <= duration_range.last;
             ++duration)
        {
            double projected_row_derivative = 0.0;
            const int first_shifted_knot = P + 1 + duration;
            for (int local_index = 0; local_index < 2 * P; ++local_index)
            {
                if (first_knot + local_index >= first_shifted_knot)
                {
                    projected_row_derivative += knot_adjoint[local_index];
                }
            }
            gradByTimes(duration) -= projected_row_derivative;
        }
    }

    // The construction correction has the same compact timing support as the
    // energy.  Evaluate each constraint row once and scatter its local jet
    // instead of repeating the same row for every affected duration.
    inline void propagateEnergyGradLocalADFused(
        const Eigen::MatrixXd &gdC,
        const Eigen::VectorXd &gdT_direct,
        Eigen::MatrixXd &gradByPoints,
        Eigen::VectorXd &gradByTimes) const
    {
        const int M = static_cast<int>(this->durations_.size());
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes = gdT_direct;

        Eigen::MatrixXd adjoint = gdC;
        this->A.solveAdj(adjoint);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjoint.row(S + i);
        }

        std::array<std::array<LocalTimingJet, P + 1>, P + 1> ders{};
        for (int row = 0; row < this->N_c; ++row)
        {
            const FixedConstraintRowInfo info = constraintRowInfoFixed(row);
            const int first_knot =
                std::min(info.eval_knot, info.span - P + 1);
            const int last_knot =
                std::max(info.eval_knot, info.span + P);
            const timing::IndexRange duration_range =
                affectedDurationRangeForKnotWindowFixed(first_knot, last_knot);
            if (duration_range.empty())
            {
                continue;
            }
            const int duration_count =
                duration_range.last - duration_range.first + 1;
            const FixedLocalTimingKnotView local_knots{
                this->knots, duration_range.first, duration_count};
            basis::dersBasisFunsFixed<LocalTimingJet, P>(
                info.derivative, info.span, local_knots(info.eval_knot),
                local_knots, ders);

            LocalTimingJet projected_row;
            for (int j = 0; j <= P; ++j)
            {
                projected_row += ders[info.derivative][j] *
                                 adjoint.row(row).dot(
                                     this->control_points.row(info.first_col + j));
            }
            for (int lane = 0; lane < duration_count; ++lane)
            {
                gradByTimes(duration_range.first + lane) -=
                    projected_row.derivative[lane];
            }
        }
    }

    // Reverse-mode construction correction.  This is the exact local
    // derivative of lambda^T A(T) C, evaluated once per row instead of once
    // per duration direction.
    inline void propagateEnergyGradLocalADReverse(
        const Eigen::MatrixXd &gdC,
        const Eigen::VectorXd &gdT_direct,
        Eigen::MatrixXd &gradByPoints,
        Eigen::VectorXd &gradByTimes) const
    {
        const int M = static_cast<int>(this->durations_.size());
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes = gdT_direct;

        Eigen::MatrixXd adjoint = gdC;
        this->A.solveAdj(adjoint);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjoint.row(S + i);
        }

        ad::ReverseTape tape(4096);
        std::array<std::array<ad::Reverse, P + 1>, P + 1> ders{};
        for (int row = 0; row < this->N_c; ++row)
        {
            const FixedConstraintRowInfo info = constraintRowInfoFixed(row);
            const int first_knot =
                std::min(info.eval_knot, info.span - P + 1);
            const int last_knot =
                std::max(info.eval_knot, info.span + P);
            const timing::IndexRange duration_range =
                affectedDurationRangeForKnotWindowFixed(first_knot, last_knot);
            if (duration_range.empty())
            {
                continue;
            }

            tape.reset();
            const int duration_count =
                duration_range.last - duration_range.first + 1;
            const FixedLocalReverseKnotView local_knots{
                tape, this->knots, first_knot, duration_range.first,
                duration_count};
            basis::dersBasisFunsFixed<ad::Reverse, P>(
                info.derivative, info.span, local_knots(info.eval_knot),
                local_knots, ders);

            ad::Reverse projected_row;
            for (int j = 0; j <= P; ++j)
            {
                projected_row += ders[info.derivative][j] *
                                 adjoint.row(row).dot(
                                     this->control_points.row(info.first_col + j));
            }
            tape.backward(projected_row);
            for (int lane = 0; lane < duration_count; ++lane)
            {
                gradByTimes(duration_range.first + lane) -=
                    local_knots.gradient(lane);
            }
        }
    }

    // Hybrid analytic propagation: interior interpolation rows use the
    // scalar basis reverse above; the O(S) endpoint derivative rows retain
    // the existing LocalTimingJet implementation.
    inline void propagateEnergyGradDerivativeControl(
        const Eigen::MatrixXd &gdC,
        const Eigen::VectorXd &gdT_direct,
        Eigen::MatrixXd &gradByPoints,
        Eigen::VectorXd &gradByTimes) const
    {
        const int M = static_cast<int>(this->durations_.size());
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes = gdT_direct;

        Eigen::MatrixXd adjoint = gdC;
        this->A.solveAdj(adjoint);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjoint.row(S + i);
        }

        const int tail_start = S + M - 1;
        std::array<std::array<LocalTimingJet, P + 1>, P + 1> ders{};
        for (int row = 0; row < this->N_c; ++row)
        {
            const FixedConstraintRowInfo info = constraintRowInfoFixed(row);
            if (row >= S && row < tail_start)
            {
                accumulateWaypointRowTimeGradient(info, adjoint, row,
                                                   gradByTimes);
                continue;
            }

            const int first_knot =
                std::min(info.eval_knot, info.span - P + 1);
            const int last_knot =
                std::max(info.eval_knot, info.span + P);
            const timing::IndexRange duration_range =
                affectedDurationRangeForKnotWindowFixed(first_knot, last_knot);
            if (duration_range.empty())
            {
                continue;
            }
            const int duration_count =
                duration_range.last - duration_range.first + 1;
            const FixedLocalTimingKnotView local_knots{
                this->knots, duration_range.first, duration_count};
            basis::dersBasisFunsFixed<LocalTimingJet, P>(
                info.derivative, info.span, local_knots(info.eval_knot),
                local_knots, ders);

            LocalTimingJet projected_row;
            for (int j = 0; j <= P; ++j)
            {
                projected_row += ders[info.derivative][j] *
                                 adjoint.row(row).dot(
                                     this->control_points.row(info.first_col + j));
            }
            for (int lane = 0; lane < duration_count; ++lane)
            {
                gradByTimes(duration_range.first + lane) -=
                    projected_row.derivative[lane];
            }
        }
    }

public:
    static constexpr int SystemOrder = S;
    static constexpr int Degree = P;

    NUBSTrajectoryT() : Base(S) {}

    inline void generate(const Eigen::MatrixXd &P_inner,
                         const Eigen::MatrixXd &headState,
                         const Eigen::MatrixXd &tailState,
                         const Eigen::VectorXd &T,
                         Eigen::MatrixXd &P_full)
    {
        this->validateGenerateInputs(P_inner, headState, tailState, T);
        const int M = static_cast<int>(T.size());
        this->N_c = this->getCtrlPtNum(M);
        this->durations_ = T;
        P_full.resize(this->N_c, Dim);
        this->knots = this->generateKnots(T, this->N_c);
        this->knotJacobian.resize(0, 0);
        buildSystemMatrixFixed(M, this->knots, this->A);

        Eigen::Matrix<double, Eigen::Dynamic, Dim> b =
            Eigen::Matrix<double, Eigen::Dynamic, Dim>::Zero(this->N_c, Dim);
        int row = 0;
        for (int d = 0; d < S; ++d)
        {
            b.row(row++) = headState.col(d).transpose();
        }
        for (int i = 1; i < M; ++i)
        {
            b.row(row++) = P_inner.row(i - 1);
        }
        for (int d = S - 1; d >= 0; --d)
        {
            b.row(row++) = tailState.col(d).transpose();
        }

        const Eigen::Matrix<double, Eigen::Dynamic, Dim> rhs = b;
        this->A.factorizeLU();
        this->A.solve(b);
        this->last_linear_solve_relative_residual_ =
            this->A.relativeResidualInfinityNorm(b, rhs);
        P_full = b;
        this->control_points = P_full;
    }

    inline void generateWithTotalDuration(const Eigen::MatrixXd &P_inner,
                                          const Eigen::MatrixXd &headState,
                                          const Eigen::MatrixXd &tailState,
                                          const Eigen::VectorXd &duration_ratios,
                                          const double total_duration,
                                          Eigen::MatrixXd &P_full)
    {
        const Eigen::VectorXd T =
            this->durationsFromRatios(duration_ratios, total_duration);
        generate(P_inner, headState, tailState, T, P_full);
    }

    inline void generateFixedRatio(const Eigen::MatrixXd &P_inner,
                                   const Eigen::MatrixXd &headState,
                                   const Eigen::MatrixXd &tailState,
                                   const Eigen::VectorXd &duration_ratios,
                                   const double total_duration,
                                   Eigen::MatrixXd &P_full)
    {
        generateWithTotalDuration(P_inner, headState, tailState,
                                  duration_ratios, total_duration, P_full);
    }

    inline void generateUniform(const Eigen::MatrixXd &P_inner,
                                const Eigen::MatrixXd &headState,
                                const Eigen::MatrixXd &tailState,
                                const double total_duration,
                                Eigen::MatrixXd &P_full)
    {
        const int M = static_cast<int>(P_inner.rows()) + 1;
        generate(P_inner, headState, tailState,
                 this->uniformDurations(M, total_duration), P_full);
    }

    inline Eigen::Matrix<double, Dim, 1> evaluate(double t,
                                                  const int d_ord = 0) const
    {
        if (d_ord > P)
        {
            return Eigen::Matrix<double, Dim, 1>::Zero();
        }
        if (t <= 0.0)
        {
            t = 0.0;
        }
        const double total_duration = this->getTotalDuration();
        if (t >= total_duration)
        {
            t = std::max(0.0, total_duration - 1.0e-12);
        }

        const int span = findSpanFixed(t, this->N_c, this->knots);
        Eigen::Matrix<double, P + 1, P + 1> ders;
        dersBasisFunsFixed(d_ord, span, t, this->knots, ders);
        Eigen::Matrix<double, Dim, 1> res =
            Eigen::Matrix<double, Dim, 1>::Zero();
        for (int j = 0; j <= P; ++j)
        {
            res += ders(d_ord, j) *
                   this->control_points.row(span - P + j).transpose();
        }
        return res;
    }

    inline Eigen::Matrix<double, Dim, 1> getDerivative(
        const double t,
        const int derivative) const
    {
        return evaluate(t, derivative);
    }

    inline Eigen::Matrix<double, Dim, 1> getPos(const double t) const
    {
        return evaluate(t, 0);
    }

    inline Eigen::Matrix<double, Dim, 1> getVel(const double t) const
    {
        return evaluate(t, 1);
    }

    inline Eigen::Matrix<double, Dim, 1> getAcc(const double t) const
    {
        return evaluate(t, 2);
    }

    inline Eigen::Matrix<double, Dim, 1> getJer(const double t) const
    {
        return evaluate(t, 3);
    }

    inline Eigen::Matrix<double, Dim, 1> getJerk(const double t) const
    {
        return getJer(t);
    }

    inline Eigen::Matrix<double, Dim, 1> getSnap(const double t) const
    {
        return evaluate(t, 4);
    }

    inline void evaluatePVA(double t,
                            Eigen::Matrix<double, Dim, 1> &pos,
                            Eigen::Matrix<double, Dim, 1> &vel,
                            Eigen::Matrix<double, Dim, 1> &acc) const
    {
        pos = evaluate(t, 0);
        vel = evaluate(t, 1);
        acc = evaluate(t, 2);
    }

    inline void evaluatePVAJ(double t,
                             Eigen::Matrix<double, Dim, 1> &pos,
                             Eigen::Matrix<double, Dim, 1> &vel,
                             Eigen::Matrix<double, Dim, 1> &acc,
                             Eigen::Matrix<double, Dim, 1> &jerk) const
    {
        evaluatePVA(t, pos, vel, acc);
        jerk = evaluate(t, 3);
    }

    inline void evaluatePVAJS(double t,
                              Eigen::Matrix<double, Dim, 1> &pos,
                              Eigen::Matrix<double, Dim, 1> &vel,
                              Eigen::Matrix<double, Dim, 1> &acc,
                              Eigen::Matrix<double, Dim, 1> &jerk,
                              Eigen::Matrix<double, Dim, 1> &snap) const
    {
        evaluatePVAJ(t, pos, vel, acc, jerk);
        snap = evaluate(t, 4);
    }

    inline double getEnergyForKnots(const Eigen::VectorXd &u_vec) const
    {
        return getEnergyForKnotsRange(
            u_vec, P, static_cast<int>(u_vec.size()) - P - 2);
    }

    inline double getEnergyForKnotsRange(const Eigen::VectorXd &u_vec,
                                         const int span_begin,
                                         const int span_end) const
    {
        const int valid_begin = P;
        const int valid_end = static_cast<int>(u_vec.size()) - P - 2;
        const int begin = std::max(span_begin, valid_begin);
        const int end = std::min(span_end, valid_end);
        if (begin > end)
        {
            return 0.0;
        }

        double cost = 0.0;
        Eigen::Matrix<double, P + 1, P + 1> ders;

        for (int span = begin; span <= end; ++span)
        {
            const double t_start = u_vec(span);
            const double t_end = u_vec(span + 1);
            if (t_end - t_start < 1.0e-12)
            {
                continue;
            }

            const double len = t_end - t_start;
            const double mid = 0.5 * (t_end + t_start);
            for (int k = 0; k < Gauss::Num; ++k)
            {
                const double t = mid + 0.5 * len * Gauss::nodes[k];
                const double w = Gauss::weights[k] * 0.5 * len;
                dersBasisFunsFixed(S, span, t, u_vec, ders);

                Eigen::Matrix<double, Dim, 1> val =
                    Eigen::Matrix<double, Dim, 1>::Zero();
                for (int j = 0; j <= P; ++j)
                {
                    val += ders(S, j) *
                           this->control_points.row(span - P + j).transpose();
                }
                cost += w * val.squaredNorm();
            }
        }
        return cost;
    }

    inline double getEnergy() const
    {
        return getEnergyExact();
    }

    inline void getEnergyPartialGradByCoeffs(double &cost,
                                             Eigen::MatrixXd &gdC) const
    {
        cost = 0.0;
        gdC.setZero(this->N_c, Dim);

        Eigen::Matrix<double, P + 1, P + 1> ders;

        for (int span = P; span < this->knots.size() - P - 1; ++span)
        {
            const double t_start = this->knots(span);
            const double t_end = this->knots(span + 1);
            if (t_end - t_start < 1.0e-12)
            {
                continue;
            }

            const double len = t_end - t_start;
            const double mid = 0.5 * (t_end + t_start);
            for (int k = 0; k < Gauss::Num; ++k)
            {
                const double t = mid + 0.5 * len * Gauss::nodes[k];
                const double w = Gauss::weights[k] * 0.5 * len;
                dersBasisFunsFixed(S, span, t, this->knots, ders);

                Eigen::Matrix<double, Dim, 1> val =
                    Eigen::Matrix<double, Dim, 1>::Zero();
                for (int j = 0; j <= P; ++j)
                {
                    val += ders(S, j) *
                           this->control_points.row(span - P + j).transpose();
                }

                cost += w * val.squaredNorm();
                for (int j = 0; j <= P; ++j)
                {
                    gdC.row(span - P + j) +=
                        2.0 * w * ders(S, j) * val.transpose();
                }
            }
        }
    }

    inline void getEnergyPartialGradByTimesAnalytic(Eigen::VectorXd &gdT_direct) const
    {
        // Dense knot sensitivity path kept mainly for validation.
        this->ensureKnotJacobian();
        const int M = static_cast<int>(this->durations_.size());
        gdT_direct.setZero(M);

        for (int span = P; span < this->knots.size() - P - 1; ++span)
        {
            const double t_start = this->knots(span);
            const double t_end = this->knots(span + 1);
            if (t_end - t_start < 1.0e-12)
            {
                continue;
            }

            const Eigen::VectorXd du0 = this->knotJacobian.row(span).transpose();
            const Eigen::VectorXd du1 = this->knotJacobian.row(span + 1).transpose();
            for (int k = 0; k < Gauss::Num; ++k)
            {
                const double xi = Gauss::nodes[k];
                const double gw = Gauss::weights[k];
                const double half = 0.5 * (t_end - t_start);
                const double mid = 0.5 * (t_end + t_start);
                const double t = mid + half * xi;
                const double w = gw * half;
                const Eigen::VectorXd dt_dT =
                    0.5 * (1.0 - xi) * du0 + 0.5 * (1.0 + xi) * du1;
                const Eigen::VectorXd dw_dT = 0.5 * gw * (du1 - du0);

                Eigen::Matrix<double, P + 1, 1> vals;
                Eigen::MatrixXd dvals_dT;
                this->evalLocalBasisAndTimeGrad(S, span, t, dt_dT, vals, dvals_dT);

                Eigen::Matrix<double, Dim, 1> acc =
                    Eigen::Matrix<double, Dim, 1>::Zero();
                Eigen::MatrixXd dacc_dT = Eigen::MatrixXd::Zero(M, Dim);
                for (int j = 0; j <= P; ++j)
                {
                    const int gidx = span - P + j;
                    acc += vals(j) * this->control_points.row(gidx).transpose();
                    for (int m = 0; m < M; ++m)
                    {
                        dacc_dT.row(m) +=
                            dvals_dT(j, m) * this->control_points.row(gidx);
                    }
                }

                const double sq = acc.squaredNorm();
                gdT_direct += dw_dT * sq;
                for (int m = 0; m < M; ++m)
                {
                    gdT_direct(m) +=
                        2.0 * w * dacc_dT.row(m).dot(acc.transpose());
                }
            }
        }
    }

    inline void getEnergyPartialGradByTimesFiniteDiff(const Eigen::VectorXd &T,
                                                      Eigen::VectorXd &gdT_direct) const
    {
        gdT_direct.resize(T.size());

        Eigen::VectorXd u_plus = this->knots;
        Eigen::VectorXd u_minus = this->knots;
        for (int i = 0; i < T.size(); ++i)
        {
            u_plus = this->knots;
            u_minus = this->knots;

            double plus_delta = 0.0;
            double minus_delta = 0.0;
            this->finiteDiffDeltasForDuration(T(i), plus_delta, minus_delta);
            this->shiftKnotsForDuration(i, plus_delta, u_plus);
            this->shiftKnotsForDuration(i, minus_delta, u_minus);

            const auto span_range = affectedEnergySpanRangeByDurationFixed(i);
            const double cost_p =
                getEnergyForKnotsRange(u_plus, span_range.first, span_range.second);
            const double cost_m =
                getEnergyForKnotsRange(u_minus, span_range.first, span_range.second);

            gdT_direct(i) = (cost_p - cost_m) / (plus_delta - minus_delta);
        }
    }

    inline void propagateGradAnalytic(const Eigen::MatrixXd &gdC,
                                      const Eigen::VectorXd &gdT_direct,
                                      Eigen::MatrixXd &gradByPoints,
                                      Eigen::VectorXd &gradByTimes) const
    {
        this->ensureKnotJacobian();
        const int M = static_cast<int>(this->durations_.size());
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes = gdT_direct;

        Eigen::MatrixXd adjGrad = gdC;
        this->A.solveAdj(adjGrad);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjGrad.row(S + i);
        }

        auto accumulate_row_adj = [&](const int row,
                                      const int d_ord,
                                      const int span,
                                      const double t,
                                      const Eigen::VectorXd &dt_dT)
        {
            Eigen::Matrix<double, P + 1, 1> vals_dummy;
            Eigen::MatrixXd dvals_dT;
            this->evalLocalBasisAndTimeGrad(d_ord, span, t, dt_dT,
                                            vals_dummy, dvals_dT);

            Eigen::MatrixXd rowVecByT = Eigen::MatrixXd::Zero(M, Dim);
            for (int j = 0; j <= P; ++j)
            {
                const int gidx = span - P + j;
                for (int m = 0; m < M; ++m)
                {
                    rowVecByT.row(m) +=
                        dvals_dT(j, m) * this->control_points.row(gidx);
                }
            }

            for (int m = 0; m < M; ++m)
            {
                gradByTimes(m) -= adjGrad.row(row).dot(rowVecByT.row(m));
            }
        };

        int row = 0;
        for (int d = 0; d < S; ++d)
        {
            accumulate_row_adj(row++, d, P, this->knots(P), Eigen::VectorXd::Zero(M));
        }
        for (int i = 1; i < M; ++i)
        {
            accumulate_row_adj(row++, 0, P + i, this->knots(P + i),
                               this->knotJacobian.row(P + i).transpose());
        }
        for (int d = S - 1; d >= 0; --d)
        {
            accumulate_row_adj(row++, d, this->N_c - 1, this->knots(this->N_c),
                               this->knotJacobian.row(this->N_c).transpose());
        }
    }

    inline void propagateGradFiniteDiff(const Eigen::MatrixXd &gdC,
                                        const Eigen::VectorXd &gdT_direct,
                                        const Eigen::VectorXd &T,
                                        Eigen::MatrixXd &gradByPoints,
                                        Eigen::VectorXd &gradByTimes) const
    {
        const int M = static_cast<int>(T.size());
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes.resize(M);

        Eigen::MatrixXd adjGrad = gdC;
        this->A.solveAdj(adjGrad);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjGrad.row(S + i);
        }

        Eigen::VectorXd u_plus = this->knots;
        Eigen::VectorXd u_minus = this->knots;
        for (int i = 0; i < M; ++i)
        {
            u_plus = this->knots;
            u_minus = this->knots;

            double plus_delta = 0.0;
            double minus_delta = 0.0;
            this->finiteDiffDeltasForDuration(T(i), plus_delta, minus_delta);
            const double denom = plus_delta - minus_delta;
            this->shiftKnotsForDuration(i, plus_delta, u_plus);
            this->shiftKnotsForDuration(i, minus_delta, u_minus);

            double plus_adj = 0.0;
            double minus_adj = 0.0;
            for (int row = 0; row < this->N_c; ++row)
            {
                if (!isConstraintRowAffectedByDurationFixed(row, i))
                {
                    continue;
                }
                plus_adj += dotSystemRowWithControlAndAdjFixed(
                    row, u_plus, this->control_points, adjGrad);
                minus_adj += dotSystemRowWithControlAndAdjFixed(
                    row, u_minus, this->control_points, adjGrad);
            }
            gradByTimes(i) = gdT_direct(i) - (plus_adj - minus_adj) / denom;
        }
    }

    inline void propagateEnergyGradFiniteDiffFull(const Eigen::MatrixXd &gdC,
                                                  const Eigen::VectorXd &T,
                                                  Eigen::MatrixXd &gradByPoints,
                                                  Eigen::VectorXd &gradByTimes) const
    {
        const int M = static_cast<int>(T.size());
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes.resize(M);

        Eigen::MatrixXd adjGrad = gdC;
        this->A.solveAdj(adjGrad);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjGrad.row(S + i);
        }

        Eigen::VectorXd u_plus = this->knots;
        Eigen::VectorXd u_minus = this->knots;
        BandedSystem A_plus;
        BandedSystem A_minus;
        for (int i = 0; i < M; ++i)
        {
            u_plus = this->knots;
            u_minus = this->knots;

            double plus_delta = 0.0;
            double minus_delta = 0.0;
            this->finiteDiffDeltasForDuration(T(i), plus_delta, minus_delta);
            const double denom = plus_delta - minus_delta;
            this->shiftKnotsForDuration(i, plus_delta, u_plus);
            this->shiftKnotsForDuration(i, minus_delta, u_minus);

            const double cost_p = getEnergyForKnots(u_plus);
            const double cost_m = getEnergyForKnots(u_minus);

            buildSystemMatrixFixed(M, u_plus, A_plus);
            buildSystemMatrixFixed(M, u_minus, A_minus);
            const double plus_adj =
                A_plus.dotMultiply(this->control_points, adjGrad);
            const double minus_adj =
                A_minus.dotMultiply(this->control_points, adjGrad);

            gradByTimes(i) =
                (cost_p - cost_m - plus_adj + minus_adj) / denom;
        }
    }

    inline void propagateEnergyGradFiniteDiff(const Eigen::MatrixXd &gdC,
                                              const Eigen::VectorXd &T,
                                              Eigen::MatrixXd &gradByPoints,
                                              Eigen::VectorXd &gradByTimes) const
    {
#ifdef NUBS_USE_FULL_FD_PROPAGATION
        propagateEnergyGradFiniteDiffFull(gdC, T, gradByPoints, gradByTimes);
        return;
#else
        const int M = static_cast<int>(T.size());
        gradByPoints.resize(std::max(0, M - 1), Dim);
        gradByTimes.resize(M);

        Eigen::MatrixXd adjGrad = gdC;
        this->A.solveAdj(adjGrad);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjGrad.row(S + i);
        }

        Eigen::VectorXd u_plus = this->knots;
        Eigen::VectorXd u_minus = this->knots;
        for (int i = 0; i < M; ++i)
        {
            u_plus = this->knots;
            u_minus = this->knots;

            double plus_delta = 0.0;
            double minus_delta = 0.0;
            this->finiteDiffDeltasForDuration(T(i), plus_delta, minus_delta);
            const double denom = plus_delta - minus_delta;
            this->shiftKnotsForDuration(i, plus_delta, u_plus);
            this->shiftKnotsForDuration(i, minus_delta, u_minus);

            const auto span_range = affectedEnergySpanRangeByDurationFixed(i);
            const double cost_p =
                getEnergyForKnotsRange(u_plus, span_range.first, span_range.second);
            const double cost_m =
                getEnergyForKnotsRange(u_minus, span_range.first, span_range.second);

            double plus_adj = 0.0;
            double minus_adj = 0.0;
            for (int row = 0; row < this->N_c; ++row)
            {
                if (!isConstraintRowAffectedByDurationFixed(row, i))
                {
                    continue;
                }
                plus_adj += dotSystemRowWithControlAndAdjFixed(
                    row, u_plus, this->control_points, adjGrad);
                minus_adj += dotSystemRowWithControlAndAdjFixed(
                    row, u_minus, this->control_points, adjGrad);
            }

#ifdef NUBS_VALIDATE_LOCAL_FD
            const double full_cost_p = getEnergyForKnots(u_plus);
            const double full_cost_m = getEnergyForKnots(u_minus);
            BandedSystem A_plus;
            BandedSystem A_minus;
            buildSystemMatrixFixed(M, u_plus, A_plus);
            buildSystemMatrixFixed(M, u_minus, A_minus);
            const double full_plus_adj =
                A_plus.dotMultiply(this->control_points, adjGrad);
            const double full_minus_adj =
                A_minus.dotMultiply(this->control_points, adjGrad);
            const double local_diff = cost_p - cost_m;
            const double full_diff = full_cost_p - full_cost_m;
            const double scale = std::max(1.0, std::abs(full_diff));
            assert(std::abs(local_diff - full_diff) <= 1.0e-8 * scale);
            const double adj_scale =
                std::max(1.0, std::abs(full_plus_adj - full_minus_adj));
            assert(std::abs((plus_adj - minus_adj) -
                            (full_plus_adj - full_minus_adj)) <=
                   1.0e-8 * adj_scale);
#endif

            gradByTimes(i) =
                (cost_p - cost_m - plus_adj + minus_adj) / denom;
        }
#endif
    }

    inline void propagateEnergyGradFixedRatioFiniteDiff(
        const Eigen::MatrixXd &gdC,
        Eigen::MatrixXd &gradByPoints,
        double &gradByTotalDuration) const
    {
        const int M = static_cast<int>(this->durations_.size());
        gradByPoints.resize(std::max(0, M - 1), Dim);

        Eigen::MatrixXd adjGrad = gdC;
        this->A.solveAdj(adjGrad);
        for (int i = 0; i < M - 1; ++i)
        {
            gradByPoints.row(i) = adjGrad.row(S + i);
        }

        double plus_delta = 0.0;
        double minus_delta = 0.0;
        this->finiteDiffDeltasForFixedRatioTotalDuration(plus_delta,
                                                         minus_delta);
        const double denom = plus_delta - minus_delta;

        Eigen::VectorXd u_plus;
        Eigen::VectorXd u_minus;
        this->scaleKnotsForTotalDurationDelta(plus_delta, u_plus);
        this->scaleKnotsForTotalDurationDelta(minus_delta, u_minus);

        const double cost_p = getEnergyForKnots(u_plus);
        const double cost_m = getEnergyForKnots(u_minus);

        BandedSystem A_plus;
        BandedSystem A_minus;
        buildSystemMatrixFixed(M, u_plus, A_plus);
        buildSystemMatrixFixed(M, u_minus, A_minus);
        const double plus_adj =
            A_plus.dotMultiply(this->control_points, adjGrad);
        const double minus_adj =
            A_minus.dotMultiply(this->control_points, adjGrad);

        gradByTotalDuration =
            (cost_p - cost_m - plus_adj + minus_adj) / denom;
    }

    inline void getEnergyAndAnalyticGrad(double &cost,
                                         Eigen::MatrixXd &gradByPoints,
                                         Eigen::VectorXd &gradByTimes) const
    {
        // This validation path uses dense knot sensitivity through BasisJetCache.
        // Prefer getEnergyAndGrad() for optimization loops.
        Eigen::MatrixXd gdC;
        getEnergyPartialGradByCoeffs(cost, gdC);
        Eigen::VectorXd gdT_direct;
        getEnergyPartialGradByTimesAnalytic(gdT_direct);
        propagateGradAnalytic(gdC, gdT_direct, gradByPoints, gradByTimes);
    }

    // Fixed-order production API. It uses the exact derivative-control energy
    // kernel with local time sensitivities, hiding the generic parent
    // implementation whose scalar Dual loop revisits a span per duration.
    inline void getEnergyPartialGradByTimesLocalAD(
        Eigen::VectorXd &gdT_direct) const
    {
        double ignored_cost = 0.0;
        Eigen::MatrixXd ignored_gdC;
        getEnergyPartialGradDerivativeControl(ignored_cost, ignored_gdC,
                                              gdT_direct);
    }

    inline void getEnergyPartialGradLocalAD(double &cost,
                                            Eigen::MatrixXd &gdC,
                                            Eigen::VectorXd &gdT_direct) const
    {
        getEnergyPartialGradDerivativeControl(cost, gdC, gdT_direct);
    }

    inline void getEnergyAndLocalADGrad(double &cost,
                                        Eigen::MatrixXd &gradByPoints,
                                        Eigen::VectorXd &gradByTimes) const
    {
        Eigen::MatrixXd gdC;
        Eigen::VectorXd gdT_direct;
        getEnergyPartialGradDerivativeControl(cost, gdC, gdT_direct);
        propagateEnergyGradDerivativeControl(gdC, gdT_direct,
                                              gradByPoints, gradByTimes);
    }

    inline void getEnergyAndGrad(double &cost,
                                 Eigen::MatrixXd &gradByPoints,
                                 Eigen::VectorXd &gradByTimes) const
    {
        getEnergyAndLocalADGrad(cost, gradByPoints, gradByTimes);
    }

    inline void getEnergyAndFiniteDiffGrad(double &cost,
                                           Eigen::MatrixXd &gradByPoints,
                                           Eigen::VectorXd &gradByTimes) const
    {
        Eigen::MatrixXd gdC;
        getEnergyPartialGradByCoeffs(cost, gdC);
        propagateEnergyGradFiniteDiff(gdC, this->durations_,
                                      gradByPoints, gradByTimes);
    }

    inline void getEnergyAndFiniteDiffGradFull(double &cost,
                                               Eigen::MatrixXd &gradByPoints,
                                               Eigen::VectorXd &gradByTimes) const
    {
        Eigen::MatrixXd gdC;
        getEnergyPartialGradByCoeffs(cost, gdC);
        propagateEnergyGradFiniteDiffFull(gdC, this->durations_,
                                          gradByPoints, gradByTimes);
    }

    inline void getEnergyAndFixedRatioGrad(double &cost,
                                           Eigen::MatrixXd &gradByPoints,
                                           double &gradByTotalDuration) const
    {
        Eigen::MatrixXd gdC;
        getEnergyPartialGradByCoeffs(cost, gdC);
        propagateEnergyGradFixedRatioFiniteDiff(gdC, gradByPoints,
                                                gradByTotalDuration);
    }

    inline void getEnergyAndTotalDurationGrad(double &cost,
                                              Eigen::MatrixXd &gradByPoints,
                                              double &gradByTotalDuration) const
    {
        getEnergyAndFixedRatioGrad(cost, gradByPoints, gradByTotalDuration);
    }

    inline void getEnergyAndUniformTimeGrad(double &cost,
                                            Eigen::MatrixXd &gradByPoints,
                                            double &gradByTotalDuration) const
    {
        getEnergyAndFixedRatioGrad(cost, gradByPoints, gradByTotalDuration);
    }
};

template <int Dim, int S>
class UniformNUBSTrajectoryT final : public NUBSTrajectoryT<Dim, S>
{
    static_assert(S >= 2 && S <= 4,
                  "UniformNUBSTrajectoryT currently supports S = 2, 3, and 4.");

private:
    using Parent = NUBSTrajectoryT<Dim, S>;
    static constexpr int P = 2 * S - 1;

    mutable bool fullSystemFactorized_ = false;
    mutable Eigen::Matrix<double, Eigen::Dynamic, Dim> uniformInnerPoints_;
    mutable Eigen::Matrix<double, Dim, S> uniformHeadState_;
    mutable Eigen::Matrix<double, Dim, S> uniformTailState_;
    mutable bool uniformProblemCached_ = false;
    mutable std::vector<Eigen::Matrix<double, P + 1, Dim>> uniformPolyCoeffs_;
    mutable bool uniformPolyCacheValid_ = false;

    struct BoundaryInverseCache
    {
        Eigen::Matrix<double, S, S> headInv;
        Eigen::Matrix<double, S, S> tailInv;
    };

    inline const BandedSystem &reducedMatrixForPieceNum(const int M) const
    {
        static std::map<int, BandedSystem> cache;
        const auto iter = cache.find(M);
        if (iter != cache.end())
        {
            return iter->second;
        }

        const int n = M - 1;
        BandedSystem reducedA;
        reducedA.create(n, S - 1, S);
        if (n > 0)
        {
            const int nc = this->getCtrlPtNum(M);
            const Eigen::VectorXd unitT = Eigen::VectorXd::Ones(M);
            const Eigen::VectorXd unitKnots = this->generateKnots(unitT, nc);
            Eigen::Matrix<double, P + 1, P + 1> ders;

            for (int i = 1; i < M; ++i)
            {
                const int row = i - 1;
                const int span = P + i;
                this->dersBasisFuns(0, span, unitKnots(span),
                                    unitKnots, ders);
                for (int j = 0; j <= P; ++j)
                {
                    const int gidx = span - P + j;
                    if (gidx >= S && gidx < nc - S)
                    {
                        reducedA(row, gidx - S) = ders(0, j);
                    }
                }
            }
            reducedA.factorizeLU();
        }

        const auto inserted = cache.emplace(M, std::move(reducedA));
        return inserted.first->second;
    }

    inline const BoundaryInverseCache &boundaryInverseForPieceNum(const int M) const
    {
        static std::map<int, BoundaryInverseCache> cache;
        const auto iter = cache.find(M);
        if (iter != cache.end())
        {
            return iter->second;
        }

        const int nc = this->getCtrlPtNum(M);
        const Eigen::VectorXd unitT = Eigen::VectorXd::Ones(M);
        const Eigen::VectorXd unitKnots = this->generateKnots(unitT, nc);
        Eigen::Matrix<double, P + 1, P + 1> ders;
        Eigen::Matrix<double, S, S> mat;
        BoundaryInverseCache value;

        mat.setZero();
        for (int d = 0; d < S; ++d)
        {
            this->dersBasisFuns(d, P, unitKnots(P), unitKnots, ders);
            for (int j = 0; j < S; ++j)
            {
                mat(d, j) = ders(d, j);
            }
        }
        value.headInv = mat.inverse();

        mat.setZero();
        const int firstLocalCol = nc - 1 - P;
        for (int d = 0; d < S; ++d)
        {
            this->dersBasisFuns(d, nc - 1, unitKnots(nc), unitKnots, ders);
            for (int j = 0; j < S; ++j)
            {
                const int gidx = nc - S + j;
                mat(d, j) = ders(d, gidx - firstLocalCol);
            }
        }
        value.tailInv = mat.inverse();

        const auto inserted = cache.emplace(M, value);
        return inserted.first->second;
    }

    inline void solveBoundaryControls(
        const Eigen::MatrixXd &headState,
        const Eigen::MatrixXd &tailState,
        Eigen::MatrixXd &controls) const
    {
        const BoundaryInverseCache &cache =
            boundaryInverseForPieceNum(static_cast<int>(this->durations_.size()));
        Eigen::Matrix<double, S, Dim> rhs;

        double h_power = 1.0;
        for (int d = 0; d < S; ++d)
        {
            rhs.row(d) = h_power * headState.col(d).transpose();
            h_power *= this->durations_(0);
        }
        controls.topRows(S) = cache.headInv * rhs;

        h_power = 1.0;
        for (int d = 0; d < S; ++d)
        {
            rhs.row(d) = h_power * tailState.col(d).transpose();
            h_power *= this->durations_(0);
        }
        controls.bottomRows(S) = cache.tailInv * rhs;
    }

    inline void ensureFullSystemFactorized() const
    {
        if (fullSystemFactorized_)
        {
            return;
        }
        auto *self = const_cast<UniformNUBSTrajectoryT *>(this);
        const int M = static_cast<int>(this->durations_.size());
        self->buildSystemMatrixA(M, this->knots, self->A);
        self->A.factorizeLU();
        fullSystemFactorized_ = true;
    }

    static inline double factorialValue(const int n)
    {
        double result = 1.0;
        for (int i = 2; i <= n; ++i)
        {
            result *= static_cast<double>(i);
        }
        return result;
    }

    inline int locateUniformPiece(double &t, double &tau) const
    {
        if (t <= 0.0)
        {
            t = 0.0;
        }

        const double total_duration = this->getTotalDuration();
        if (t >= total_duration)
        {
            t = std::max(0.0, total_duration - 1.0e-12);
        }

        const int M = static_cast<int>(this->durations_.size());
        const double h = total_duration / static_cast<double>(M);
        int piece = static_cast<int>(t / h);
        piece = std::min(std::max(piece, 0), M - 1);
        tau = t - h * static_cast<double>(piece);
        return piece;
    }

    inline void buildUniformPolynomialCache() const
    {
        const int M = static_cast<int>(this->durations_.size());
        uniformPolyCoeffs_.assign(
            M, Eigen::Matrix<double, P + 1, Dim>::Zero());

        Eigen::Matrix<double, P + 1, P + 1> ders;
        for (int piece = 0; piece < M; ++piece)
        {
            const int span = P + piece;
            const double t0 = this->knots(span);
            this->dersBasisFuns(P, span, t0, this->knots, ders);

            auto &coeffs = uniformPolyCoeffs_[piece];
            coeffs.setZero();
            for (int k = 0; k <= P; ++k)
            {
                const double inv_factorial = 1.0 / factorialValue(k);
                for (int j = 0; j <= P; ++j)
                {
                    coeffs.row(k) +=
                        inv_factorial * ders(k, j) *
                        this->control_points.row(span - P + j);
                }
            }
        }

        uniformPolyCacheValid_ = true;
    }

    inline std::array<Eigen::Matrix<double, Dim, 1>, 5>
    evaluateCachedDerivativesToSnap(double t) const
    {
        std::array<Eigen::Matrix<double, Dim, 1>, 5> values;
        for (auto &value : values)
        {
            value.setZero();
        }

        double tau = 0.0;
        const int piece = locateUniformPiece(t, tau);
        const auto &coeffs = uniformPolyCoeffs_[piece];

        values[0] = coeffs.row(P).transpose();
        for (int k = P - 1; k >= 0; --k)
        {
            for (int d = 4; d >= 1; --d)
            {
                values[d] =
                    tau * values[d] +
                    static_cast<double>(d) * values[d - 1];
            }
            values[0] = tau * values[0] + coeffs.row(k).transpose();
        }
        return values;
    }

public:
    static constexpr int SystemOrder = S;
    static constexpr int Degree = P;

    UniformNUBSTrajectoryT() : Parent() {}

    using Parent::generateFixedRatio;
    using Parent::generateWithTotalDuration;

    inline double getReducedUniformSystemConditionNumber() const
    {
        const int M = static_cast<int>(this->durations_.size());
        if (M <= 1)
        {
            return 1.0;
        }
        return Parent::conditionNumber2(reducedMatrixForPieceNum(M));
    }

    inline void generate(const Eigen::MatrixXd &P_inner,
                         const Eigen::MatrixXd &headState,
                         const Eigen::MatrixXd &tailState,
                         const Eigen::VectorXd &T,
                         Eigen::MatrixXd &P_full)
    {
        Parent::generate(P_inner, headState, tailState, T, P_full);
        fullSystemFactorized_ = true;
        uniformProblemCached_ = false;
        uniformPolyCacheValid_ = false;
    }

    inline void generateUniform(const Eigen::MatrixXd &P_inner,
                                const Eigen::MatrixXd &headState,
                                const Eigen::MatrixXd &tailState,
                                const double total_duration,
                                Eigen::MatrixXd &P_full)
    {
        const int M = static_cast<int>(P_inner.rows()) + 1;
        const Eigen::VectorXd T = this->uniformDurations(M, total_duration);
        this->validateGenerateInputs(P_inner, headState, tailState, T);

        this->N_c = this->getCtrlPtNum(M);
        this->durations_ = T;
        this->knots = this->generateKnots(T, this->N_c);
        this->knotJacobian.resize(0, 0);
        P_full.setZero(this->N_c, Dim);

        solveBoundaryControls(headState, tailState, P_full);
        const BandedSystem &reducedA = reducedMatrixForPieceNum(M);

        const int n = M - 1;
        if (n > 0)
        {
            Eigen::Matrix<double, Eigen::Dynamic, Dim> b(n, Dim);
            Eigen::Matrix<double, P + 1, P + 1> ders;
            for (int i = 1; i < M; ++i)
            {
                const int row = i - 1;
                const int span = P + i;
                b.row(row) = P_inner.row(row);
                this->dersBasisFuns(0, span, this->knots(span),
                                    this->knots, ders);
                for (int j = 0; j <= P; ++j)
                {
                    const int gidx = span - P + j;
                    if (gidx < S || gidx >= this->N_c - S)
                    {
                        b.row(row) -= ders(0, j) * P_full.row(gidx);
                    }
                }
            }

            reducedA.solve(b);
            P_full.middleRows(S, n) = b;
        }

        this->control_points = P_full;
        fullSystemFactorized_ = false;
        uniformProblemCached_ = true;
        uniformInnerPoints_ = P_inner;
        uniformHeadState_ = headState.template leftCols<S>();
        uniformTailState_ = tailState.template leftCols<S>();
        uniformPolyCacheValid_ = false;
    }

    inline void prepareEvaluationCache() const
    {
        if (!uniformPolyCacheValid_)
        {
            buildUniformPolynomialCache();
        }
    }

    inline int findUniformSpan(double &t) const
    {
        if (t <= 0.0)
        {
            t = 0.0;
            return P;
        }

        const double total_duration = this->getTotalDuration();
        if (t >= total_duration)
        {
            t = std::max(0.0, total_duration - 1.0e-12);
        }

        const int M = static_cast<int>(this->durations_.size());
        const double h = total_duration / static_cast<double>(M);
        int piece = static_cast<int>(t / h);
        piece = std::min(std::max(piece, 0), M - 1);
        return P + piece;
    }

    inline Eigen::Matrix<double, Dim, 1> evaluate(double t,
                                                  const int d_ord = 0) const
    {
        if (d_ord > P)
        {
            return Eigen::Matrix<double, Dim, 1>::Zero();
        }
        if (uniformPolyCacheValid_ && d_ord <= P)
        {
            double tau = 0.0;
            const int piece = locateUniformPiece(t, tau);
            const auto &coeffs = uniformPolyCoeffs_[piece];

            Eigen::Matrix<double, Dim, 1> res =
                Eigen::Matrix<double, Dim, 1>::Zero();
            if (d_ord <= P)
            {
                res = coeffs.row(P).transpose();
                double multiplier = 1.0;
                for (int m = 0; m < d_ord; ++m)
                {
                    multiplier *= static_cast<double>(P - m);
                }
                res *= multiplier;

                for (int k = P - 1; k >= d_ord; --k)
                {
                    multiplier = 1.0;
                    for (int m = 0; m < d_ord; ++m)
                    {
                        multiplier *= static_cast<double>(k - m);
                    }
                    res = tau * res + multiplier * coeffs.row(k).transpose();
                }
            }
            return res;
        }

        if (d_ord <= P)
        {
            prepareEvaluationCache();
            return evaluate(t, d_ord);
        }

        const int span = findUniformSpan(t);
        Eigen::Matrix<double, P + 1, P + 1> ders;
        this->dersBasisFuns(d_ord, span, t, this->knots, ders);

        Eigen::Matrix<double, Dim, 1> res =
            Eigen::Matrix<double, Dim, 1>::Zero();
        for (int j = 0; j <= P; ++j)
        {
            res += ders(d_ord, j) *
                   this->control_points.row(span - P + j).transpose();
        }
        return res;
    }

    inline void evaluatePVA(double t,
                            Eigen::Matrix<double, Dim, 1> &pos,
                            Eigen::Matrix<double, Dim, 1> &vel,
                            Eigen::Matrix<double, Dim, 1> &acc) const
    {
        if (!uniformPolyCacheValid_)
        {
            prepareEvaluationCache();
        }

        const auto values = evaluateCachedDerivativesToSnap(t);
        pos = values[0];
        vel = values[1];
        acc = values[2];
    }

    inline void evaluatePVAJ(double t,
                             Eigen::Matrix<double, Dim, 1> &pos,
                             Eigen::Matrix<double, Dim, 1> &vel,
                             Eigen::Matrix<double, Dim, 1> &acc,
                             Eigen::Matrix<double, Dim, 1> &jerk) const
    {
        if (!uniformPolyCacheValid_)
        {
            prepareEvaluationCache();
        }

        const auto values = evaluateCachedDerivativesToSnap(t);
        pos = values[0];
        vel = values[1];
        acc = values[2];
        jerk = values[3];
    }

    inline void evaluatePVAJS(double t,
                              Eigen::Matrix<double, Dim, 1> &pos,
                              Eigen::Matrix<double, Dim, 1> &vel,
                              Eigen::Matrix<double, Dim, 1> &acc,
                              Eigen::Matrix<double, Dim, 1> &jerk,
                              Eigen::Matrix<double, Dim, 1> &snap) const
    {
        if (!uniformPolyCacheValid_)
        {
            prepareEvaluationCache();
        }

        const auto values = evaluateCachedDerivativesToSnap(t);
        pos = values[0];
        vel = values[1];
        acc = values[2];
        jerk = values[3];
        snap = values[4];
    }

    inline Eigen::Matrix<double, Dim, 1> getDerivative(
        const double t,
        const int derivative) const
    {
        return evaluate(t, derivative);
    }

    inline Eigen::Matrix<double, Dim, 1> getPos(const double t) const
    {
        return evaluate(t, 0);
    }

    inline Eigen::Matrix<double, Dim, 1> getVel(const double t) const
    {
        return evaluate(t, 1);
    }

    inline Eigen::Matrix<double, Dim, 1> getAcc(const double t) const
    {
        return evaluate(t, 2);
    }

    inline Eigen::Matrix<double, Dim, 1> getJer(const double t) const
    {
        return evaluate(t, 3);
    }

    inline Eigen::Matrix<double, Dim, 1> getJerk(const double t) const
    {
        return getJer(t);
    }

    inline Eigen::Matrix<double, Dim, 1> getSnap(const double t) const
    {
        return evaluate(t, 4);
    }

    inline void getEnergyAndFiniteDiffGrad(double &cost,
                                           Eigen::MatrixXd &gradByPoints,
                                           Eigen::VectorXd &gradByTimes) const
    {
        ensureFullSystemFactorized();
        Parent::getEnergyAndFiniteDiffGrad(cost, gradByPoints, gradByTimes);
    }

    inline void getEnergyAndLocalADGrad(double &cost,
                                        Eigen::MatrixXd &gradByPoints,
                                        Eigen::VectorXd &gradByTimes) const
    {
        ensureFullSystemFactorized();
        Parent::getEnergyAndLocalADGrad(cost, gradByPoints, gradByTimes);
    }

    inline void getEnergyAndGrad(double &cost,
                                 Eigen::MatrixXd &gradByPoints,
                                 Eigen::VectorXd &gradByTimes) const
    {
        getEnergyAndLocalADGrad(cost, gradByPoints, gradByTimes);
    }

    inline void getEnergyAndFiniteDiffGradFull(double &cost,
                                               Eigen::MatrixXd &gradByPoints,
                                               Eigen::VectorXd &gradByTimes) const
    {
        ensureFullSystemFactorized();
        Parent::getEnergyAndFiniteDiffGradFull(cost, gradByPoints, gradByTimes);
    }

    inline void getEnergyAndFixedRatioGrad(double &cost,
                                           Eigen::MatrixXd &gradByPoints,
                                           double &gradByTotalDuration) const
    {
        ensureFullSystemFactorized();
        Parent::getEnergyAndFixedRatioGrad(cost, gradByPoints,
                                           gradByTotalDuration);
    }

    inline void getEnergyAndTotalDurationGrad(double &cost,
                                              Eigen::MatrixXd &gradByPoints,
                                              double &gradByTotalDuration) const
    {
        getEnergyAndUniformTimeGrad(cost, gradByPoints, gradByTotalDuration);
    }

    inline void getEnergyAndUniformTimeGrad(double &cost,
                                            Eigen::MatrixXd &gradByPoints,
                                            double &gradByTotalDuration) const
    {
        if (!uniformProblemCached_)
        {
            ensureFullSystemFactorized();
            Parent::getEnergyAndFixedRatioGrad(cost, gradByPoints,
                                               gradByTotalDuration);
            return;
        }

        Eigen::MatrixXd gdC;
        this->getEnergyPartialGradByCoeffs(cost, gdC);
        const int M = static_cast<int>(this->durations_.size());
        const int n = M - 1;
        gradByPoints.resize(n, Dim);
        if (n > 0)
        {
            Eigen::MatrixXd adj = gdC.middleRows(S, n);
            reducedMatrixForPieceNum(M).solveAdj(adj);
            gradByPoints = adj;
        }

        double plus_delta = 0.0;
        double minus_delta = 0.0;
        this->finiteDiffDeltasForFixedRatioTotalDuration(plus_delta,
                                                         minus_delta);

        UniformNUBSTrajectoryT plusTraj;
        UniformNUBSTrajectoryT minusTraj;
        Eigen::MatrixXd tmp;
        const double total_duration = this->getTotalDuration();
        plusTraj.generateUniform(uniformInnerPoints_, uniformHeadState_,
                                 uniformTailState_,
                                 total_duration + plus_delta, tmp);
        minusTraj.generateUniform(uniformInnerPoints_, uniformHeadState_,
                                  uniformTailState_,
                                  total_duration + minus_delta, tmp);
        gradByTotalDuration =
            (plusTraj.getEnergy() - minusTraj.getEnergy()) /
            (plus_delta - minus_delta);
    }
};

template <int Dim>
using CubicNUBS = NUBSTrajectoryT<Dim, 2>;

template <int Dim>
using QuinticNUBS = NUBSTrajectoryT<Dim, 3>;

template <int Dim>
using SepticNUBS = NUBSTrajectoryT<Dim, 4>;

template <int Dim>
using NUBSTrajectoryS2 = CubicNUBS<Dim>;

template <int Dim>
using NUBSTrajectoryS3 = QuinticNUBS<Dim>;

template <int Dim>
using NUBSTrajectoryS4 = SepticNUBS<Dim>;

template <int Dim, int MaxP = 7>
using UniformNUBSTrajectory = NUBSTrajectory<Dim, MaxP>;

template <int Dim, int S>
using UBSTrajectoryT = UniformNUBSTrajectoryT<Dim, S>;

template <int Dim>
using UniformCubicNUBS = UniformNUBSTrajectoryT<Dim, 2>;

template <int Dim>
using UniformQuinticNUBS = UniformNUBSTrajectoryT<Dim, 3>;

template <int Dim>
using UniformSepticNUBS = UniformNUBSTrajectoryT<Dim, 4>;

template <int Dim>
using CubicUBS = UBSTrajectoryT<Dim, 2>;

template <int Dim>
using QuinticUBS = UBSTrajectoryT<Dim, 3>;

template <int Dim>
using SepticUBS = UBSTrajectoryT<Dim, 4>;

template <int Dim>
using UniformNUBSTrajectoryS2 = UniformCubicNUBS<Dim>;

template <int Dim>
using UniformNUBSTrajectoryS3 = UniformQuinticNUBS<Dim>;

template <int Dim>
using UniformNUBSTrajectoryS4 = UniformSepticNUBS<Dim>;

} // namespace nubs

#endif
