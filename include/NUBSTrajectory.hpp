#ifndef NUBS_TRAJECTORY_HPP
#define NUBS_TRAJECTORY_HPP

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace nubs
{

class BandedSystem
{
private:
    int N = 0;
    int lowerBw = 0;
    int upperBw = 0;
    std::vector<double> data;

public:
    inline void create(const int &n, const int &p, const int &q)
    {
        if (N == n && lowerBw == p && upperBw == q)
        {
            std::fill(data.begin(), data.end(), 0.0);
            return;
        }

        N = n;
        lowerBw = p;
        upperBw = q;
        data.assign(N * (lowerBw + upperBw + 1), 0.0);
    }

    inline double &operator()(const int &i, const int &j)
    {
        return data[(i - j + upperBw) * N + j];
    }

    inline const double &operator()(const int &i, const int &j) const
    {
        return data[(i - j + upperBw) * N + j];
    }

    template <typename EIGENMAT>
    inline EIGENMAT multiply(const EIGENMAT &x) const
    {
        EIGENMAT res = EIGENMAT::Zero(N, x.cols());
        for (int i = 0; i < N; ++i)
        {
            const int j_start = std::max(0, i - lowerBw);
            const int j_end = std::min(N - 1, i + upperBw);
            for (int j = j_start; j <= j_end; ++j)
            {
                res.row(i) += operator()(i, j) * x.row(j);
            }
        }
        return res;
    }

    template <typename EIGENMATX, typename EIGENMATY>
    inline double dotMultiply(const EIGENMATX &x, const EIGENMATY &y) const
    {
        double res = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const int j_start = std::max(0, i - lowerBw);
            const int j_end = std::min(N - 1, i + upperBw);
            for (int j = j_start; j <= j_end; ++j)
            {
                res += operator()(i, j) * y.row(i).dot(x.row(j));
            }
        }
        return res;
    }

    inline void factorizeLU()
    {
        int iM, jM;
        double cVl;
        for (int k = 0; k <= N - 2; ++k)
        {
            iM = std::min(k + lowerBw, N - 1);
            cVl = operator()(k, k);
            if (std::abs(cVl) < 1.0e-14)
            {
                throw std::runtime_error("BandedSystem::factorizeLU(): near-zero pivot.");
            }
            for (int i = k + 1; i <= iM; ++i)
            {
                if (operator()(i, k) != 0.0)
                {
                    operator()(i, k) /= cVl;
                }
            }
            jM = std::min(k + upperBw, N - 1);
            for (int j = k + 1; j <= jM; ++j)
            {
                cVl = operator()(k, j);
                if (cVl != 0.0)
                {
                    for (int i = k + 1; i <= iM; ++i)
                    {
                        if (operator()(i, k) != 0.0)
                        {
                            operator()(i, j) -= operator()(i, k) * cVl;
                        }
                    }
                }
            }
        }
        if (N > 0 && std::abs(operator()(N - 1, N - 1)) < 1.0e-14)
        {
            throw std::runtime_error("BandedSystem::factorizeLU(): near-zero final pivot.");
        }
    }

    template <typename EIGENMAT>
    inline void solve(EIGENMAT &b) const
    {
        int iM;
        for (int j = 0; j <= N - 1; ++j)
        {
            iM = std::min(j + lowerBw, N - 1);
            for (int i = j + 1; i <= iM; ++i)
            {
                if (operator()(i, j) != 0.0)
                {
                    b.row(i) -= operator()(i, j) * b.row(j);
                }
            }
        }
        for (int j = N - 1; j >= 0; --j)
        {
            b.row(j) /= operator()(j, j);
            iM = std::max(0, j - upperBw);
            for (int i = iM; i <= j - 1; ++i)
            {
                if (operator()(i, j) != 0.0)
                {
                    b.row(i) -= operator()(i, j) * b.row(j);
                }
            }
        }
    }

    template <typename EIGENMAT>
    inline void solveAdj(EIGENMAT &b) const
    {
        int iM;
        for (int j = 0; j <= N - 1; ++j)
        {
            b.row(j) /= operator()(j, j);
            iM = std::min(j + upperBw, N - 1);
            for (int i = j + 1; i <= iM; ++i)
            {
                if (operator()(j, i) != 0.0)
                {
                    b.row(i) -= operator()(j, i) * b.row(j);
                }
            }
        }
        for (int j = N - 1; j >= 0; --j)
        {
            iM = std::max(0, j - lowerBw);
            for (int i = iM; i <= j - 1; ++i)
            {
                if (operator()(j, i) != 0.0)
                {
                    b.row(i) -= operator()(j, i) * b.row(j);
                }
            }
        }
    }
};

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

        A.factorizeLU();
        A.solve(b);
        P_full = b;
        control_points = P_full;
    }

    inline Eigen::Matrix<double, Dim, 1> evaluate(double t,
                                                  const int d_ord = 0) const
    {
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

    inline void getEnergyAndAnalyticGrad(double &cost,
                                         Eigen::MatrixXd &gradByPoints,
                                         Eigen::VectorXd &gradByTimes) const
    {
        // This validation path uses dense knot sensitivity through BasisJetCache.
        // Prefer getEnergyAndFiniteDiffGrad() for optimization loops.
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

    inline void getEnergyAndFiniteDiffGradFull(double &cost,
                                               Eigen::MatrixXd &gradByPoints,
                                               Eigen::VectorXd &gradByTimes) const
    {
        Eigen::MatrixXd gdC;
        getEnergyPartialGradByCoeffs(cost, gdC);
        propagateEnergyGradFiniteDiffFull(gdC, durations_, gradByPoints, gradByTimes);
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
class NUBSTrajectoryT final : public NUBSTrajectory<Dim, 2 * S - 1>
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

        this->A.factorizeLU();
        this->A.solve(b);
        P_full = b;
        this->control_points = P_full;
    }

    inline Eigen::Matrix<double, Dim, 1> evaluate(double t,
                                                  const int d_ord = 0) const
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
        return getEnergyForKnots(this->knots);
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

    inline void getEnergyAndAnalyticGrad(double &cost,
                                         Eigen::MatrixXd &gradByPoints,
                                         Eigen::VectorXd &gradByTimes) const
    {
        // This validation path uses dense knot sensitivity through BasisJetCache.
        // Prefer getEnergyAndFiniteDiffGrad() for optimization loops.
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

} // namespace nubs

#endif
