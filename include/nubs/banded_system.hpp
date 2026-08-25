#ifndef NUBS_BANDED_SYSTEM_HPP
#define NUBS_BANDED_SYSTEM_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nubs
{

class BandedSystem
{
public:
    struct FactorizationStats
    {
        double minimum_abs_pivot = std::numeric_limits<double>::infinity();
        int near_zero_pivot_count = 0;
    };

private:
    int N = 0;
    int lowerBw = 0;
    int upperBw = 0;
    std::vector<double> data;
    std::vector<double> unfactorized_data;
    FactorizationStats factorization_stats_;

public:
    inline void create(const int &n, const int &p, const int &q)
    {
        if (n < 0 || p < 0 || q < 0)
        {
            throw std::runtime_error("BandedSystem::create(): invalid dimensions.");
        }
        if (N == n && lowerBw == p && upperBw == q)
        {
            std::fill(data.begin(), data.end(), 0.0);
        }
        else
        {
            N = n;
            lowerBw = p;
            upperBw = q;
            data.assign(N * (lowerBw + upperBw + 1), 0.0);
        }
        factorization_stats_ = FactorizationStats{};
        unfactorized_data.clear();
    }

    inline int rows() const { return N; }
    inline int lowerBandwidth() const { return lowerBw; }
    inline int upperBandwidth() const { return upperBw; }
    inline const FactorizationStats &factorizationStats() const
    {
        return factorization_stats_;
    }

    // Access the coefficient prior to LU factorisation when that snapshot is
    // available; otherwise return the current matrix coefficient. This keeps
    // numerical diagnostics independent of the in-place factor storage.
    inline double originalCoefficient(const int i, const int j) const
    {
        if (i < 0 || i >= N || j < 0 || j >= N ||
            i - j < -upperBw || i - j > lowerBw)
        {
            return 0.0;
        }
        const std::size_t offset =
            static_cast<std::size_t>(i - j + upperBw) * N + j;
        return unfactorized_data.empty() ? data[offset] : unfactorized_data[offset];
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
        factorization_stats_ = FactorizationStats{};
        unfactorized_data = data;
        int iM, jM;
        double cVl;
        for (int k = 0; k <= N - 2; ++k)
        {
            iM = std::min(k + lowerBw, N - 1);
            cVl = operator()(k, k);
            factorization_stats_.minimum_abs_pivot =
                std::min(factorization_stats_.minimum_abs_pivot, std::abs(cVl));
            if (std::abs(cVl) < 1.0e-14)
            {
                ++factorization_stats_.near_zero_pivot_count;
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
        if (N > 0)
        {
            const double final_pivot = std::abs(operator()(N - 1, N - 1));
            factorization_stats_.minimum_abs_pivot =
                std::min(factorization_stats_.minimum_abs_pivot, final_pivot);
            if (final_pivot < 1.0e-14)
            {
                ++factorization_stats_.near_zero_pivot_count;
                throw std::runtime_error("BandedSystem::factorizeLU(): near-zero final pivot.");
            }
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

    // Infinity-norm relative residual of an original system after factorizeLU.
    // Keeping one banded copy is intentional: it makes numerical diagnostics
    // available without converting the system to a dense matrix.
    template <typename EIGENMATX, typename EIGENMATB>
    inline double relativeResidualInfinityNorm(const EIGENMATX &x,
                                               const EIGENMATB &b) const
    {
        if (unfactorized_data.empty() && N > 0)
        {
            throw std::runtime_error(
                "BandedSystem::relativeResidualInfinityNorm(): factorizeLU() has not been called.");
        }
        if (x.rows() != N || b.rows() != N || x.cols() != b.cols())
        {
            throw std::runtime_error(
                "BandedSystem::relativeResidualInfinityNorm(): incompatible matrix dimensions.");
        }

        double residual_norm = 0.0;
        double rhs_norm = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const int j_start = std::max(0, i - lowerBw);
            const int j_end = std::min(N - 1, i + upperBw);
            for (int column = 0; column < x.cols(); ++column)
            {
                double ax = 0.0;
                for (int j = j_start; j <= j_end; ++j)
                {
                    ax += unfactorized_data[(i - j + upperBw) * N + j] *
                          x(j, column);
                }
                residual_norm = std::max(residual_norm,
                                         std::abs(ax - b(i, column)));
                rhs_norm = std::max(rhs_norm, std::abs(b(i, column)));
            }
        }
        return residual_norm / std::max(1.0, rhs_norm);
    }
};

} // namespace nubs

#endif
