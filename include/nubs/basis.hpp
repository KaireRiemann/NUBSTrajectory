#ifndef NUBS_BASIS_HPP
#define NUBS_BASIS_HPP

#include "nubs/dual.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace nubs
{
namespace basis
{

// Algorithm A2.3 from The NURBS Book, expressed over a scalar type.  Branch
// decisions use the primal value so the differentiated program has the same
// knot-multiplicity behaviour as the double implementation.
template <typename Scalar, int MaxP, typename KnotAccessor>
inline void dersBasisFuns(const int degree,
                          const int derivative_order,
                          const int span,
                          const Scalar &t,
                          const KnotAccessor &knots,
                          std::array<std::array<Scalar, MaxP + 1>, MaxP + 1> &ders)
{
    const int n = std::min(derivative_order, degree);
    for (auto &row : ders)
    {
        row.fill(Scalar{});
    }

    std::array<std::array<Scalar, MaxP + 1>, MaxP + 1> ndu{};
    std::array<Scalar, MaxP + 1> left{};
    std::array<Scalar, MaxP + 1> right{};
    ndu[0][0] = Scalar(1.0);

    for (int j = 1; j <= degree; ++j)
    {
        left[j] = t - knots(span + 1 - j);
        right[j] = knots(span + j) - t;
        Scalar saved{};
        for (int r = 0; r < j; ++r)
        {
            ndu[j][r] = right[r + 1] + left[j - r];
            const Scalar temp =
                std::abs(ad::primal(ndu[j][r])) < 1.0e-15
                    ? Scalar{}
                    : ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    for (int j = 0; j <= degree; ++j)
    {
        ders[0][j] = ndu[j][degree];
    }

    std::array<std::array<Scalar, MaxP + 1>, 2> a{};
    for (int r = 0; r <= degree; ++r)
    {
        int s1 = 0;
        int s2 = 1;
        a[0].fill(Scalar{});
        a[0][0] = Scalar(1.0);
        for (int k = 1; k <= n; ++k)
        {
            Scalar d{};
            const int rk = r - k;
            const int pk = degree - k;
            a[s2].fill(Scalar{});
            if (r >= k)
            {
                const Scalar den = ndu[pk + 1][rk];
                a[s2][0] = std::abs(ad::primal(den)) < 1.0e-15
                               ? Scalar{}
                               : a[s1][0] / den;
                d = a[s2][0] * ndu[rk][pk];
            }
            const int j1 = rk >= -1 ? 1 : -rk;
            const int j2 = r - 1 <= pk ? k - 1 : degree - r;
            for (int j = j1; j <= j2; ++j)
            {
                const Scalar den = ndu[pk + 1][rk + j];
                a[s2][j] = std::abs(ad::primal(den)) < 1.0e-15
                               ? Scalar{}
                               : (a[s1][j] - a[s1][j - 1]) / den;
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk)
            {
                const Scalar den = ndu[pk + 1][r];
                a[s2][k] = std::abs(ad::primal(den)) < 1.0e-15
                               ? Scalar{}
                               : -a[s1][k - 1] / den;
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            std::swap(s1, s2);
        }
    }

    double factor = static_cast<double>(degree);
    for (int k = 1; k <= n; ++k)
    {
        for (int j = 0; j <= degree; ++j)
        {
            ders[k][j] *= Scalar(factor);
        }
        factor *= static_cast<double>(degree - k);
    }
}

// Fixed-degree counterpart used by NUBSTrajectoryT.  Keeping the degree in
// the template parameter permits the compiler to specialize the small
// (P+1)-sized recurrence for the production S=2/3/4 cases.
template <typename Scalar, int Degree, typename KnotAccessor>
inline void dersBasisFunsFixed(
    const int derivative_order,
    const int span,
    const Scalar &t,
    const KnotAccessor &knots,
    std::array<std::array<Scalar, Degree + 1>, Degree + 1> &ders)
{
    const int n = std::min(derivative_order, Degree);
    for (auto &row : ders)
    {
        row.fill(Scalar{});
    }

    std::array<std::array<Scalar, Degree + 1>, Degree + 1> ndu{};
    std::array<Scalar, Degree + 1> left{};
    std::array<Scalar, Degree + 1> right{};
    ndu[0][0] = Scalar(1.0);

    for (int j = 1; j <= Degree; ++j)
    {
        left[j] = t - knots(span + 1 - j);
        right[j] = knots(span + j) - t;
        Scalar saved{};
        for (int r = 0; r < j; ++r)
        {
            ndu[j][r] = right[r + 1] + left[j - r];
            const Scalar temp =
                std::abs(ad::primal(ndu[j][r])) < 1.0e-15
                    ? Scalar{}
                    : ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    for (int j = 0; j <= Degree; ++j)
    {
        ders[0][j] = ndu[j][Degree];
    }

    std::array<std::array<Scalar, Degree + 1>, 2> a{};
    for (int r = 0; r <= Degree; ++r)
    {
        int s1 = 0;
        int s2 = 1;
        a[0].fill(Scalar{});
        a[0][0] = Scalar(1.0);
        for (int k = 1; k <= n; ++k)
        {
            Scalar d{};
            const int rk = r - k;
            const int pk = Degree - k;
            a[s2].fill(Scalar{});
            if (r >= k)
            {
                const Scalar den = ndu[pk + 1][rk];
                a[s2][0] = std::abs(ad::primal(den)) < 1.0e-15
                               ? Scalar{}
                               : a[s1][0] / den;
                d = a[s2][0] * ndu[rk][pk];
            }
            const int j1 = rk >= -1 ? 1 : -rk;
            const int j2 = r - 1 <= pk ? k - 1 : Degree - r;
            for (int j = j1; j <= j2; ++j)
            {
                const Scalar den = ndu[pk + 1][rk + j];
                a[s2][j] = std::abs(ad::primal(den)) < 1.0e-15
                               ? Scalar{}
                               : (a[s1][j] - a[s1][j - 1]) / den;
                d += a[s2][j] * ndu[rk + j][pk];
            }
            if (r <= pk)
            {
                const Scalar den = ndu[pk + 1][r];
                a[s2][k] = std::abs(ad::primal(den)) < 1.0e-15
                               ? Scalar{}
                               : -a[s1][k - 1] / den;
                d += a[s2][k] * ndu[r][pk];
            }
            ders[k][r] = d;
            std::swap(s1, s2);
        }
    }

    double factor = static_cast<double>(Degree);
    for (int k = 1; k <= n; ++k)
    {
        for (int j = 0; j <= Degree; ++j)
        {
            ders[k][j] *= Scalar(factor);
        }
        factor *= static_cast<double>(Degree - k);
    }
}

} // namespace basis
} // namespace nubs

#endif
