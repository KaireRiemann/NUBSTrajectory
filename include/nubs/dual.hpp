#ifndef NUBS_DUAL_HPP
#define NUBS_DUAL_HPP

#include <cmath>

namespace nubs
{
namespace ad
{

// A scalar dual number used for one local timing direction at a time.  The
// production timing gradient never allocates an M-dimensional derivative.
struct Dual
{
    double value = 0.0;
    double derivative = 0.0;

    constexpr Dual() = default;
    constexpr Dual(const double v, const double d = 0.0)
        : value(v), derivative(d)
    {
    }

    constexpr Dual &operator+=(const Dual &rhs)
    {
        value += rhs.value;
        derivative += rhs.derivative;
        return *this;
    }

    constexpr Dual &operator-=(const Dual &rhs)
    {
        value -= rhs.value;
        derivative -= rhs.derivative;
        return *this;
    }

    constexpr Dual &operator*=(const Dual &rhs)
    {
        derivative = derivative * rhs.value + value * rhs.derivative;
        value *= rhs.value;
        return *this;
    }

    constexpr Dual &operator/=(const Dual &rhs)
    {
        const double denominator = rhs.value * rhs.value;
        derivative = (derivative * rhs.value - value * rhs.derivative) /
                     denominator;
        value /= rhs.value;
        return *this;
    }
};

constexpr inline Dual operator+(Dual lhs, const Dual &rhs)
{
    lhs += rhs;
    return lhs;
}

constexpr inline Dual operator-(Dual lhs, const Dual &rhs)
{
    lhs -= rhs;
    return lhs;
}

constexpr inline Dual operator*(Dual lhs, const Dual &rhs)
{
    lhs *= rhs;
    return lhs;
}

constexpr inline Dual operator/(Dual lhs, const Dual &rhs)
{
    lhs /= rhs;
    return lhs;
}

constexpr inline Dual operator-(const Dual &x)
{
    return Dual(-x.value, -x.derivative);
}

constexpr inline bool operator==(const Dual &lhs, const Dual &rhs)
{
    return lhs.value == rhs.value && lhs.derivative == rhs.derivative;
}

constexpr inline bool operator!=(const Dual &lhs, const Dual &rhs)
{
    return !(lhs == rhs);
}

inline Dual sqrt(const Dual &x)
{
    const double root = std::sqrt(x.value);
    return Dual(root, x.derivative / (2.0 * root));
}

inline double primal(const double x)
{
    return x;
}

constexpr inline double primal(const Dual &x)
{
    return x.value;
}

} // namespace ad
} // namespace nubs

#endif
