#ifndef NUBS_DUAL_HPP
#define NUBS_DUAL_HPP

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

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

// A fixed-width forward-mode jet for the small timing stencil attached to a
// single B-spline span or system row.  Unlike an Eigen vector derivative, its
// storage is stack allocated and its width is independent of the trajectory
// length.  It is deliberately separate from Dual: Dual remains the compact
// one-direction reference implementation used by the generic runtime-order
// class, while LocalJet lets the fixed-order classes evaluate a local stencil
// once and scatter all of its timing derivatives afterwards.
template <int Width>
struct LocalJet
{
    static_assert(Width > 0, "LocalJet must carry at least one derivative.");

    double value = 0.0;
    std::array<double, Width> derivative{};

    constexpr LocalJet() = default;
    constexpr LocalJet(const double v) : value(v) {}

    constexpr LocalJet &operator+=(const LocalJet &rhs)
    {
        value += rhs.value;
        for (int i = 0; i < Width; ++i)
        {
            derivative[i] += rhs.derivative[i];
        }
        return *this;
    }

    constexpr LocalJet &operator-=(const LocalJet &rhs)
    {
        value -= rhs.value;
        for (int i = 0; i < Width; ++i)
        {
            derivative[i] -= rhs.derivative[i];
        }
        return *this;
    }

    constexpr LocalJet &operator*=(const LocalJet &rhs)
    {
        const double lhs_value = value;
        for (int i = 0; i < Width; ++i)
        {
            derivative[i] = derivative[i] * rhs.value +
                            lhs_value * rhs.derivative[i];
        }
        value = lhs_value * rhs.value;
        return *this;
    }

    constexpr LocalJet &operator/=(const LocalJet &rhs)
    {
        const double lhs_value = value;
        const double denominator = rhs.value * rhs.value;
        for (int i = 0; i < Width; ++i)
        {
            derivative[i] = (derivative[i] * rhs.value -
                             lhs_value * rhs.derivative[i]) /
                            denominator;
        }
        value = lhs_value / rhs.value;
        return *this;
    }
};

template <int Width>
constexpr inline LocalJet<Width> operator+(LocalJet<Width> lhs,
                                            const LocalJet<Width> &rhs)
{
    lhs += rhs;
    return lhs;
}

template <int Width>
constexpr inline LocalJet<Width> operator-(LocalJet<Width> lhs,
                                            const LocalJet<Width> &rhs)
{
    lhs -= rhs;
    return lhs;
}

template <int Width>
constexpr inline LocalJet<Width> operator*(LocalJet<Width> lhs,
                                            const LocalJet<Width> &rhs)
{
    lhs *= rhs;
    return lhs;
}

template <int Width>
constexpr inline LocalJet<Width> operator/(LocalJet<Width> lhs,
                                            const LocalJet<Width> &rhs)
{
    lhs /= rhs;
    return lhs;
}

template <int Width>
constexpr inline LocalJet<Width> operator*(LocalJet<Width> lhs,
                                            const double rhs)
{
    lhs.value *= rhs;
    for (int i = 0; i < Width; ++i)
    {
        lhs.derivative[i] *= rhs;
    }
    return lhs;
}

template <int Width>
constexpr inline LocalJet<Width> operator*(const double lhs,
                                            LocalJet<Width> rhs)
{
    return rhs * lhs;
}

template <int Width>
constexpr inline LocalJet<Width> operator/(LocalJet<Width> lhs,
                                            const double rhs)
{
    lhs.value /= rhs;
    for (int i = 0; i < Width; ++i)
    {
        lhs.derivative[i] /= rhs;
    }
    return lhs;
}

template <int Width>
constexpr inline LocalJet<Width> operator/(const double lhs,
                                            const LocalJet<Width> &rhs)
{
    return LocalJet<Width>(lhs) / rhs;
}

template <int Width>
constexpr inline LocalJet<Width> operator-(const LocalJet<Width> &x)
{
    LocalJet<Width> result;
    result.value = -x.value;
    for (int i = 0; i < Width; ++i)
    {
        result.derivative[i] = -x.derivative[i];
    }
    return result;
}

template <int Width>
constexpr inline bool operator==(const LocalJet<Width> &lhs,
                                 const LocalJet<Width> &rhs)
{
    return lhs.value == rhs.value && lhs.derivative == rhs.derivative;
}

template <int Width>
constexpr inline bool operator!=(const LocalJet<Width> &lhs,
                                 const LocalJet<Width> &rhs)
{
    return !(lhs == rhs);
}

template <int Width>
constexpr inline double primal(const LocalJet<Width> &x)
{
    return x.value;
}

class ReverseTape;

// A scalar reverse-mode value backed by a small, reusable local tape.  It is
// intended for one B-spline span or one constraint row: the number of output
// derivatives is then at most 2P-1, while the differentiated expression is
// scalar.  Reverse mode obtains that entire local gradient in one backward
// sweep instead of carrying a 2P-1-wide derivative through every operation.
struct Reverse
{
    ReverseTape *tape = nullptr;
    int node = -1;
    double constant = 0.0;

    constexpr Reverse() = default;
    constexpr Reverse(const double value) : constant(value) {}
    Reverse(ReverseTape *owner, int node_index) : tape(owner), node(node_index) {}

    double value() const;

    Reverse &operator+=(const Reverse &rhs);
    Reverse &operator-=(const Reverse &rhs);
    Reverse &operator*=(const Reverse &rhs);
    Reverse &operator/=(const Reverse &rhs);
};

class ReverseTape
{
    struct Node
    {
        double value = 0.0;
        int lhs = -1;
        int rhs = -1;
        double d_lhs = 0.0;
        double d_rhs = 0.0;
    };

    std::vector<Node> nodes_;
    std::vector<double> adjoints_;

    inline void validateOwner(const Reverse &value) const
    {
        if (value.tape != nullptr && value.tape != this)
        {
            throw std::runtime_error("Reverse values belong to different tapes.");
        }
    }

    inline Reverse append(const double value,
                          const int lhs,
                          const int rhs,
                          const double d_lhs,
                          const double d_rhs)
    {
        nodes_.push_back({value, lhs, rhs, d_lhs, d_rhs});
        return Reverse(this, static_cast<int>(nodes_.size()) - 1);
    }

public:
    explicit ReverseTape(const std::size_t reserve_nodes = 0)
    {
        if (reserve_nodes > 0)
        {
            nodes_.reserve(reserve_nodes);
            adjoints_.reserve(reserve_nodes);
        }
    }

    inline void reset()
    {
        nodes_.clear();
        adjoints_.clear();
    }

    inline Reverse variable(const double value = 0.0)
    {
        return append(value, -1, -1, 0.0, 0.0);
    }

    inline double value(const Reverse &input) const
    {
        validateOwner(input);
        return input.node >= 0 ? nodes_[input.node].value : input.constant;
    }

    inline Reverse add(const Reverse &lhs, const Reverse &rhs)
    {
        validateOwner(lhs);
        validateOwner(rhs);
        return append(value(lhs) + value(rhs), lhs.node, rhs.node, 1.0, 1.0);
    }

    inline Reverse subtract(const Reverse &lhs, const Reverse &rhs)
    {
        validateOwner(lhs);
        validateOwner(rhs);
        return append(value(lhs) - value(rhs), lhs.node, rhs.node, 1.0, -1.0);
    }

    inline Reverse multiply(const Reverse &lhs, const Reverse &rhs)
    {
        validateOwner(lhs);
        validateOwner(rhs);
        const double lhs_value = value(lhs);
        const double rhs_value = value(rhs);
        return append(lhs_value * rhs_value, lhs.node, rhs.node,
                      rhs_value, lhs_value);
    }

    inline Reverse divide(const Reverse &lhs, const Reverse &rhs)
    {
        validateOwner(lhs);
        validateOwner(rhs);
        const double lhs_value = value(lhs);
        const double rhs_value = value(rhs);
        return append(lhs_value / rhs_value, lhs.node, rhs.node,
                      1.0 / rhs_value,
                      -lhs_value / (rhs_value * rhs_value));
    }

    inline Reverse negate(const Reverse &input)
    {
        validateOwner(input);
        return append(-value(input), input.node, -1, -1.0, 0.0);
    }

    inline void backward(const Reverse &output)
    {
        validateOwner(output);
        adjoints_.assign(nodes_.size(), 0.0);
        if (output.node < 0)
        {
            return;
        }
        adjoints_[output.node] = 1.0;
        for (int node = static_cast<int>(nodes_.size()) - 1; node >= 0; --node)
        {
            const double adjoint = adjoints_[node];
            const Node &current = nodes_[node];
            if (current.lhs >= 0)
            {
                adjoints_[current.lhs] += adjoint * current.d_lhs;
            }
            if (current.rhs >= 0)
            {
                adjoints_[current.rhs] += adjoint * current.d_rhs;
            }
        }
    }

    inline double gradient(const Reverse &input) const
    {
        validateOwner(input);
        return input.node >= 0 && input.node < static_cast<int>(adjoints_.size())
                   ? adjoints_[input.node]
                   : 0.0;
    }
};

inline double Reverse::value() const
{
    return tape == nullptr ? constant : tape->value(*this);
}

inline Reverse selectTapeAndAdd(const Reverse &lhs, const Reverse &rhs)
{
    if (lhs.tape == nullptr && rhs.tape == nullptr)
    {
        return Reverse(lhs.constant + rhs.constant);
    }
    ReverseTape *tape = lhs.tape != nullptr ? lhs.tape : rhs.tape;
    return tape->add(lhs, rhs);
}

inline Reverse selectTapeAndSubtract(const Reverse &lhs, const Reverse &rhs)
{
    if (lhs.tape == nullptr && rhs.tape == nullptr)
    {
        return Reverse(lhs.constant - rhs.constant);
    }
    ReverseTape *tape = lhs.tape != nullptr ? lhs.tape : rhs.tape;
    return tape->subtract(lhs, rhs);
}

inline Reverse selectTapeAndMultiply(const Reverse &lhs, const Reverse &rhs)
{
    if (lhs.tape == nullptr && rhs.tape == nullptr)
    {
        return Reverse(lhs.constant * rhs.constant);
    }
    ReverseTape *tape = lhs.tape != nullptr ? lhs.tape : rhs.tape;
    return tape->multiply(lhs, rhs);
}

inline Reverse selectTapeAndDivide(const Reverse &lhs, const Reverse &rhs)
{
    if (lhs.tape == nullptr && rhs.tape == nullptr)
    {
        return Reverse(lhs.constant / rhs.constant);
    }
    ReverseTape *tape = lhs.tape != nullptr ? lhs.tape : rhs.tape;
    return tape->divide(lhs, rhs);
}

inline Reverse operator+(const Reverse &lhs, const Reverse &rhs)
{
    return selectTapeAndAdd(lhs, rhs);
}

inline Reverse operator-(const Reverse &lhs, const Reverse &rhs)
{
    return selectTapeAndSubtract(lhs, rhs);
}

inline Reverse operator*(const Reverse &lhs, const Reverse &rhs)
{
    return selectTapeAndMultiply(lhs, rhs);
}

inline Reverse operator/(const Reverse &lhs, const Reverse &rhs)
{
    return selectTapeAndDivide(lhs, rhs);
}

inline Reverse operator-(const Reverse &input)
{
    return input.tape == nullptr ? Reverse(-input.constant) : input.tape->negate(input);
}

inline Reverse &Reverse::operator+=(const Reverse &rhs)
{
    *this = *this + rhs;
    return *this;
}

inline Reverse &Reverse::operator-=(const Reverse &rhs)
{
    *this = *this - rhs;
    return *this;
}

inline Reverse &Reverse::operator*=(const Reverse &rhs)
{
    *this = *this * rhs;
    return *this;
}

inline Reverse &Reverse::operator/=(const Reverse &rhs)
{
    *this = *this / rhs;
    return *this;
}

inline Reverse operator*(const Reverse &lhs, const double rhs)
{
    return lhs * Reverse(rhs);
}

inline Reverse operator*(const double lhs, const Reverse &rhs)
{
    return Reverse(lhs) * rhs;
}

inline Reverse operator/(const Reverse &lhs, const double rhs)
{
    return lhs / Reverse(rhs);
}

inline double primal(const Reverse &input)
{
    return input.value();
}

} // namespace ad
} // namespace nubs

#endif
