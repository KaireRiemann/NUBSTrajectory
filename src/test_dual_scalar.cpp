#include "nubs/basis.hpp"

#include <cmath>
#include <exception>
#include <iostream>

namespace
{

void requireNear(const double actual,
                 const double expected,
                 const char *label,
                 const double tolerance = 1.0e-12)
{
    if (std::abs(actual - expected) > tolerance)
    {
        throw std::runtime_error(std::string(label) + " mismatch");
    }
}

struct ShiftedKnotView
{
    double shift = 0.0;

    double operator()(const int index) const
    {
        static constexpr double knots[] = {0.0, 0.0, 0.0, 0.0,
                                           1.0, 2.0, 2.0, 2.0};
        return knots[index] + (index >= 4 ? shift : 0.0);
    }
};

struct DualKnotView
{
    nubs::ad::Dual operator()(const int index) const
    {
        static constexpr double knots[] = {0.0, 0.0, 0.0, 0.0,
                                           1.0, 2.0, 2.0, 2.0};
        return nubs::ad::Dual(knots[index], index >= 4 ? 1.0 : 0.0);
    }
};

} // namespace

int main()
{
    try
    {
        using nubs::ad::Dual;
        const Dual x(2.0, 1.0);
        const Dual y(3.0, -0.5);
        const Dual expression = (x * y + x / y - y) / x;
        requireNear(expression.value, 1.8333333333333333, "dual value");
        requireNear(expression.derivative, 0.5555555555555556,
                    "dual derivative");

        std::array<std::array<double, 4>, 4> double_ders{};
        std::array<std::array<Dual, 4>, 4> dual_ders{};
        nubs::basis::dersBasisFuns<double, 3>(3, 3, 3, 0.5,
                                               ShiftedKnotView{}, double_ders);
        nubs::basis::dersBasisFuns<Dual, 3>(3, 3, 3, Dual(0.5, 0.0),
                                             DualKnotView{}, dual_ders);
        for (int d = 0; d <= 3; ++d)
        {
            for (int j = 0; j <= 3; ++j)
            {
                requireNear(dual_ders[d][j].value, double_ders[d][j],
                            "scalar-generic basis value");

                constexpr double h = 1.0e-6;
                std::array<std::array<double, 4>, 4> plus{};
                std::array<std::array<double, 4>, 4> minus{};
                nubs::basis::dersBasisFuns<double, 3>(
                    3, 3, 3, 0.5, ShiftedKnotView{h}, plus);
                nubs::basis::dersBasisFuns<double, 3>(
                    3, 3, 3, 0.5, ShiftedKnotView{-h}, minus);
                const double finite_difference =
                    (plus[d][j] - minus[d][j]) / (2.0 * h);
                requireNear(dual_ders[d][j].derivative, finite_difference,
                            "scalar-generic basis derivative", 1.0e-6);
            }
        }
        std::cout << "[dual scalar and scalar-generic basis] PASS" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_dual_scalar failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
