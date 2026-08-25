#include "nubs/timing_stencil.hpp"

#include <algorithm>
#include <exception>
#include <iostream>
#include <set>

namespace
{

bool shifted(const int knot, const int degree, const int duration)
{
    return knot >= degree + 1 + duration;
}

bool referenceEnergyAffected(const int physical_span,
                             const int degree,
                             const int duration)
{
    const int span = degree + physical_span;
    bool all_shifted = true;
    bool all_unshifted = true;
    const auto update = [&](const int knot)
    {
        all_shifted = all_shifted && shifted(knot, degree, duration);
        all_unshifted = all_unshifted && !shifted(knot, degree, duration);
    };
    update(span);
    update(span + 1);
    for (int knot = span - degree + 1; knot <= span + degree; ++knot)
    {
        update(knot);
    }
    return !all_shifted && !all_unshifted;
}

std::set<int> referenceConstraintRows(const int degree,
                                      const int order,
                                      const int duration,
                                      const int pieces)
{
    const int control_count = pieces + degree;
    std::set<int> result;
    const int rows = control_count;
    for (int row = 0; row < rows; ++row)
    {
        int span = 0;
        int eval = 0;
        if (row < order)
        {
            span = degree;
            eval = degree;
        }
        else if (row < order + pieces - 1)
        {
            const int waypoint = row - order + 1;
            span = degree + waypoint;
            eval = span;
        }
        else
        {
            span = control_count - 1;
            eval = control_count;
        }
        bool all_shifted = shifted(eval, degree, duration);
        bool all_unshifted = !shifted(eval, degree, duration);
        for (int knot = span - degree + 1; knot <= span + degree; ++knot)
        {
            all_shifted = all_shifted && shifted(knot, degree, duration);
            all_unshifted = all_unshifted && !shifted(knot, degree, duration);
        }
        if (!all_shifted && !all_unshifted)
        {
            result.insert(row);
        }
    }
    return result;
}

} // namespace

int main()
{
    try
    {
        for (const int order : {2, 3, 4})
        {
            const int degree = 2 * order - 1;
            for (int pieces = 1; pieces <= 16; ++pieces)
            {
                for (int duration = 0; duration < pieces; ++duration)
                {
                    const auto range = nubs::timing::affectedPhysicalSpans(
                        degree, duration, pieces);
                    for (int span = 0; span < pieces; ++span)
                    {
                        const bool in_stencil =
                            !range.empty() && span >= range.first && span <= range.last;
                        if (in_stencil != referenceEnergyAffected(
                                              span, degree, duration))
                        {
                            throw std::runtime_error("energy stencil mismatch");
                        }
                    }

                    const auto rows = nubs::timing::affectedConstraintRows(
                        degree, order, duration, pieces);
                    const std::set<int> actual(rows.begin(), rows.end());
                    if (actual != referenceConstraintRows(degree, order,
                                                          duration, pieces))
                    {
                        throw std::runtime_error("constraint stencil mismatch");
                    }
                }
            }
        }
        std::cout << "[exact timing stencils] PASS" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_timing_stencil failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
