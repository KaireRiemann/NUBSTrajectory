#ifndef NUBS_TIMING_STENCIL_HPP
#define NUBS_TIMING_STENCIL_HPP

#include <algorithm>
#include <vector>

namespace nubs
{
namespace timing
{

struct IndexRange
{
    int first = 1;
    int last = 0;

    constexpr bool empty() const { return first > last; }
};

// Physical span r is affected by duration T_k precisely when
// r - p + 1 <= k <= r + p - 1.  The returned range is in physical-span
// coordinates [0, M-1], not global knot indices.
inline IndexRange affectedPhysicalSpans(const int degree,
                                        const int duration_index,
                                        const int piece_num)
{
    return {std::max(0, duration_index - degree + 1),
            std::min(piece_num - 1, duration_index + degree - 1)};
}

// Constraint matrix rows affected by T_k.  Row layout is: s start-state rows,
// M-1 waypoint rows, and s terminal-state rows.
inline std::vector<int> affectedConstraintRows(const int degree,
                                                const int continuity_order,
                                                const int duration_index,
                                                const int piece_num)
{
    std::vector<int> rows;
    if (piece_num <= 0)
    {
        return rows;
    }
    if (duration_index <= degree - 1)
    {
        for (int row = 0; row < continuity_order; ++row)
        {
            rows.push_back(row);
        }
    }

    const int waypoint_begin = std::max(1, duration_index - degree + 1);
    const int waypoint_end = std::min(piece_num - 1,
                                      duration_index + degree - 1);
    for (int waypoint = waypoint_begin; waypoint <= waypoint_end; ++waypoint)
    {
        rows.push_back(continuity_order + waypoint - 1);
    }

    if (duration_index >= piece_num - degree)
    {
        const int tail_start = continuity_order + piece_num - 1;
        for (int row = 0; row < continuity_order; ++row)
        {
            rows.push_back(tail_start + row);
        }
    }
    return rows;
}

} // namespace timing
} // namespace nubs

#endif
