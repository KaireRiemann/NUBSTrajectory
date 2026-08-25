#include "NUBSTrajectory.hpp"
#include "tools/minco_adapter.hpp"
#include "tools/optimization_test_cases.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr int Dim = 3;
constexpr int PieceNum = 6;
constexpr int SamplesPerTrajectory = 240;

template <int S>
class SpanEvaluator final : public nubs::NUBSTrajectoryT<Dim, S>
{
    using GenericBase = nubs::NUBSTrajectory<Dim, 2 * S - 1>;
    static constexpr int P = 2 * S - 1;

public:
    Eigen::Matrix<double, Dim, 1> evaluateOnSpan(const double t,
                                                  const int derivative,
                                                  const int span) const
    {
        Eigen::Matrix<double, P + 1, P + 1> ders;
        GenericBase::dersBasisFuns(derivative, span, t, this->knots, ders);
        Eigen::Matrix<double, Dim, 1> value =
            Eigen::Matrix<double, Dim, 1>::Zero();
        for (int local = 0; local <= P; ++local)
        {
            value += ders(derivative, local) *
                     this->control_points.row(span - P + local).transpose();
        }
        return value;
    }
};

struct RobustnessRow
{
    int order = 0;
    int piece_num = 0;
    int runs = 0;
    double duration_ratio = 0.0;
    double total_duration = 0.0;
    double min_duration = 0.0;
    double max_duration = 0.0;
    double max_position_error_vs_minco = 0.0;
    double energy_abs_error_vs_minco = 0.0;
    double energy_rel_error_vs_minco = 0.0;
    double max_boundary_residual = 0.0;
    double max_waypoint_residual = 0.0;
    double max_waypoint_derivative_error_vs_minco = 0.0;
    double max_continuity_jump = 0.0;
    std::array<double, 5> continuity_jump_by_derivative{};
    double full_system_condition_number_2 = 0.0;
    double reduced_uniform_system_condition_number_2 = 0.0;
    double min_abs_pivot = std::numeric_limits<double>::infinity();
    double max_relative_linear_residual = 0.0;
};

inline Eigen::VectorXd stronglyNonuniformDurations(const double ratio)
{
    Eigen::VectorXd durations(PieceNum);
    for (int piece = 0; piece < PieceNum; ++piece)
    {
        durations(piece) = piece % 2 == 0 ? 1.0 : ratio;
    }
    // Keep the total duration fixed so the experiment isolates allocation
    // non-uniformity rather than a global time scaling effect.
    durations *= static_cast<double>(PieceNum) / durations.sum();
    return durations;
}

inline void updateMaximum(double &target, const double value)
{
    if (!std::isfinite(value))
    {
        throw std::runtime_error("non-finite robustness metric");
    }
    target = std::max(target, value);
}

template <int S>
void evaluateOneCase(const nubs_test::ProblemData<Dim> &data,
                     RobustnessRow &row)
{
    constexpr int P = 2 * S - 1;
    SpanEvaluator<S> trajectory;
    Eigen::MatrixXd controls;
    trajectory.generate(data.inner_points, data.head_state, data.tail_state,
                        data.durations, controls);
    nubs_test::requireFiniteMatrix(controls, "non-finite NUBS controls");

    nubs_test::SplineMincoAdapter<Dim> minco;
    minco.generate(data.inner_points, data.head_state, data.tail_state,
                   data.durations, S);

    updateMaximum(row.energy_abs_error_vs_minco,
                  std::abs(trajectory.getEnergy() - minco.getEnergy()));
    updateMaximum(row.energy_rel_error_vs_minco,
                  nubs_test::relativeError(trajectory.getEnergy(), minco.getEnergy()));

    const double total_duration = data.durations.sum();
    for (int sample = 0; sample <= SamplesPerTrajectory; ++sample)
    {
        double t = total_duration * static_cast<double>(sample) /
                   static_cast<double>(SamplesPerTrajectory);
        if (sample == SamplesPerTrajectory)
        {
            t = nubs_test::previousTime(total_duration);
        }
        updateMaximum(row.max_position_error_vs_minco,
                      (trajectory.evaluate(t) - minco.evaluate(t, 0)).norm());
    }

    const int control_num = trajectory.getCtrlPtNum(PieceNum);
    for (int derivative = 0; derivative < S; ++derivative)
    {
        updateMaximum(row.max_boundary_residual,
                      (trajectory.evaluateOnSpan(0.0, derivative, P) -
                       data.head_state.col(derivative)).norm());
        updateMaximum(row.max_boundary_residual,
                      (trajectory.evaluateOnSpan(total_duration, derivative,
                                                 control_num - 1) -
                       data.tail_state.col(derivative)).norm());
    }

    for (int waypoint = 1; waypoint < PieceNum; ++waypoint)
    {
        const double t = nubs_test::cumulativeTime(data.durations, waypoint);
        const int left_span = P + waypoint - 1;
        const int right_span = P + waypoint;
        const Eigen::Matrix<double, Dim, 1> expected =
            data.inner_points.row(waypoint - 1).transpose();
        updateMaximum(row.max_waypoint_residual,
                      (trajectory.evaluateOnSpan(t, 0, left_span) - expected).norm());
        updateMaximum(row.max_waypoint_residual,
                      (trajectory.evaluateOnSpan(t, 0, right_span) - expected).norm());

        for (int derivative = 1; derivative <= S; ++derivative)
        {
            updateMaximum(row.max_waypoint_derivative_error_vs_minco,
                          (trajectory.evaluateOnSpan(t, derivative, right_span) -
                           minco.evaluate(t, derivative)).norm());
        }
        for (int derivative = 0; derivative <= S; ++derivative)
        {
            const double jump =
                (trajectory.evaluateOnSpan(t, derivative, left_span) -
                 trajectory.evaluateOnSpan(t, derivative, right_span)).norm();
            updateMaximum(row.max_continuity_jump, jump);
            updateMaximum(row.continuity_jump_by_derivative[derivative], jump);
        }
    }

    updateMaximum(row.full_system_condition_number_2,
                  trajectory.getFullSystemConditionNumber());
    row.min_abs_pivot = std::min(
        row.min_abs_pivot, trajectory.A.factorizationStats().minimum_abs_pivot);
    updateMaximum(row.max_relative_linear_residual,
                  trajectory.getLastLinearSolveRelativeResidual());

    nubs::UniformNUBSTrajectoryT<Dim, S> uniform_trajectory;
    Eigen::MatrixXd uniform_controls;
    uniform_trajectory.generateUniform(data.inner_points, data.head_state,
                                       data.tail_state, total_duration,
                                       uniform_controls);
    updateMaximum(row.reduced_uniform_system_condition_number_2,
                  uniform_trajectory.getReducedUniformSystemConditionNumber());
}

template <int S>
RobustnessRow runOneRatio(const double ratio, const int runs)
{
    RobustnessRow row;
    row.order = S;
    row.piece_num = PieceNum;
    row.runs = runs;
    row.duration_ratio = ratio;
    const Eigen::VectorXd durations = stronglyNonuniformDurations(ratio);
    row.total_duration = durations.sum();
    row.min_duration = durations.minCoeff();
    row.max_duration = durations.maxCoeff();

    for (int run = 0; run < runs; ++run)
    {
        const auto data = nubs_test::makeRandomProblem<Dim>(
            S, PieceNum,
            static_cast<unsigned int>(910000 + 10000 * S +
                                      100 * static_cast<int>(ratio) + run));
        auto case_data = data;
        case_data.durations = durations;
        evaluateOneCase<S>(case_data, row);
    }

    std::cout << "S=" << row.order << ", R=" << row.duration_ratio
              << std::scientific << std::setprecision(3)
              << ", pos_vs_minco=" << row.max_position_error_vs_minco
              << ", energy_rel_vs_minco=" << row.energy_rel_error_vs_minco
              << ", boundary=" << row.max_boundary_residual
              << ", waypoint=" << row.max_waypoint_residual
              << ", derivative=" << row.max_waypoint_derivative_error_vs_minco
              << ", continuity=" << row.max_continuity_jump
              << ", cond_full=" << row.full_system_condition_number_2
              << ", cond_reduced=" << row.reduced_uniform_system_condition_number_2
              << ", residual=" << row.max_relative_linear_residual
              << std::endl;
    return row;
}

inline void writeCsv(const std::string &path,
                     const std::vector<RobustnessRow> &rows)
{
    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error("failed to open CSV output: " + path);
    }
    out << "order,piece_count,runs,duration_ratio,total_duration,min_duration,"
           "max_duration,max_position_error_vs_minco,energy_abs_error_vs_minco,"
           "energy_rel_error_vs_minco,max_boundary_residual,"
           "max_waypoint_residual,max_waypoint_derivative_error_vs_minco,"
           "max_continuity_jump,continuity_jump_d0,continuity_jump_d1,"
           "continuity_jump_d2,continuity_jump_d3,continuity_jump_d4,"
           "full_system_condition_number_2,reduced_uniform_system_condition_number_2,"
           "min_abs_pivot,max_relative_linear_residual\n";
    out << std::setprecision(17);
    for (const RobustnessRow &row : rows)
    {
        out << row.order << ',' << row.piece_num << ',' << row.runs << ','
            << row.duration_ratio << ',' << row.total_duration << ','
            << row.min_duration << ',' << row.max_duration << ','
            << row.max_position_error_vs_minco << ','
            << row.energy_abs_error_vs_minco << ','
            << row.energy_rel_error_vs_minco << ','
            << row.max_boundary_residual << ','
            << row.max_waypoint_residual << ','
            << row.max_waypoint_derivative_error_vs_minco << ','
            << row.max_continuity_jump;
        for (const double jump : row.continuity_jump_by_derivative)
        {
            out << ',' << jump;
        }
        out << ',' << row.full_system_condition_number_2 << ','
            << row.reduced_uniform_system_condition_number_2 << ','
            << row.min_abs_pivot << ',' << row.max_relative_linear_residual
            << '\n';
    }
}

} // namespace

int main(int argc, char **argv)
{
    const int runs = argc > 1 ? std::max(1, std::atoi(argv[1])) : 20;
    const std::string csv_path = argc > 2
                                     ? argv[2]
                                     : "nonuniform_robustness.csv";
    try
    {
        const std::vector<double> ratios = {1.0, 10.0, 100.0, 1000.0};
        std::vector<RobustnessRow> rows;
        rows.reserve(3 * ratios.size());
        std::cout << "[strongly non-uniform duration robustness]" << std::endl;
        std::cout << "fixed total duration=" << PieceNum
                  << ", M=" << PieceNum << ", runs=" << runs << std::endl;
        for (const double ratio : ratios)
        {
            rows.push_back(runOneRatio<2>(ratio, runs));
            rows.push_back(runOneRatio<3>(ratio, runs));
            rows.push_back(runOneRatio<4>(ratio, runs));
        }
        writeCsv(csv_path, rows);
        std::cout << "wrote CSV: " << csv_path << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "bench_nonuniform_robustness failed: " << e.what()
                  << std::endl;
        return 1;
    }
    return 0;
}
