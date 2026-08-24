#include "NUBSTrajectory.hpp"
#include "tools/minco_adapter.hpp"
#include "tools/test_common.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr int Dim = 3;

template <int S>
struct CaseSummary
{
    int order = S;
    int piece_num = 0;
    std::array<double, S + 1> derivative_max{};
    double nubs_energy = 0.0;
    double minco_energy = 0.0;
    double energy_abs_error = 0.0;
    double energy_rel_error = 0.0;
};

template <int S>
CaseSummary<S> evaluateCase(const int piece_num,
                            const int sample_num)
{
    const auto data = nubs_test::makeProblem<Dim>(S, piece_num);

    nubs::NUBSTrajectoryT<Dim, S> nubs;
    Eigen::MatrixXd control_points;
    nubs.generate(data.inner_points, data.head_state, data.tail_state,
                  data.durations, control_points);

    nubs_test::SplineMincoAdapter<Dim> minco;
    minco.generate(data.inner_points, data.head_state, data.tail_state,
                   data.durations, S);

    CaseSummary<S> summary;
    summary.piece_num = piece_num;

    const double total_t = data.durations.sum();
    for (int sample = 0; sample <= sample_num; ++sample)
    {
        double t = total_t * static_cast<double>(sample) /
                   static_cast<double>(sample_num);
        if (sample == sample_num)
        {
            t = nubs_test::previousTime(total_t);
        }

        for (int der = 0; der <= S; ++der)
        {
            const double err =
                (nubs.evaluate(t, der) - minco.evaluate(t, der)).norm();
            summary.derivative_max[der] =
                std::max(summary.derivative_max[der], err);
        }
    }

    summary.nubs_energy = nubs.getEnergy();
    summary.minco_energy = minco.getEnergy();
    summary.energy_abs_error =
        std::abs(summary.nubs_energy - summary.minco_energy);
    summary.energy_rel_error =
        nubs_test::relativeError(summary.nubs_energy, summary.minco_energy);
    return summary;
}

template <int S>
void exportSampleCsv(const std::filesystem::path &path,
                     const int piece_num,
                     const int sample_num)
{
    const auto data = nubs_test::makeProblem<Dim>(S, piece_num);

    nubs::NUBSTrajectoryT<Dim, S> nubs;
    Eigen::MatrixXd control_points;
    nubs.generate(data.inner_points, data.head_state, data.tail_state,
                  data.durations, control_points);

    nubs_test::SplineMincoAdapter<Dim> minco;
    minco.generate(data.inner_points, data.head_state, data.tail_state,
                   data.durations, S);

    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error("failed to open output file: " + path.string());
    }

    out << "s,M,total_time,t,nubs_x,nubs_y,nubs_z,minco_x,minco_y,minco_z";
    for (int der = 0; der <= S; ++der)
    {
        out << ",err_d" << der;
    }
    out << "\n";
    out << std::setprecision(17);

    const double total_t = data.durations.sum();
    for (int sample = 0; sample <= sample_num; ++sample)
    {
        double t = total_t * static_cast<double>(sample) /
                   static_cast<double>(sample_num);
        if (sample == sample_num)
        {
            t = nubs_test::previousTime(total_t);
        }

        const Eigen::Vector3d nubs_pos = nubs.evaluate(t, 0);
        const Eigen::Vector3d minco_pos = minco.evaluate(t, 0);

        out << S << ',' << piece_num << ',' << total_t << ',' << t << ','
            << nubs_pos(0) << ',' << nubs_pos(1) << ',' << nubs_pos(2) << ','
            << minco_pos(0) << ',' << minco_pos(1) << ',' << minco_pos(2);
        for (int der = 0; der <= S; ++der)
        {
            out << ',' << (nubs.evaluate(t, der) - minco.evaluate(t, der)).norm();
        }
        out << "\n";
    }
}

template <int S>
void exportWaypointCsv(const std::filesystem::path &path,
                       const int piece_num)
{
    const auto data = nubs_test::makeProblem<Dim>(S, piece_num);

    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error("failed to open output file: " + path.string());
    }

    out << "s,M,index,t,x,y,z,next_duration\n";
    out << std::setprecision(17);

    double t = 0.0;
    for (int i = 0; i <= piece_num; ++i)
    {
        const double next_duration = i < piece_num ? data.durations(i) : 0.0;
        out << S << ',' << piece_num << ',' << i << ',' << t << ','
            << data.waypoints(i, 0) << ',' << data.waypoints(i, 1) << ','
            << data.waypoints(i, 2) << ',' << next_duration << "\n";
        if (i < piece_num)
        {
            t += data.durations(i);
        }
    }
}

template <int S>
void appendSummaryRows(std::ofstream &out,
                       const std::vector<int> &piece_counts,
                       const int sample_num)
{
    for (const int piece_num : piece_counts)
    {
        const auto summary = evaluateCase<S>(piece_num, sample_num);
        out << S << ',' << piece_num;
        for (int der = 0; der <= 4; ++der)
        {
            if (der <= S)
            {
                out << ',' << summary.derivative_max[der];
            }
            else
            {
                out << ',';
            }
        }
        out << ',' << summary.nubs_energy
            << ',' << summary.minco_energy
            << ',' << summary.energy_abs_error
            << ',' << summary.energy_rel_error << "\n";
    }
}

void exportSummaryCsv(const std::filesystem::path &path)
{
    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error("failed to open output file: " + path.string());
    }
    out << std::setprecision(17);
    out << "s,M,max_position,max_velocity,max_acceleration,max_jerk,max_snap,"
           "nubs_energy,minco_energy,energy_abs_error,energy_rel_error\n";

    const std::vector<int> piece_counts{4, 6, 8, 12, 16};
    constexpr int sample_num = 320;
    appendSummaryRows<3>(out, piece_counts, sample_num);
    appendSummaryRows<4>(out, piece_counts, sample_num);
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        std::filesystem::path output_dir = "docs/data";
        if (argc >= 2)
        {
            output_dir = argv[1];
        }
        std::filesystem::create_directories(output_dir);

        exportSampleCsv<4>(output_dir / "fig2_equivalence_samples_s4.csv", 8, 420);
        exportWaypointCsv<4>(output_dir / "fig2_equivalence_waypoints_s4.csv", 8);
        exportSummaryCsv(output_dir / "fig2_equivalence_summary.csv");

        std::cout << "wrote " << (output_dir / "fig2_equivalence_samples_s4.csv")
                  << "\n";
        std::cout << "wrote " << (output_dir / "fig2_equivalence_waypoints_s4.csv")
                  << "\n";
        std::cout << "wrote " << (output_dir / "fig2_equivalence_summary.csv")
                  << "\n";
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "export_equivalence_figure_data failed: " << e.what()
                  << std::endl;
        return 1;
    }
}
