#include "NUBSTrajectory.hpp"
#include "tools/minco_adapter.hpp"
#include "tools/optimization_test_cases.hpp"

#include <algorithm>
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

struct ErrorAccumulator
{
    double max_abs = 0.0;
    double max_rel = 0.0;
    double sum_sq = 0.0;
    std::size_t count = 0;

    void add(const double analytic, const double numeric)
    {
        const double abs_err = std::abs(analytic - numeric);
        const double scale = std::max({1.0, std::abs(analytic),
                                       std::abs(numeric)});
        max_abs = std::max(max_abs, abs_err);
        max_rel = std::max(max_rel, abs_err / scale);
        sum_sq += abs_err * abs_err;
        ++count;
    }

    template <typename DerivedA, typename DerivedB>
    void addMatrix(const Eigen::MatrixBase<DerivedA> &analytic,
                   const Eigen::MatrixBase<DerivedB> &numeric)
    {
        if (analytic.rows() != numeric.rows() ||
            analytic.cols() != numeric.cols())
        {
            throw std::runtime_error("gradient shape mismatch");
        }
        for (Eigen::Index r = 0; r < analytic.rows(); ++r)
        {
            for (Eigen::Index c = 0; c < analytic.cols(); ++c)
            {
                add(analytic(r, c), numeric(r, c));
            }
        }
    }

    double rms() const
    {
        if (count == 0)
        {
            return 0.0;
        }
        return std::sqrt(sum_sq / static_cast<double>(count));
    }
};

struct GradientRow
{
    int order = 0;
    int piece_num = 0;
    int runs = 0;
    int point_vars = 0;
    int time_vars = 0;
    int total_vars = 0;
    double mean_energy = 0.0;
    double external_point_max_abs_error = 0.0;
    double external_point_rms_abs_error = 0.0;
    double external_point_max_rel_error = 0.0;
    double external_time_max_abs_error = 0.0;
    double external_time_rms_abs_error = 0.0;
    double external_time_max_rel_error = 0.0;
    double minco_energy_max_abs_error = 0.0;
    double minco_energy_max_rel_error = 0.0;
    double minco_point_max_abs_error = 0.0;
    double minco_point_rms_abs_error = 0.0;
    double minco_point_max_rel_error = 0.0;
    double minco_time_max_abs_error = 0.0;
    double minco_time_rms_abs_error = 0.0;
    double minco_time_max_rel_error = 0.0;
};

template <int Dim, int S>
double energyForPerturbedPoint(
    const nubs_test::ProblemData<Dim> &data,
    const int row,
    const int col,
    const double delta)
{
    auto points = data.inner_points;
    points(row, col) += delta;
    return nubs_test::energyForFixed<Dim, S>(
        points, data.head_state, data.tail_state, data.durations);
}

template <int Dim, int S>
double energyForPerturbedTime(
    const nubs_test::ProblemData<Dim> &data,
    const int index,
    const double delta)
{
    auto durations = data.durations;
    durations(index) += delta;
    if (!(durations(index) > 1.0e-8))
    {
        durations(index) = 1.0e-8;
    }
    return nubs_test::energyForFixed<Dim, S>(
        data.inner_points, data.head_state, data.tail_state, durations);
}

template <int Dim, int S>
void computeExternalFiniteDiffGradient(
    const nubs_test::ProblemData<Dim> &data,
    Eigen::Matrix<double, Eigen::Dynamic, Dim> &grad_points,
    Eigen::VectorXd &grad_times)
{
    const int M = static_cast<int>(data.durations.size());
    grad_points.setZero(std::max(0, M - 1), Dim);
    grad_times.setZero(M);

    for (int r = 0; r < data.inner_points.rows(); ++r)
    {
        for (int c = 0; c < Dim; ++c)
        {
            const double h = 1.0e-4 *
                             std::max(1.0, std::abs(data.inner_points(r, c)));
            const double plus_2_cost =
                energyForPerturbedPoint<Dim, S>(data, r, c, 2.0 * h);
            const double plus_1_cost =
                energyForPerturbedPoint<Dim, S>(data, r, c, h);
            const double minus_1_cost =
                energyForPerturbedPoint<Dim, S>(data, r, c, -h);
            const double minus_2_cost =
                energyForPerturbedPoint<Dim, S>(data, r, c, -2.0 * h);
            grad_points(r, c) =
                (-plus_2_cost + 8.0 * plus_1_cost -
                 8.0 * minus_1_cost + minus_2_cost) /
                (12.0 * h);
        }
    }

    for (int i = 0; i < M; ++i)
    {
        double h = 1.0e-4 * std::max(1.0, std::abs(data.durations(i)));
        h = std::min(h, 0.20 * data.durations(i));
        const double plus_2_cost =
            energyForPerturbedTime<Dim, S>(data, i, 2.0 * h);
        const double plus_1_cost =
            energyForPerturbedTime<Dim, S>(data, i, h);
        const double minus_1_cost =
            energyForPerturbedTime<Dim, S>(data, i, -h);
        const double minus_2_cost =
            energyForPerturbedTime<Dim, S>(data, i, -2.0 * h);
        grad_times(i) =
            (-plus_2_cost + 8.0 * plus_1_cost -
             8.0 * minus_1_cost + minus_2_cost) /
            (12.0 * h);
    }
}

template <int S>
GradientRow runOne(const int piece_num, const int runs)
{
    constexpr int Dim = 3;

    ErrorAccumulator point_errors;
    ErrorAccumulator time_errors;
    ErrorAccumulator minco_energy_errors;
    ErrorAccumulator minco_point_errors;
    ErrorAccumulator minco_time_errors;

    double energy_sum = 0.0;

    for (int run = 0; run < runs; ++run)
    {
        const unsigned int seed =
            static_cast<unsigned int>(710000 + 10000 * S +
                                      100 * piece_num + run);
        const auto data =
            nubs_test::makeRandomProblem<Dim>(S, piece_num, seed);

        nubs::NUBSTrajectoryT<Dim, S> trajectory;
        Eigen::MatrixXd control_points;
        trajectory.generate(data.inner_points, data.head_state,
                            data.tail_state, data.durations, control_points);

        double nubs_cost = 0.0;
        Eigen::MatrixXd nubs_grad_points;
        Eigen::VectorXd nubs_grad_times;
        trajectory.getEnergyAndFiniteDiffGrad(nubs_cost, nubs_grad_points,
                                              nubs_grad_times);
        energy_sum += nubs_cost;

        Eigen::Matrix<double, Eigen::Dynamic, Dim> numeric_grad_points;
        Eigen::VectorXd numeric_grad_times;
        computeExternalFiniteDiffGradient<Dim, S>(
            data, numeric_grad_points, numeric_grad_times);

        point_errors.addMatrix(nubs_grad_points, numeric_grad_points);
        time_errors.addMatrix(nubs_grad_times, numeric_grad_times);

        nubs_test::SplineMincoAdapter<Dim> minco;
        minco.generate(data.inner_points, data.head_state, data.tail_state,
                       data.durations, S);
        double minco_cost = 0.0;
        Eigen::MatrixXd minco_coeff_grad;
        Eigen::Matrix<double, Eigen::Dynamic, Dim> minco_grad_points;
        Eigen::VectorXd minco_grad_times;
        Eigen::VectorXd minco_partial_grad_times;
        minco.getEnergyAndGrad(minco_cost, minco_coeff_grad,
                               minco_grad_points, minco_grad_times,
                               minco_partial_grad_times);

        minco_energy_errors.add(nubs_cost, minco_cost);
        minco_point_errors.addMatrix(nubs_grad_points, minco_grad_points);
        minco_time_errors.addMatrix(nubs_grad_times, minco_grad_times);
    }

    GradientRow row;
    row.order = S;
    row.piece_num = piece_num;
    row.runs = runs;
    row.point_vars = std::max(0, piece_num - 1) * Dim;
    row.time_vars = piece_num;
    row.total_vars = row.point_vars + row.time_vars;
    row.mean_energy = energy_sum / static_cast<double>(runs);
    row.external_point_max_abs_error = point_errors.max_abs;
    row.external_point_rms_abs_error = point_errors.rms();
    row.external_point_max_rel_error = point_errors.max_rel;
    row.external_time_max_abs_error = time_errors.max_abs;
    row.external_time_rms_abs_error = time_errors.rms();
    row.external_time_max_rel_error = time_errors.max_rel;
    row.minco_energy_max_abs_error = minco_energy_errors.max_abs;
    row.minco_energy_max_rel_error = minco_energy_errors.max_rel;
    row.minco_point_max_abs_error = minco_point_errors.max_abs;
    row.minco_point_rms_abs_error = minco_point_errors.rms();
    row.minco_point_max_rel_error = minco_point_errors.max_rel;
    row.minco_time_max_abs_error = minco_time_errors.max_abs;
    row.minco_time_rms_abs_error = minco_time_errors.rms();
    row.minco_time_max_rel_error = minco_time_errors.max_rel;

    std::cout << "S=" << S
              << ", M=" << piece_num
              << ", runs=" << runs
              << std::scientific << std::setprecision(3)
              << ", external_point_max="
              << row.external_point_max_abs_error
              << ", external_time_max="
              << row.external_time_max_abs_error
              << ", minco_energy_max="
              << row.minco_energy_max_abs_error
              << ", minco_point_max=" << row.minco_point_max_abs_error
              << ", minco_time_max=" << row.minco_time_max_abs_error
              << std::endl;
    return row;
}

std::vector<GradientRow> runBenchmark(const int runs)
{
    const std::vector<int> piece_counts = {2, 4, 8, 16, 32, 64};
    std::vector<GradientRow> rows;
    rows.reserve(2 * piece_counts.size());

    std::cout << "[energy gradient validation]" << std::endl;
    std::cout << "orders: s=3, s=4" << std::endl;
    std::cout << "piece counts: 2, 4, 8, 16, 32, 64" << std::endl;
    std::cout << "runs per row: " << runs << std::endl;

    for (const int piece_num : piece_counts)
    {
        rows.push_back(runOne<3>(piece_num, runs));
        rows.push_back(runOne<4>(piece_num, runs));
    }
    return rows;
}

void writeCsv(const std::string &path,
              const std::vector<GradientRow> &rows)
{
    std::ofstream out(path);
    if (!out)
    {
        throw std::runtime_error("failed to open CSV output: " + path);
    }

    out << "order,piece_count,runs,point_vars,time_vars,total_vars,"
        << "mean_energy,"
        << "external_point_max_abs_error,external_point_rms_abs_error,"
        << "external_point_max_rel_error,"
        << "external_time_max_abs_error,external_time_rms_abs_error,"
        << "external_time_max_rel_error,"
        << "minco_energy_max_abs_error,minco_energy_max_rel_error,"
        << "minco_point_max_abs_error,minco_point_rms_abs_error,"
        << "minco_point_max_rel_error,"
        << "minco_time_max_abs_error,minco_time_rms_abs_error,"
        << "minco_time_max_rel_error\n";
    out << std::setprecision(17);
    for (const auto &row : rows)
    {
        out << row.order << ','
            << row.piece_num << ','
            << row.runs << ','
            << row.point_vars << ','
            << row.time_vars << ','
            << row.total_vars << ','
            << row.mean_energy << ','
            << row.external_point_max_abs_error << ','
            << row.external_point_rms_abs_error << ','
            << row.external_point_max_rel_error << ','
            << row.external_time_max_abs_error << ','
            << row.external_time_rms_abs_error << ','
            << row.external_time_max_rel_error << ','
            << row.minco_energy_max_abs_error << ','
            << row.minco_energy_max_rel_error << ','
            << row.minco_point_max_abs_error << ','
            << row.minco_point_rms_abs_error << ','
            << row.minco_point_max_rel_error << ','
            << row.minco_time_max_abs_error << ','
            << row.minco_time_rms_abs_error << ','
            << row.minco_time_max_rel_error << '\n';
    }
}

} // namespace

int main(int argc, char **argv)
{
    int runs = 100;
    std::string csv_path;
    if (argc > 1)
    {
        runs = std::max(1, std::atoi(argv[1]));
    }
    if (argc > 2)
    {
        csv_path = argv[2];
    }

    try
    {
        const auto rows = runBenchmark(runs);
        if (!csv_path.empty())
        {
            writeCsv(csv_path, rows);
            std::cout << "wrote CSV: " << csv_path << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "bench_energy_gradient_validation failed: "
                  << e.what() << std::endl;
        return 1;
    }
    return 0;
}
