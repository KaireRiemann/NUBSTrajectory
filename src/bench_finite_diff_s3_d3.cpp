#include "tools/optimization_test_cases.hpp"

#include <exception>
#include <iostream>

int main()
{
    try
    {
        nubs_test::runFiniteDiffBenchmarkS3D3();
    }
    catch (const std::exception &e)
    {
        std::cerr << "bench_finite_diff_s3_d3 failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
