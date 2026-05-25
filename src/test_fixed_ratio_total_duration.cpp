#include "tools/optimization_test_cases.hpp"

#include <iostream>

int main()
{
    try
    {
        nubs_test::runFixedRatioTotalDurationCase();
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_fixed_ratio_total_duration failed: "
                  << e.what() << std::endl;
        return 1;
    }
    return 0;
}
