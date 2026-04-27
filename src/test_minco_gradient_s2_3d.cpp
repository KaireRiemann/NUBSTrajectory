#include "tools/test_cases.hpp"

#include <exception>
#include <iostream>

int main()
{
    try
    {
        nubs_test::runMincoEnergyGradient3DCase<2>();
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_minco_gradient_s2_3d failed: "
                  << e.what() << std::endl;
        return 1;
    }
    return 0;
}
