#include "tools/test_cases.hpp"

#include <exception>
#include <iostream>

int main()
{
    try
    {
        nubs_test::runMincoEnergyGradient3DCase<4>();
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_minco_gradient_s4_3d failed: "
                  << e.what() << std::endl;
        return 1;
    }
    return 0;
}
