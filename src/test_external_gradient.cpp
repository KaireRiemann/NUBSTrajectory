#include "tools/optimization_test_cases.hpp"

#include <exception>
#include <iostream>

int main()
{
    try
    {
        nubs_test::runExternalGradientCase();
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_external_gradient failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
