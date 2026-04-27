#include "tools/optimization_test_cases.hpp"

#include <exception>
#include <iostream>

int main()
{
    try
    {
        nubs_test::runBoundaryGridCase();
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_boundary_grid failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
