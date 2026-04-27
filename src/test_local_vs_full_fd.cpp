#include "tools/optimization_test_cases.hpp"

#include <exception>
#include <iostream>

int main()
{
    try
    {
        nubs_test::runLocalVsFullFiniteDiffCase();
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_local_vs_full_fd failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
