#include "tools/optimization_test_cases.hpp"

#include <iostream>

int main()
{
    try
    {
        nubs_test::runUniformTimeCase();
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_uniform_time failed: "
                  << e.what() << std::endl;
        return 1;
    }
    return 0;
}
