#include "tools/optimization_test_cases.hpp"

#include <exception>
#include <iostream>

int main()
{
    try
    {
        nubs_test::runGenericSpecializedEquivalenceCase();
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_generic_specialized_equivalence failed: "
                  << e.what() << std::endl;
        return 1;
    }
    return 0;
}
