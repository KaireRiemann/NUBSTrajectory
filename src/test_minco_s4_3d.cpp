#include "tools/test_cases.hpp"

int main()
{
    try
    {
        nubs_test::runMincoComparison3DCase<4>();
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_minco_s4_3d failed: " << e.what() << std::endl;
        return 1;
    }
}

