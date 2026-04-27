#include "tools/test_cases.hpp"

int main()
{
    try
    {
        nubs_test::runBasicConstructionCase<3, 4>();
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_basic_s4_d3 failed: " << e.what() << std::endl;
        return 1;
    }
}

