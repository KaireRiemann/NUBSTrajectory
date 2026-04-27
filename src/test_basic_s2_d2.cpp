#include "tools/test_cases.hpp"

int main()
{
    try
    {
        nubs_test::runBasicConstructionCase<2, 2>();
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_basic_s2_d2 failed: " << e.what() << std::endl;
        return 1;
    }
}

