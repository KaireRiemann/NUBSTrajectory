#include "tools/test_cases.hpp"

int main()
{
    try
    {
        nubs_test::runBasicConstructionCase<3, 3>();
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_basic_s3_d3 failed: " << e.what() << std::endl;
        return 1;
    }
}

