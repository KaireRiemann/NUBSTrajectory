#include "tools/test_cases.hpp"

int main()
{
    try
    {
        nubs_test::runCenteredGradientCase<2, 2>();
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "test_centered_gradient_s2_d2 failed: " << e.what() << std::endl;
        return 1;
    }
}

