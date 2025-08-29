/**
 * A collection of spherical functions
 *
 * Copyright (c) 2025, Michal Vlnas
 */

#include <iostream>
#include <functions.h>

int main()
{
    auto val1 = sphc::cf_f1<float>(0.23f, 0.42f);
    std::cout << val1 << std::endl;

    auto val2 = sphc::fornberg_f1(0.2, 0.1);
    std::cout << val2 << std::endl;

    auto val3 = sphc::beentjes_f4<double>(0.5, 1.0);
    std::cout << val3 << std::endl;

    return 0;
}
