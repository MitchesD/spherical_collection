/**
 * A collection of spherical functions
 *
 * Copyright (c) 2025-2026, Michal Vlnas
 */

#include <iostream>
#include <functions.h>

int main()
{
    auto const val1 = sphc::zsymnetric::z1<float>(0.23f, 0.42f);
    std::cout << val1 << std::endl;

    auto const val2 = sphc::absolute_values::a1<float>(0.2, 0.1);
    std::cout << val2 << std::endl;

    auto const val3 = sphc::oscillatory::o5<double>(0.5, 1.0);
    std::cout << val3 << std::endl;

    auto const val4 = sphc::eval_function("l1", 0.123, 0.345);
    std::cout << val4 << std::endl;

    return 0;
}
