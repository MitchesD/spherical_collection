/**
 * A collection of spherical functions
 *
 * Copyright (c) 2025-2026, Michal Vlnas
 */

#ifndef SPHERICAL_COLLECTION_FUNCTIONS_H
#define SPHERICAL_COLLECTION_FUNCTIONS_H

#include <tuple>
#include <cmath>
#include <functional>
#include <unordered_map>

namespace sphc
{

#define F_PI            3.14159265358979323846f
#define F_PI_2          1.570796327f
#define F_INV_PI        0.31830988618379067154f
#define F_INV_TWOPI     0.15915494309189533577f

/**
 * Internal helper functions
 */
namespace
{

template<typename Float>
std::tuple<Float, Float, Float> spherical_to_xyz(Float const& theta, Float const& phi)
{
    Float x = std::sin(theta) * std::cos(phi);
    Float y = std::sin(theta) * std::sin(phi);
    Float z = std::cos(theta);
    return { x, y, z };
}

template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

template<typename Float>
Float dot(Float x1, Float y1, Float z1, Float x2, Float y2, Float z2)
{
    return x1 * x2 + y1 * y2 + z1 * z2;
}

}

/**
 * A collection of spherical functions
 */

namespace polynomial
{

// From: "On spherical harmonics based numerical quadrature over the surface of a sphere"
// fornberg_f1
template<typename Float>
Float p1(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 1.0f + x + y * y + x * x * y + x * x * x * x + y * y * y * y * y + x * x * y * y * z * z;
}

} // namespace polynomial

namespace discontinuous
{

// From: "On spherical harmonics based numerical quadrature over the surface of a sphere"
// fornberg_f4
template<typename Float>
Float d1(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f + (Float)sgn(-9.0f * x - 9.0f * y + 9.0f * z)) / 9.0f;
}

// beentjes_f4
template<typename Float>
Float d2(Float const theta, Float const phi)
{
    Float constexpr alpha = 9.0f;
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f - (Float)sgn(x + y - z)) / alpha;
}

// beentjes_f5
template<typename Float>
Float d3(Float const theta, Float const phi)
{
    Float constexpr alpha = 9.0f;
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f - (Float)sgn(F_PI * x + y)) / alpha;
}

// 10. From: "Spherical Harmonics Collocation: A Computational Intercomparison of Several Grids"
// bellet_f4
template<typename Float>
Float d4(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 0.5f * (1.0f + (Float)sgn(x - 0.5f));
}

} // namespace discontinuous

namespace smooth_approx
{

// beentjes_f3
template<typename Float>
Float s1(Float const theta, Float const phi)
{
    Float constexpr alpha = 9.0f;
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f + std::tanh(-alpha * x - alpha * y + alpha * z)) / alpha;
}

// 9. From: "Numerical Quadrature over the Surface of a Sphere"
// reegar_f3
template<typename Float>
Float s2(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (F_PI_2 + std::atan(300.0f * (z - 9999.0f / 10000.0f))) / F_PI;
}

// 12. From: "Numerical quadrature over smooth surfaces with boundaries"
// reegar_f4
template<typename Float>
Float s3(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 0.5f + std::atan(1000.0f * (z - 9999.0f / (10000.0f * 2.0f * std::sqrt(2)))) / F_PI;
}

} // namespace smooth_approx

namespace oscillatory
{

// 6. renka_f3
template<typename Float>
Float o1(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.25f + std::cos(5.4f * y)) * std::cos(6.0f * z) / (6.0f + 6.0f * (3.0f * x - 1) * (3.0f * x - 1));
}

// 23. cf_f10
template<typename Float>
Float o2(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return x * x + y * y + z * z + 5 + 2.5f * std::cos((std::acos(z) - F_PI) / 2.0f) * std::sin(16.0f * std::acos(z));
}

// 25. cf_f12
template<typename Float>
Float o3(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::sin(10.0f * x) + std::cos(12.0f * y) - std::sin(15.0f * z) + 0.2f * std::cos(18.0f * x) + 3.0f;
}

// 26. cf_f13
template<typename Float>
Float o4(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::exp(-std::sin(5.0f * x) - std::cos(6.0f * y)) + 0.3f * std::sin(10.0f * z);
}

// 17. cf_f4
template<typename Float>
Float o5(Float const theta, Float const phi)
{
    return 1.0f + std::cos(5.0f * phi) / 5.0f + std::sin(5.0f * theta);
}

// 18. cf_f5
template<typename Float>
Float o6(Float const theta, Float const phi)
{
    return 1.5f * std::exp(0.8f * (std::sin(theta) * std::cos(phi) + 0.5f * std::cos(theta))) +
           1.2f * std::exp(0.6f * (-std::sin(theta) * std::sin(phi) + 0.3f * std::cos(theta))) +
           0.8f * std::exp(0.5f * std::cos(theta)) + 0.5f * (1.0f + std::cos(6.0f * theta) * std::sin(4.0f * phi));
}

// 19. cf_f6
template<typename Float>
Float o7(Float const theta, Float const phi)
{
    return 1.0f + 0.5f * std::cos(theta) + 0.3f * std::cos(2.0f * phi);
}

} // namespace oscillatory

namespace lobes
{

// 7. renka_f4
template<typename Float>
Float l1(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::exp(-(81.0f / 16.0f) *
        (std::pow(x - 0.5f, 2.0f) + std::pow(y - 0.5f, 2.0f) + std::pow(z - 0.5f, 2.0f))) / 3.0f;
}

// 8. renka_f5
template<typename Float>
Float l2(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::exp(-(81.0f / 4.0f) *
        (std::pow(x - 0.5f, 2.0f) + std::pow(y - 0.5f, 2.0f) + std::pow(z - 0.5f, 2.0f))) / 3.0f;
}

// 13. From: "On spherical harmonics based numerical quadrature over the surface of a sphere"
// fornberg_f2
template<typename Float>
Float l3(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 0.75f * std::exp(-(9.0f * x - 2.0f) * (9.0f * x - 2.0f) / 4.0f -
                (9.0f * y - 2.0f) * (9.0f * y - 2.0f) / 4.0f -
                (9.0f * z - 2.0f) * (9.0f * z - 2.0f) / 4.0f) +
            0.75f * std::exp(-(9.0f * x + 1.0f) * (9.0f * x + 1.0f) / 49.0f -
                (9.0f * y + 1.0f) / 10.0f -
                (9.0f * z + 1.0f) / 10.0f) +
            0.5f * std::exp(-(9.0f * x - 7.0f) * (9.0f * x - 7.0f) / 4.0f -
                (9.0f * y - 3.0f) * (9.0f * y - 3.0f) / 4.0f -
                (9.0f * z - 5.0f) * (9.0f * z - 5.0f) / 4.0f) -
            0.2f * std::exp(-(9.0f * x - 4.0f) * (9.0f * x - 4.0f) -
                (9.0f * y - 7.0f) * (9.0f * y - 7.0f) -
                (9.0f * z - 5.0f) * (9.0f * z - 5.0f));
}

} // namespace lobes

namespace absolute_values
{

// 14. cf_f1
template<typename Float>
Float a1(Float const theta, Float const phi)
{
    return std::abs(std::sin(std::cos(2.0f * phi) - 2.0f * theta)) + std::abs(std::cos(2.0f * theta));
}

// 15. cf_f2
template<typename Float>
Float a2(Float const theta, Float const phi)
{
    return std::abs(std::sin(2.0f * phi - theta)) + std::abs(std::cos(2.0f * theta));
}

// 20. cf_f7
template<typename Float>
Float a3(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::abs(std::cos(3.0f * x) + std::sin(2.0f * y) + 0.5f * z * z);
}

// 21. cf_f8
template<typename Float>
Float a4(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::abs(std::sin(2.0f * x) * std::cos(3.0f * y) + 0.5f * z * z + 0.3f * std::sin(5.0f * x) * std::cos(4.0f * z));
}

// 22. cf_f9
template<typename Float>
Float a5(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::abs(x * x - y * y + 0.5f * x * z - 0.3f * y * z);
}

// 24. cf_f11
template<typename Float>
Float a6(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::abs(std::sin(10.0f * x) * std::cos(12.0f * y) * std::sin(15.0f * z) + std::cos(20.0f * x));
}

} // namespace absolute_values

namespace zsymnetric
{

// 16. cf_f3
template<typename Float>
Float z1(Float const /*theta*/, Float const phi)
{
    return 1.0f + std::sin(5.0f * phi) / 5.0f;
}

// 27. cf_f14
template<typename Float>
Float z2(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::exp(-2.0f * (x * x + y * y)) * std::sin(4.0f * z);
}

// 28. cf_15
template<typename Float>
Float z3(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (x * x + y * y) * std::exp(-3.0f * z * z);
}

} // namespace zsymmetric

/**
 * Get a function by its identifier
 *
 * @param id The identifier of the function (e.g., "p1", "d1", "s1", etc.)
 * @return A std::function that takes two Float arguments (theta and phi) and returns a Float
 */
template<typename Float>
std::function<Float(Float, Float)> get_function(std::string const& id)
{
    static std::unordered_map<std::string, std::function<Float(Float, Float)>> functions =
    {
        { "p1", &sphc::polynomial::p1<Float> },
        { "d1", &sphc::discontinuous::d1<Float> },
        { "d2", &sphc::discontinuous::d2<Float> },
        { "d3", &sphc::discontinuous::d3<Float> },
        { "d4", &sphc::discontinuous::d4<Float> },
        { "s1", &sphc::smooth_approx::s1<Float> },
        { "s2", &sphc::smooth_approx::s2<Float> },
        { "s3", &sphc::smooth_approx::s3<Float> },
        { "o1", &sphc::oscillatory::o1<Float> },
        { "o2", &sphc::oscillatory::o2<Float> },
        { "o3", &sphc::oscillatory::o3<Float> },
        { "o4", &sphc::oscillatory::o4<Float> },
        { "o5", &sphc::oscillatory::o5<Float> },
        { "o6", &sphc::oscillatory::o6<Float> },
        { "o7", &sphc::oscillatory::o7<Float> },
        { "l1", &sphc::lobes::l1<Float> },
        { "l2", &sphc::lobes::l2<Float> },
        { "l3", &sphc::lobes::l3<Float> },
        { "a1", &sphc::absolute_values::a1<Float> },
        { "a2", &sphc::absolute_values::a2<Float> },
        { "a3", &sphc::absolute_values::a3<Float> },
        { "a4", &sphc::absolute_values::a4<Float> },
        { "a5", &sphc::absolute_values::a5<Float> },
        { "a6", &sphc::absolute_values::a6<Float> },
        { "z1", &sphc::zsymnetric::z1<Float> },
        { "z2", &sphc::zsymnetric::z2<Float> },
        { "z3", &sphc::zsymnetric::z3<Float> }
    };

    return functions.at(id);
}

template<typename Float>
Float eval_function(std::string const& id, Float theta, Float phi)
{
    return (get_function<Float>(id))(theta, phi);
}

/**
 * Get function integral by its identifier
 * @param id The identifier of the function (e.g., "p1", "d1", "s1", etc.)
 * @return Surface integral value
 */
template<typename Float>
Float get_integral(std::string const& id)
{
    static std::unordered_map<std::string, Float> integrals =
    {
        { "p1", 19.3881 },
        { "d1", 1.3962 },
        { "d2", 4.0 * M_PI / 9.0 },
        { "d3", 4.0 * M_PI / 9.0 },
        { "d4", M_PI },
        { "s1", 4.0 * M_PI / 9.0 },
        { "s2", 0.0496 },
        { "s3", 4.0634 },
        { "o1", 0.1292 },
        { "o2", 75.3944 },
        { "o3", 37.0324 },
        { "o4", 20.79 },
        { "o5", 12.5664 },
        { "o6", 54.3111 },
        { "o7", 12.5664 },
        { "l1", 0.2181 },
        { "l2", 0.0415 },
        { "l3", 6.6961 },
        { "a1", 15.7323 },
        { "a2", 15.6589 },
        { "a3", 11.5484 },
        { "a4", 5.7017 },
        { "a5", 5.5630 },
        { "a6", 8.5842 },
        { "z1", 12.5664 },
        { "z2", 0.0 },
        { "z3", 5.3857 }
    };

    return static_cast<Float>(integrals.at(id));
}

/**
 * Get function maximum by its identifier
 * @param id The identifier of the function (e.g., "p1", "d1", "s1", etc.)
 * @return Global maximum value
 */
template<typename Float>
Float get_maximum(std::string const& id)
{
    static std::unordered_map<std::string, Float> maximums =
    {
        { "p1", 3.1476 },
        { "d1", 0.2222 },
        { "d2", 0.2222 },
        { "d3", 0.2222 },
        { "d4", 1.0 },
        { "s1", 0.2222 },
        { "s2", 0.5095 },
        { "s3", 0.9995 },
        { "o1", 0.3168 },
        { "o2", 8.4730 },
        { "o3", 5.9997 },
        { "o4", 7.6885 },
        { "o5", 2.2 },
        { "o6", 6.9121 },
        { "o7", 1.8 },
        { "l1", 0.3043 },
        { "l2", 0.2317 },
        { "l3", 2.1802 },
        { "a1", 1.9191 },
        { "a2", 2 },
        { "a3", 2.2536 },
        { "a4", 1.3673 },
        { "a5", 1.0596 },
        { "a6", 1.9963 },
        { "z1", 1.2 },
        { "z2", 0.7568 },
        { "z3", 1 }
    };

    return static_cast<Float>(maximums.at(id));
}

} // namespace sphc

#endif // SPHERICAL_COLLECTION_FUNCTIONS_H
