/**
 * A collection of spherical functions
 *
 * Copyright (c) 2025, Michal Vlnas
 */

#ifndef SPHERICAL_COLLECTION_FUNCTIONS_H
#define SPHERICAL_COLLECTION_FUNCTIONS_H

#include <tuple>
#include <cmath>

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

template<typename Float>
Float fornberg_f1(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 1.0f + x + y * y + x * x * y + x * x * x * x + y * y * y * y * y + x * x * y * y * z * z;
}

template<typename Float>
Float fornberg_f4(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f + (Float)sgn(-9.0f * x - 9.0f * y + 9.0f * z)) / 9.0f;
}

template<typename Float>
Float beentjes_f3(Float const theta, Float const phi)
{
    Float constexpr alpha = 9.0f;
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f + std::tanh(-alpha * x - alpha * y + alpha * z)) / alpha;
}

template<typename Float>
Float beentjes_f4(Float const theta, Float const phi)
{
    Float constexpr alpha = 9.0f;
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f - (Float)sgn(x + y - z)) / alpha;
}

template<typename Float>
Float beentjes_f5(Float const theta, Float const phi)
{
    Float constexpr alpha = 9.0f;
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f - (Float)sgn(F_PI * x + y)) / alpha;
}

template<typename Float>
Float renka_f3(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::abs((1.25f + std::cos(5.4f * y)) * std::cos(6.0f * z) / (6.0f + 6.0f * (3.0f * x - 1) * (3.0f * x - 1)));
}

template<typename Float>
Float renka_f4(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    // exp[-(81/16)((X - .5)^2 + (Y - .5)^2 + (Z - .5)^2)]/3;
    return std::exp(-(81.0f / 16.0f) *
        (std::pow(x - 0.5f, 2.0f) + std::pow(y - 0.5f, 2.0f) + std::pow(z - 0.5f, 2.0f))) / 3.0f;
}

template<typename Float>
Float renka_f5(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    // exp[-(81/4)((X - .5)^2 + (Y - .5)^2 + (Z - .5)^2)]/3;
    return std::exp(-(81.0f / 4.0f) *
        (std::pow(x - 0.5f, 2.0f) + std::pow(y - 0.5f, 2.0f) + std::pow(z - 0.5f, 2.0f))) / 3.0f;
}

// From: "Numerical Quadrature over the Surface of a Sphere"
template<typename Float>
Float reegar_f3(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (F_PI_2 + std::atan(300.0f * (z - 9999.0f / 10000.0f))) / F_PI;
}

// From: "Spherical Harmonics Collocation: A Computational Intercomparison of Several Grids"
template<typename Float>
Float bellet_f4(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 0.5f * (1.0f + (Float)sgn(x - 0.5f));
}

// From: "Numerical quadrature over smooth surfaces with boundaries"
template<typename Float>
Float reegar_f2(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 2.0f / F_PI * std::atan(z);
}

// From: "Numerical quadrature over smooth surfaces with boundaries"
template<typename Float>
Float reegar_f4(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 0.5f + std::atan(1000.0f * (z - 9999.0f / (10000.0f * 2.0f * std::sqrt(2)))) / F_PI;
}

template<typename Float>
Float franke(Float const theta, Float const phi)
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

/**
 * Custom designed functions - partially present in Vlnas2025 et al.
 */

template<typename Float>
Float cf_f1(Float const theta, Float const phi)
{
    // (Abs[Sin[Cos[2 \[Phi]] - 2 \[Theta]]] + Abs[Cos[2 \[Theta]]])
    return std::abs(std::sin(std::cos(2.0f * phi) - 2.0f * theta)) + std::abs(std::cos(2.0f * theta));
}

template<typename Float>
Float cf_f2(Float const theta, Float const phi)
{
    // Abs[Sin[2 \[Phi] - \[Theta]]] + Abs[Cos[2 \[Theta]]]
    return std::abs(std::sin(2.0f * phi - theta)) + std::abs(std::cos(2.0f * theta));
}

template<typename Float>
Float cf_f3(Float const /*theta*/, Float const phi)
{
    // 1 + Sin[5 \[Phi]]/5
    return 1.0f + std::sin(5.0f * phi) / 5.0f;
}

template<typename Float>
Float cf_f4(Float const theta, Float const phi)
{
    // 1 + Cos[5 \[Phi]]/5 + Sin[5 \[Theta]]
    return 1.0f + std::cos(5.0f * phi) / 5.0f + std::sin(5.0f * theta);
}

template<typename Float>
Float cf_f5(Float const theta, Float const phi)
{
    return std::exp(2.0f * dot(
                std::sin(theta) * std::cos(phi),
                std::sin(theta) * std::sin(phi),
                std::cos(theta), -1.0f, -1.0f, 0.8f)) +
            std::exp(1.5f * dot(
                std::sin(theta) * std::cos(phi),
                std::sin(theta) * std::sin(phi),
                std::cos(theta), 1.0f, -1.0f, 0.8f)) +
            std::exp(theta) + 10.0f * std::exp(dot(
                    std::sin(theta) * std::cos(phi),
                    std::sin(theta) * std::sin(phi),
                    std::cos(theta), 0.8f, 0.3f, -4.0f) - 1.0f) +
            4.0f * std::abs(std::cos(45.0f * theta + 45.0f * phi));
}

template<typename Float>
Float cf_f6(Float const theta, Float const phi)
{
    return 1.0f + 0.5f * std::cos(theta) + 0.3f * std::cos(2.0f * phi);
}

template<typename Float>
Float cf_f7(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::abs(std::cos(3.0f * x) + std::sin(2.0f * y) + 0.5f * z * z);
}

template<typename Float>
Float cf_f8(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::abs(std::sin(2.0f * x) * std::cos(3.0f * y) + 0.5f * z * z + 0.3f * std::sin(5.0f * x) * std::cos(4.0f * z));
}

template<typename Float>
Float cf_f9(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::abs(x * x - y * y + 0.5f * x * z - 0.3f * y * z);
}

template<typename Float>
Float cf_f10(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return x * x + y * y + z * z + 5 + 2.5f * std::cos((theta - F_PI) / 2.0f) * std::sin(16.0f * theta);
}

template<typename Float>
Float cf_f11(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    // f(x,y,z)=|sin(10x)cos(12y)sin(15z)+cos(20x)|
    return std::abs(std::sin(10.0f * x) * std::cos(12.0f * y) * std::sin(15.0f * z) + std::cos(20.0f * x));
}

template<typename Float>
Float cf_f12(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    // f(x,y,z)=sin(10x)+cos(12y)−sin(15z)+0.2cos(18x) + 3
    return std::sin(10.0f * x) + std::cos(12.0f * y) - std::sin(15.0f * z) + 0.2f * std::cos(18.0f * x) + 3.0f;
}

template<typename Float>
Float cf_f13(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::exp(-std::sin(5.0f * x) - std::cos(6.0f * y)) + 0.3f * std::sin(10.0f * z);
}

template<typename Float>
Float cf_f14(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::exp(-2.0f * (x * x + y * y)) * std::sin(4.0f * z);
}

template<typename Float>
Float cf_15(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (x * x + y * y) * std::exp(-3.0f * z * z);
}

} // namespace sphc

#endif // SPHERICAL_COLLECTION_FUNCTIONS_H
