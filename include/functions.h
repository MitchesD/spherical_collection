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

// 1. From: "On spherical harmonics based numerical quadrature over the surface of a sphere"
template<typename Float>
Float fornberg_f1(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 1.0f + x + y * y + x * x * y + x * x * x * x + y * y * y * y * y + x * x * y * y * z * z;
}

// 2. From: "On spherical harmonics based numerical quadrature over the surface of a sphere"
template<typename Float>
Float fornberg_f4(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f + (Float)sgn(-9.0f * x - 9.0f * y + 9.0f * z)) / 9.0f;
}

// 3.
template<typename Float>
Float beentjes_f3(Float const theta, Float const phi)
{
    Float constexpr alpha = 9.0f;
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f + std::tanh(-alpha * x - alpha * y + alpha * z)) / alpha;
}

// 4.
template<typename Float>
Float beentjes_f4(Float const theta, Float const phi)
{
    Float constexpr alpha = 9.0f;
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f - (Float)sgn(x + y - z)) / alpha;
}

// 5.
template<typename Float>
Float beentjes_f5(Float const theta, Float const phi)
{
    Float constexpr alpha = 9.0f;
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.0f - (Float)sgn(F_PI * x + y)) / alpha;
}

// 6.
template<typename Float>
Float renka_f3(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (1.25f + std::cos(5.4f * y)) * std::cos(6.0f * z) / (6.0f + 6.0f * (3.0f * x - 1) * (3.0f * x - 1));
}

// 7.
template<typename Float>
Float renka_f4(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    // exp[-(81/16)((X - .5)^2 + (Y - .5)^2 + (Z - .5)^2)]/3;
    return std::exp(-(81.0f / 16.0f) *
        (std::pow(x - 0.5f, 2.0f) + std::pow(y - 0.5f, 2.0f) + std::pow(z - 0.5f, 2.0f))) / 3.0f;
}

// 8.
template<typename Float>
Float renka_f5(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    // exp[-(81/4)((X - .5)^2 + (Y - .5)^2 + (Z - .5)^2)]/3;
    return std::exp(-(81.0f / 4.0f) *
        (std::pow(x - 0.5f, 2.0f) + std::pow(y - 0.5f, 2.0f) + std::pow(z - 0.5f, 2.0f))) / 3.0f;
}

// 9. From: "Numerical Quadrature over the Surface of a Sphere"
template<typename Float>
Float reegar_f3(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (F_PI_2 + std::atan(300.0f * (z - 9999.0f / 10000.0f))) / F_PI;
}

// 10. From: "Spherical Harmonics Collocation: A Computational Intercomparison of Several Grids"
template<typename Float>
Float bellet_f4(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 0.5f * (1.0f + (Float)sgn(x - 0.5f));
}

// 11. From: "Numerical quadrature over smooth surfaces with boundaries"
template<typename Float>
Float reegar_f2(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 2.0f / F_PI * std::atan(z);
}

// 12. From: "Numerical quadrature over smooth surfaces with boundaries"
template<typename Float>
Float reegar_f4(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return 0.5f + std::atan(1000.0f * (z - 9999.0f / (10000.0f * 2.0f * std::sqrt(2)))) / F_PI;
}

// 13. From: "On spherical harmonics based numerical quadrature over the surface of a sphere"
template<typename Float>
Float fornberg_f2(Float const theta, Float const phi)
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
 * Custom designed functions - partially present in [Vlnas2025].
 */

// 14.
template<typename Float>
Float cf_f1(Float const theta, Float const phi)
{
    // (Abs[Sin[Cos[2 \[Phi]] - 2 \[Theta]]] + Abs[Cos[2 \[Theta]]])
    return std::abs(std::sin(std::cos(2.0f * phi) - 2.0f * theta)) + std::abs(std::cos(2.0f * theta));
}

// 15.
template<typename Float>
Float cf_f2(Float const theta, Float const phi)
{
    // Abs[Sin[2 \[Phi] - \[Theta]]] + Abs[Cos[2 \[Theta]]]
    return std::abs(std::sin(2.0f * phi - theta)) + std::abs(std::cos(2.0f * theta));
}

// 16.
template<typename Float>
Float cf_f3(Float const /*theta*/, Float const phi)
{
    // 1 + Sin[5 \[Phi]]/5
    return 1.0f + std::sin(5.0f * phi) / 5.0f;
}

// 17.
template<typename Float>
Float cf_f4(Float const theta, Float const phi)
{
    // 1 + Cos[5 \[Phi]]/5 + Sin[5 \[Theta]]
    return 1.0f + std::cos(5.0f * phi) / 5.0f + std::sin(5.0f * theta);
}

// 18.
template<typename Float>
Float cf_f5(Float const theta, Float const phi)
{
    // f_{n,5}(\theta, \phi) = 1.5\exp\left(0.8\left(\sin\theta \cos\phi + 0.5 \cos\theta\right)\right)
    // + 1.2\exp\left(0.6\left(-\sin\theta \sin\phi + 0.3 \cos\theta\right)\right)
    // + 0.8\exp\left(0.5 \cos\theta\right) + 0.5\left(1 + \cos(6\theta)\sin(4\phi)\right)
    return 1.5f * std::exp(0.8f * (std::sin(theta) * std::cos(phi) + 0.5f * std::cos(theta))) +
           1.2f * std::exp(0.6f * (-std::sin(theta) * std::sin(phi) + 0.3f * std::cos(theta))) +
           0.8f * std::exp(0.5f * std::cos(theta)) + 0.5f * (1.0f + std::cos(6.0f * theta) * std::sin(4.0f * phi));
}

// 19.
template<typename Float>
Float cf_f6(Float const theta, Float const phi)
{
    return 1.0f + 0.5f * std::cos(theta) + 0.3f * std::cos(2.0f * phi);
}

// 20.
template<typename Float>
Float cf_f7(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::abs(std::cos(3.0f * x) + std::sin(2.0f * y) + 0.5f * z * z);
}

// 21.
template<typename Float>
Float cf_f8(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::abs(std::sin(2.0f * x) * std::cos(3.0f * y) + 0.5f * z * z + 0.3f * std::sin(5.0f * x) * std::cos(4.0f * z));
}

// 22.
template<typename Float>
Float cf_f9(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::abs(x * x - y * y + 0.5f * x * z - 0.3f * y * z);
}

// 23.
template<typename Float>
Float cf_f10(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return x * x + y * y + z * z + 5 + 2.5f * std::cos((theta - F_PI) / 2.0f) * std::sin(16.0f * theta);
}

// 24.
template<typename Float>
Float cf_f11(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    // f(x,y,z)=|sin(10x)cos(12y)sin(15z)+cos(20x)|
    return std::abs(std::sin(10.0f * x) * std::cos(12.0f * y) * std::sin(15.0f * z) + std::cos(20.0f * x));
}

// 25.
template<typename Float>
Float cf_f12(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    // f(x,y,z)=sin(10x)+cos(12y)−sin(15z)+0.2cos(18x) + 3
    return std::sin(10.0f * x) + std::cos(12.0f * y) - std::sin(15.0f * z) + 0.2f * std::cos(18.0f * x) + 3.0f;
}

// 26.
template<typename Float>
Float cf_f13(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::exp(-std::sin(5.0f * x) - std::cos(6.0f * y)) + 0.3f * std::sin(10.0f * z);
}

// 27.
template<typename Float>
Float cf_f14(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return std::exp(-2.0f * (x * x + y * y)) * std::sin(4.0f * z);
}

// 28.
template<typename Float>
Float cf_15(Float const theta, Float const phi)
{
    auto [x, y, z] = spherical_to_xyz(theta, phi);
    return (x * x + y * y) * std::exp(-3.0f * z * z);
}

enum FunctionIndexing : uint8_t
{
    FEX_MAX = 13,
    FCF_MAX = FEX_MAX + 15 + 1
};

template<typename Float>
std::function<Float(Float, Float)> get_function(uint8_t id)
{
    static std::unordered_map<uint8_t, std::function<Float(Float, Float)>> functions =
    {
        // Functions based on literature
        { 1, &sphc::fornberg_f1<Float> },
        { 2, &sphc::fornberg_f4<Float> },
        { 3, &sphc::beentjes_f3<Float> },
        { 4, &sphc::beentjes_f4<Float> },
        { 5, &sphc::beentjes_f5<Float> },
        { 6, &sphc::renka_f3<Float> },
        { 7, &sphc::renka_f4<Float> },
        { 8, &sphc::renka_f5<Float> },
        { 9, &sphc::reegar_f3<Float> },
        { 10, &sphc::bellet_f4<Float> },
        { 11, &sphc::reegar_f2<Float> },
        { 12, &sphc::reegar_f4<Float> },
        { 13, &sphc::fornberg_f2<Float> },
        // Custom functions
        { FEX_MAX + 1, &sphc::cf_f1<Float> },
        { FEX_MAX + 2, &sphc::cf_f2<Float> },
        { FEX_MAX + 3, &sphc::cf_f3<Float> },
        { FEX_MAX + 4, &sphc::cf_f4<Float> },
        { FEX_MAX + 5, &sphc::cf_f5<Float> },
        { FEX_MAX + 6, &sphc::cf_f6<Float> },
        { FEX_MAX + 7, &sphc::cf_f7<Float> },
        { FEX_MAX + 8, &sphc::cf_f8<Float> },
        { FEX_MAX + 9, &sphc::cf_f9<Float> },
        { FEX_MAX + 10, &sphc::cf_f10<Float> },
        { FEX_MAX + 11, &sphc::cf_f11<Float> },
        { FEX_MAX + 12, &sphc::cf_f12<Float> },
        { FEX_MAX + 13, &sphc::cf_f13<Float> },
        { FEX_MAX + 14, &sphc::cf_f14<Float> },
        { FEX_MAX + 15, &sphc::cf_15<Float> }
    };

    return functions[id];
}

} // namespace sphc

#endif // SPHERICAL_COLLECTION_FUNCTIONS_H
