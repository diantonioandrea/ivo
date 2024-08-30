/**
 * @file Concepts.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Main concepts.
 * @date 2024-07-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASE_CONCEPTS
#define BASE_CONCEPTS

#include "./Includes.hpp"
#include "./Primitives.hpp"

namespace ivo {

    // Numerical objects.

    /**
     * @brief Summable objects.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Summable = requires(T x, T y) {
        {x += y} -> std::convertible_to<T>;
        {x -= y} -> std::convertible_to<T>;
    };

    /**
     * @brief Multipliable objects.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Multipliable = requires(T x, T y) {
        {x *= y} -> std::convertible_to<T>;
        {x /= y} -> std::convertible_to<T>;
    };

    /**
     * @brief "Normable" objects.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Normable = requires(T x) {
        {std::abs(x)} -> std::convertible_to<Real>;
    };

    /**
     * @brief Numerical object.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Numerical = Summable<T> && Multipliable<T> && Normable<T> && std::convertible_to<T, Real>;

    // Other concepts.

    /**
     * @brief Conjugable objects.
     * 
     * @tparam T 
     */
    template<typename T>
    concept Conjugable = requires(T x) {
        {std::conj(x)} -> std::convertible_to<T>;
    };

}

#endif