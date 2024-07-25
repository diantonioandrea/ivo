/**
 * @file Vector.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Vector related methods.
 * @date 2024-07-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ALGEBRA_METHODS_VECTOR
#define ALGEBRA_METHODS_VECTOR

#include "../Includes.hpp"
#include "../Vector.hpp"

namespace ivo {

    /**
     * @brief Vectorial dot product.
     * 
     * @tparam T Numerical type.
     * @param x First vector.
     * @param y Second vector.
     * @return T 
     */
    template<Numerical T>
    T dot(const Vector<T> &x, const Vector<T> &y) {
        std::vector<T> x_entries = x.entries();
        std::vector<T> y_entries = y.entries();

        if constexpr (Conjugable<T>)
            return std::transform_reduce(x_entries.begin(), x_entries.end(), y_entries.begin(), static_cast<T>(0), std::plus{}, [](const T &x_entry, const T &y_entry){ return x_entry * std::conj(y_entry); });

        return std::transform_reduce(x_entries.begin(), x_entries.end(), y_entries.begin(), static_cast<T>(0), std::plus{}, [](const T &x_entry, const T &y_entry){ return x_entry * y_entry; });
    }

    /**
     * @brief Vectorial norm.
     * 
     * @tparam T Numerical type.
     * @param x Vector.
     * @return T 
     */
    template<Numerical T>
    Real norm(const Vector<T> &x) {
        std::vector<T> x_entries = x.entries();
        return std::sqrt(std::transform_reduce(x_entries.begin(), x_entries.end(), static_cast<T>(0), std::plus{}, [](const T &entry){ return std::abs(entry) * std::abs(entry); }));
    }

}

#endif