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

#include "../Vector.hpp"

namespace ivo {

    /**
     * @brief Dot product.
     * 
     * @tparam T Numerical type.
     * @param x vector.
     * @param y vector.
     * @return T 
     */
    template<Numerical T>
    T dot(const Vector<T> &x, const Vector<T> &y) {
        #ifndef NDEBUG // Integrity check.
        assert(x.size() == y.size());
        #endif

        std::vector<T> x_entries = x.entries();
        std::vector<T> y_entries = y.entries();

        if constexpr (Conjugable<T>)
            return std::transform_reduce(x_entries.begin(), x_entries.end(), y_entries.begin(), static_cast<T>(0), std::plus{}, [](const T &x_entry, const T &y_entry){ return x_entry * std::conj(y_entry); });

        return std::transform_reduce(x_entries.begin(), x_entries.end(), y_entries.begin(), static_cast<T>(0), std::plus{}, [](const T &x_entry, const T &y_entry){ return x_entry * y_entry; });
    }

    /**
     * @brief Cross product.
     * 
     * @tparam T Numerical type.
     * @param x First vector.
     * @param y Second vector.
     * @return Vector<T> 
     */
    template<Numerical T>
    Vector<T> cross(const Vector<T> &x, const Vector<T> &y) {
        #ifndef NDEBUG // Integrity check.
        assert(x.size() == 3);
        assert(y.size() == 3);
        #endif

        return Vector<T>{{x(1) * y(2) - x(2) * y(1), x(2) * y(0) - x(0) * y(2), x(0) * y(1) - x(1) * y(0)}};
    }

    /**
     * @brief Vectorial norm.
     * 
     * @tparam T Numerical type.
     * @param x Vector.
     * @return Real 
     */
    template<Numerical T>
    Real norm(const Vector<T> &x) {
        std::vector<T> x_entries = x.entries();
        return std::sqrt(std::transform_reduce(x_entries.begin(), x_entries.end(), static_cast<T>(0), std::plus{}, [](const T &entry){ return std::abs(entry) * std::abs(entry); }));
    }
    
    /**
     * @brief Min{vector}.
     * 
     * @tparam T Numerical type.
     * @param x Vector.
     * @return T 
     */
    template<Numerical T>
    T min(const Vector<T> &x) {
        std::vector<T> x_entries = x.entries();
        return *std::min_element(x_entries.begin(), x_entries.end());
    }

    /**
     * @brief Max{vector}.
     * 
     * @tparam T Numerical type.
     * @param x Vector.
     * @return T 
     */
    template<Numerical T>
    T max(const Vector<T> &x) {
        std::vector<T> x_entries = x.entries();
        return *std::max_element(x_entries.begin(), x_entries.end());
    }

    /**
     * @brief Flips a vector.
     * 
     * @tparam T Numerical type.
     * @param x Vector.
     * @return Vector<T> 
     */
    template<Numerical T>
    Vector<T> flipped(const Vector<T> &x) {
        Vector<T> _flipped{x.size()};

        for(Natural j = 0; j < x.size(); ++j)
            _flipped(j, x(x.size() - j - 1));

        return _flipped;
    }

    /**
     * @brief Stacks two vectors.
     * 
     * @tparam T Numerical type.
     * @param x First vector.
     * @param y Second vector.
     * @return Vector<T> 
     */
    template<Numerical T>
    Vector<T> stacked(const Vector<T> &x, const Vector<T> &y) {
        Vector<T> _stacked{x.size() + y.size()};

        for(Natural j = 0; j < x.size(); ++j)
            _stacked(j, x(j));
        
        for(Natural j = 0; j < y.size(); ++j)
            _stacked(x.size() + j, y(j));

        return _stacked;
    }

    /**
     * @brief Stepped vector [a, a + step, ..., b].
     * 
     * @tparam T Numerical type.
     * @param a Start.
     * @param b End.
     * @param step Step.
     * @return Vector<Real> 
     */
    template<Numerical T>
    Vector<T> stepped(const T &a, const T &b, const T &step = static_cast<T>(1)) {
        #ifndef NDEBUG // Integrity check.
        assert(((a < b) && (step > NUMERICAL_ZERO)) || ((a > b) && (step < -NUMERICAL_ZERO)));
        #endif

        std::vector<T> _stepped;

        if(step > NUMERICAL_ZERO)
            for(T j = a; j <= b + NUMERICAL_ZERO; j += step)
                _stepped.emplace_back(j);
        else
            for(T j = a; j >= b - NUMERICAL_ZERO; j += step)
                _stepped.emplace_back(j);

        return Vector{_stepped};
    }

}

namespace std {

    /**
     * @brief Vectorial std::abs.
     * 
     * @tparam T Numerical type.
     * @param x Vector.
     * @return Vector<T> 
     */
    template<ivo::Numerical T>
    ivo::Vector<T> abs(const ivo::Vector<T> &x) {
        ivo::Vector<T> result(x.size());

        for(ivo::Natural j = 0; j < x.size(); ++j)
            result(j, std::abs(x(j)));

        return result;
    }

}

#endif