/**
 * @file Fem_Legendre.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Fem/Legendre.hpp implementation.
 * @date 2024-08-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Utilities.

    /**
     * @brief Binomial coefficient.
     * 
     * @param n Natural.
     * @param k Natural.
     * @return Natural 
     */
    Natural __binomial(const Natural &n, const Natural &k) {
        #ifndef NDEBUG // Integrity check.
        assert(n >= k);
        #endif

        if((n == k) || (k == 0))
            return 1;

        // Recursive formula.
        return __binomial(n - 1, k - 1) + __binomial(n - 1, k);
    }

    // Legendre polynomials.

    /**
     * @brief 1D Legendre polynomials evaluation over a vector of points.
     * 
     * @param x Real points.
     * @param n Degree.
     * @return Vector<Real> 
     */
    Vector<Real> legendre1(const Vector<Real> &x, const Natural &n) {
        Vector<Real> legendre{x.size()};

        // Recursive formula.
        for(Natural k = 0; k <= n; ++k) {
            Vector<Real> power{x.size(), 1.0L};

            for(Natural j = 0; j < k; ++j)
                power *= 0.5L * (x - 1.0L);

            legendre += __binomial(n, k) * __binomial(n + k, k) * power;
        }

        return legendre;
    }

    /**
     * @brief 1D Legendre polynomials evaluation over a vector of points. First derivative.
     * 
     * @param x Real points.
     * @param n Degree.
     * @return Vector<Real> 
     */
    Vector<Real> legendre_grad1(const Vector<Real> &x, const Natural &n) {
        Vector<Real> legendre{x.size()};

        // Recursive formula.
        for(Natural k = 1; k <= n; ++k) {
            Vector<Real> power{x.size(), 1.0L};

            for(Natural j = 1; j < k; ++j)
                power *= 0.5L * (x - 1.0L);

            legendre += 0.5L * __binomial(n, k) * __binomial(n + k, k) * power;
        }

        return legendre;
    }

}