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

    namespace internal {

        /**
         * @brief Binomial coefficient.
         * 
         * @param n Natural.
         * @param k Natural.
         * @return Natural 
         */
        Natural binomial(const Natural &n, const Natural &k) {
            #ifndef NDEBUG // Integrity check.
            assert(n >= k);
            #endif

            if((n == k) || (k == 0))
                return 1;

            // Recursive formula.
            return internal::binomial(n - 1, k - 1) + internal::binomial(n - 1, k);
        }

        /**
         * @brief Legendre polynomials evaluation over a vector of points.
         * 
         * @param x Real points.
         * @param n Degree.
         * @param k Derivative.
         * @return Vector<Real> 
         */
        Vector<Real> legendre(const Vector<Real> &x, const Natural &n, const Natural &k) {
            Vector<Real> y{x.size()};

            // Recursive formula.
            for(Natural j = k; j <= n; ++j) {
                Vector<Real> y_j{x.size(), 1.0L};

                for(Natural h = 0; h < j - k; ++h)
                    y_j *= 0.5L * (x - 1.0L);

                for(Natural h = 0; h < k; ++h)
                    y_j *= 0.5L * static_cast<Real>(j - h);

                y += internal::binomial(n, j) * internal::binomial(n + j, j) * y_j;
            }

            return y;
        }

    }

}