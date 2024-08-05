/**
 * @file Fem_Quadrature.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Fem/Quadrature.hpp implementation.
 * @date 2024-08-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Utilities.

    Vector<Real> __vcos(const Vector<Real> &vector) {
        Vector<Real> result{vector.size()};

        for(Natural j = 0; j < vector.size(); ++j)
            result(j, std::cos(vector(j)));

        return result;
    }

    // Gauss-Legendre nodes.

    /**
     * @brief 1D Gauss-Legendre nodes and weights over [a, b].
     * 
     * @param n Order.
     * @param a Interval's a.
     * @param b Interval's b.
     * @return std::array<Vector<Real>, 2> 
     */
    std::array<Vector<Real>, 2> gauss1(const Natural &n, const Real &a, const Real &b) {
        #ifndef NDEBUG // Integrity check.
        assert(a < b);
        assert(n % 2 == 1); // n must be odd.
        #endif

        // Nodes to be computed.
        Natural m = (n + 1) / 2;

        // Initialization. Operation vector.
        Vector<Real> z = __vcos(M_PI * (stepped<Real>(1, m) - 0.25L) / (n + 0.5L));

        // Error.
        Real error_scalar = 1.0L + QUADRATURE_ZERO;
        Vector<Real> error_vector{m, error_scalar};

        // Temp vector.
        Vector<Real> temp{m};

        // Algorithm.
        while(error_scalar > QUADRATURE_ZERO) {

            // Operation matrix.
            Matrix<Real> p{m, 3};
            p.column(0, Vector<Real>{m, 1.0L});

            // Iterations.
            for(Real j = 1; j <= n; j += 1.0L) {

                // Copy.
                p.column(2, p.column(1));
                p.column(1, p.column(0));

                // Update.
                p.column(0, (2.0L * j - 1.0L) / j * z * p.column(1) - (j - 1.0L) / j * p.column(2));
            }

            // Update step.
            temp = n * (z * p.column(0) - p.column(1)) / (z * z - 1.0L);
            Vector<Real> z_old{z};
            
            for(Natural j = 0; j < m; ++j)
                if((error_vector(j) > QUADRATURE_ZERO) && (temp(j) > NUMERICAL_ZERO))
                    z(j, z_old(j) - p(j, 1) / temp(j));

            // Error update.
            error_vector = std::abs(z - z_old);
            error_scalar = max(error_vector);
        }

        // Nodes.
        Mask mask{m, true};
        mask(m - 1, false);

        Vector<Real> nodes = ((b + a) / 2.0L - (b - a) / 2.0L) * stacked(z(mask), -flipped(z));

        // Reflection.
        Vector<Real> z_ref = stacked(z(mask), flipped(z));
        Vector<Real> temp_ref = stacked(temp(mask), flipped(temp));

        // Weights.
        Vector<Real> weights = (b - a) / (1.0L - z_ref * z_ref) / (temp_ref * temp_ref);

        return {nodes, weights};
    }

}