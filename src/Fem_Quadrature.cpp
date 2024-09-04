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
    std::array<Vector<Real>, 2> __gauss1(const Natural &n, const Real &a, const Real &b) {
        #ifndef NDEBUG // Integrity check.
        assert(a < b);
        assert(n % 2 == 1); // n must be odd.
        #endif

        // Nodes to be computed.
        Natural m = (n + 1) / 2;

        // Initialization. Operation vector.
        Vector<Real> z = __vcos(M_PI * (stepped<Real>(1, m) - 0.25L) / (n + 0.5L));

        // Error.
        Real error_scalar = 1.0L + constants::quadrature_zero;
        Vector<Real> error_vector{m, error_scalar};

        // Temp vector.
        Vector<Real> temp{m};

        // Algorithm.
        while(error_scalar > constants::quadrature_zero) {

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
            
            for(Natural j = 0; j < m; ++j) {
                Real update = (error_vector(j) > constants::quadrature_zero) ? p(j, 0) / temp(j) : 0.0L;
                z(j, z_old(j) - update);
            }

            // Error update.
            error_vector = std::abs(z - z_old);
            error_scalar = max(error_vector);
        }

        // Nodes.
        Mask mask{m, true};
        mask(m - 1, false);

        // Nodes.
        Vector<Real> nodes = ((b + a) / 2.0L - (b - a) / 2.0L) * stacked(z(mask), -flipped(z));

        // Reflection.
        Vector<Real> z_ref = stacked(z(mask), flipped(z));
        Vector<Real> temp_ref = stacked(temp(mask), flipped(temp));

        // Weights.
        Vector<Real> weights = (b - a) / (1.0L - z_ref * z_ref) / (temp_ref * temp_ref);

        return {nodes, weights};
    }

    // Gauss-Legendre nodes and weights over reference structures.

    /**
     * @brief Gauss-Legendre nodes and weights over the reference interval [-1, 1].
     * 
     * @param n Order.
     * @return std::array<Vector<Real>, 2> 
     */
    std::array<Vector<Real>, 2> quadrature1(const Natural &n) { return __gauss1(n, -1.0L, 1.0L); }

    /**
     * @brief Gauss-Legendre nodes and weights over the reference triangle {(0, 0), (0, 1), (1, 0)}.
     * 
     * @param n 
     * @return std::array<Vector<Real>, 3> 
     */
    std::array<Vector<Real>, 3> quadrature2(const Natural &n) {

        // Square nodes and weights.
        auto [nodes1, weights1] = quadrature1(n);

        // Temporary nodes and weights.
        Vector<Real> nodes2_x{n * n};
        Vector<Real> weights2_x{n * n};

        Vector<Real> nodes2_y{n * n};
        Vector<Real> weights2_y{n * n};

        for(Natural j = 0; j < n; ++j)
            for(Natural k = 0; k < n; ++k) {
                
                // Nodes.
                nodes2_x(j * n + k, nodes1(j));
                nodes2_y(j * n + k, nodes1(k));

                // Weights.
                weights2_x(j * n + k, weights1(j));
                weights2_y(j * n + k, weights1(k));
            }
        
        // Triangle nodes and weights.
        Vector<Real> nodes2t_x = (1.0L + nodes2_x) / 2.0L;
        Vector<Real> nodes2t_y = (1.0L - nodes2_x) * (1.0L + nodes2_y) / 4.0L;
        Vector<Real> weights2t = (1.0L - nodes2_x) * weights2_x * weights2_y / 8.0L;

        return {nodes2t_x, nodes2t_y, weights2t};
    }

}