/**
 * @file Square.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Square functions and data.
 * @date 2024-10-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef SRC_SQUARE
#define SRC_SQUARE

// Testing.
#include "./Test.hpp"

#include <Ivo.hpp>

namespace ivo {
    namespace square {

        // Domain.
        
        const ivo::Point21 a{0.0, 0.0};
        const ivo::Point21 b{1.0, 0.0};
        const ivo::Point21 c{1.0, 1.0};
        const ivo::Point21 d{0.0, 1.0};

        const ivo::Polygon21 abcd = {a, b, c, d};

        // Coefficients.

        /**
         * @brief Convection coefficient.
         * 
         * @param t 
         * @return std::array<Real, 2> 
         */
        std::array<Real, 2> convection(const Real &x, const Real &y, const Real &t) {
            return {1.0L, 1.0L};
        }

        /**
         * @brief Diffusion coefficient.
         * 
         */
        const Real diffusion = 0.005L;

        /**
         * @brief Reaction coefficient.
         * 
         * @param t 
         * @return Real 
         */
        Real reaction(const Real &x, const Real &y, const Real &t) {
            return 0.5L;
        }

        /**
         * @brief Boundary layer coefficient.
         * 
         */
        const Real boundary = 0.05L;

        // Solution.
        
        /**
         * @brief Exact solution.
         * 
         * @param x 
         * @param y 
         * @param t 
         * @return Real 
         */
        Real u(const Real &x, const Real &y, const Real &t) {
            const Real b = boundary;
            const auto [c_x, c_y] = convection(x, y, t);

            return (1.0L - std::exp(-t)) * (2.0L * x + 2.0L * y - x * y * (1.0L - std::exp(c_x * (x - 1) / b)) * (1.0L - std::exp(c_y * (y - 1) / b)));
        }

        /**
         * @brief Exact solution's gradient.
         * 
         * @param x 
         * @param y 
         * @param t 
         * @return std::array<Real, 2> 
         */
        std::array<Real, 2> u_xy(const Real &x, const Real &y, const Real &t) {
            const Real b = boundary;
            const auto [c_x, c_y] = convection(x, y, t);

            return {(1.0L - std::exp(-t)) * (2.0L - (1.0L - std::exp(c_y * (y - 1) / b)) * (y * (1.0L - std::exp(c_x * (x - 1) / b)) + x * y * (-c_x / b * std::exp(c_x * (x - 1) / b)))), (1.0L - std::exp(-t)) * (2.0L - (1.0L - std::exp(c_x * (x - 1) / b)) * (x * (1.0L - std::exp(c_y * (y - 1) / b)) + x * y * (-c_y / b * std::exp(c_y * (y - 1) / b))))};
        }
        
        /**
         * @brief Exact solution's time derivative.
         * 
         * @param x 
         * @param y 
         * @param t 
         * @return Real 
         */
        Real u_t(const Real &x, const Real &y, const Real &t) {
            const Real b = boundary;
            const auto [c_x, c_y] = convection(x, y, t);

            return std::exp(-t) * (2.0L * x + 2.0L * y - x * y * (1.0L - std::exp(c_x * (x - 1) / b)) * (1.0L - std::exp(c_y * (y - 1) / b)));
        }

        /**
         * @brief Exact solution's laplacian.
         * 
         * @param x 
         * @param y 
         * @param t 
         * @return Real 
         */
        Real u_xxyy(const Real &x, const Real &y, const Real &t) {
            const Real b = boundary;
            const auto [c_x, c_y] = convection(x, y, t);
            const Real cb_x = c_x / b;
            const Real cb_y = c_y / b;

            return (std::exp(-t) - 1.0L) * (x * (1.0L - std::exp(cb_x * (x - 1.0L))) * (2.0L * (-cb_y * std::exp(cb_y * (y - 1.0L))) + y * (-cb_y * cb_y * std::exp(cb_y * (y - 1.0L)))) + y * (1.0L - std::exp(cb_y * (y - 1.0L))) * (2.0L * (-cb_x * std::exp(cb_x * (x - 1.0L))) + x * (-cb_x * cb_x * std::exp(cb_x * (x - 1.0L)))));
        }

        // Data.

        /**
         * @brief Initial condition.
         * 
         * @param x 
         * @param y 
         * @return Real 
         */
        Real u0(const Real &x, const Real &y) {
            return u(x, y, 0.0L);
        }

        /**
         * @brief Dirichlet boundary condition.
         * 
         * @param x 
         * @param y 
         * @param t 
         * @return Real 
         */
        Real gd(const Real &x, const Real &y, const Real &t) {
            return u(x, y, t);
        }

        /**
         * @brief Neumann boundary condition.
         * 
         * @param x 
         * @param y 
         * @param t 
         * @return Real 
         */
        Real gn(const Real &x, const Real &y, const Real &t) {
            const auto [u_x, u_y] = u_xy(x, y, t);

            if(x <= ivo::constants::algebra_zero)
                return -1.0L * diffusion * u_x;

            if(x >= 1.0L - ivo::constants::algebra_zero)
                return diffusion * u_x;
            
            if(y <= ivo::constants::algebra_zero)
                return -1.0L * diffusion * u_y;
            
            return diffusion * u_y;
        }

        // Source.

        /**
         * @brief Source.
         * 
         * @param x 
         * @param y 
         * @param t 
         * @return Real 
         */
        Real g(const Real &x, const Real &y, const Real &t) {
            const auto [u_x, u_y] = u_xy(x, y, t);
            const auto [c_x, c_y] = convection(x, y, t);

            return u_t(x, y, t) - diffusion * u_xxyy(x, y, t) + c_x * u_x + c_y * u_y + reaction(x, y, t) * u(x, y, t);
        }

    }
}

#endif