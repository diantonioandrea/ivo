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

#include <Ivo.hpp>

namespace ivo {
    namespace square {

        // Domain.
        
        const ivo::Point21 a{0.0, 0.0};
        const ivo::Point21 b{1.0, 0.0};
        const ivo::Point21 c{1.0, 1.0};
        const ivo::Point21 d{0.0, 1.0};

        const ivo::Polygon21 abcd = {a, b, c, d};

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
            return std::sin(x) * std::sin(y) * std::sin(t);
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
            return {std::cos(x) * std::sin(y) * std::sin(t), std::sin(x) * std::cos(y) * std::sin(t)};
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
            return std::sin(x) * std::sin(y) * std::cos(t);
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
            return -2.0L * u(x, y, t);
        }

        // Coefficients.

        /**
         * @brief Convection coefficient.
         * 
         * @param t 
         * @return std::array<Real, 2> 
         */
        std::array<Real, 2> convection(const Real &t) {
            return {0.5L, 0.5L};
        }

        /**
         * @brief Diffusion coefficient.
         * 
         * @param t 
         * @return Real 
         */
        Real diffusion(const Real &t) {
            return 1.0L;
        }

        /**
         * @brief Reaction coefficient.
         * 
         * @param t 
         * @return Real 
         */
        Real reaction(const Real &t) {
            return 1.0L;
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
            auto [u_x, u_y] = u_xy(x, y, t);

            if(x <= ivo::constants::algebra_zero)
                return -1.0L * diffusion(t) * u_x;

            if(x >= 1.0L - ivo::constants::algebra_zero)
                return diffusion(t) * u_x;
            
            if(y <= ivo::constants::algebra_zero)
                return -1.0L * diffusion(t) * u_y;
            
            return diffusion(t) * u_y;
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
            auto [u_x, u_y] = u_xy(x, y, t);
            auto [convection_x, convection_y] = convection(t);

            return u_t(x, y, t) - diffusion(t) * u_xxyy(x, y, t) + convection_x * u_x + convection_y * u_y + reaction(t) * u(x, y, t);
        }

    }
}

#endif