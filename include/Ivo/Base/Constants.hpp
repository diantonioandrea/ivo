/**
 * @file Constants.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Constants.
 * @date 2024-07-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASE_CONSTANTS
#define BASE_CONSTANTS

#include "./Primitives.hpp"

namespace ivo {

    namespace constants {

        // Zeros.

        /**
         * @brief Zero tolerance.
         * 
         */
        constexpr Real zero = 1E-20;

        /**
         * @brief Algebra zero tolerance (solvers).
         * 
         */
        constexpr Real algebra_zero = 1E-12;

        /**
         * @brief Geometrical zero tolerance.
         * 
         */
        constexpr Real geometry_zero = 1E-12;
        
        /**
         * @brief Quadrature zero tolerance.
         * 
         */
        constexpr Real quadrature_zero = 1E-14;

        // Diagrams.

        /**
         * @brief Lloyd stopping multiplier.
         * 
         */
        constexpr Real diagram_stop = 1E-4;

        /**
         * @brief Collapse multiplier.
         * 
         */
        constexpr Real diagram_collapse = 1E-1;

        // Constants.

        /**
         * @brief Quadrature order.
         * 
         */
        constexpr Natural quadrature = 5;

        // Solvers.

        /**
         * @brief Solvers' maximum number of iterations.
         * 
         */
        constexpr Natural solvers_stop = 1E4;

        /**
         * @brief Restarted GMRES threshold.
         * 
         */
        constexpr Natural gmres_restart = 2E2;

    }

}

#endif