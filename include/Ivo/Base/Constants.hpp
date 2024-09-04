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

    // Zeros.

    /**
     * @brief Zero tolerance.
     * 
     */
    constexpr Real ___zero = 1E-16;

    /**
     * @brief Geometrical zero tolerance.
     * 
     */
    constexpr Real ___geometrical_zero = 1E-12;
    
    /**
     * @brief Quadrature zero tolerance.
     * 
     */
    constexpr Real ___quadrature_zero = 1E-14;

    // Diagrams.

    /**
     * @brief Lloyd stopping multiplier.
     * 
     */
    constexpr Real ___diagram_stop = 1E-6;

    /**
     * @brief Collapse multiplier.
     * 
     */
    constexpr Real ___diagram_collapse = 1E-1;

    // Constants.

    /**
     * @brief Quadrature order.
     * 
     */
    constexpr Natural ___quadrature = 7;

}

#endif