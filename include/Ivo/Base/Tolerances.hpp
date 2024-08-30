/**
 * @file Tolerances.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Tolerances.
 * @date 2024-07-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BASE_TOLERANCES
#define BASE_TOLERANCES

// Base zero tolerance.
#ifndef NUMERICAL_ZERO
#define NUMERICAL_ZERO 1E-16
#endif

// Geometrical zero tolerance.
#ifndef GEOMETRICAL_ZERO
#define GEOMETRICAL_ZERO 1E-12
#endif

// Quadrature zero tolerance.
#ifndef QUADRATURE_ZERO
#define QUADRATURE_ZERO 1E-14
#endif

// Lloyd stopping multiplier.
#ifndef LLOYD_STOP
#define LLOYD_STOP 1E-8
#endif

// Collapse multiplier.
#ifndef COLLAPSE
#define COLLAPSE 1E-1
#endif

#endif