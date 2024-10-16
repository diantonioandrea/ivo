/**
 * @file Legendre.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Legendre polynomials evaluation.
 * @date 2024-08-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FEM_LEGENDRE
#define FEM_LEGENDRE

#include "./Includes.hpp"

namespace ivo {

    // Legendre polynomials.

    namespace internal {

        Vector<Real> legendre(const Vector<Real> &, const Natural &, const Natural &k = 0);

    }

}

#endif