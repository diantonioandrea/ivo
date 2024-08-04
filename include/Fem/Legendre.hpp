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

    Vector<Real> legendre1(const Vector<Real> &, const Natural &);
    Vector<Real> legendre_grad1(const Vector<Real> &, const Natural &);

}

#endif