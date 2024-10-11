/**
 * @file Quadrature.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Quadrature methods.
 * @date 2024-08-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FEM_QUADRATURE
#define FEM_QUADRATURE

#include "./Includes.hpp"

namespace ivo {

    namespace internal {

        // Gauss-Legendre nodes and weights.
    
        std::array<Vector<Real>, 2> gauss1(const Natural &, const Real &, const Real &);

    }

    // Gauss-Legendre nodes and weights over reference structures.

    std::array<Vector<Real>, 2> quadrature1t(const Natural &);
    std::array<Vector<Real>, 2> quadrature1x(const Natural &);
    std::array<Vector<Real>, 3> quadrature2xy(const Natural &);

}

#endif