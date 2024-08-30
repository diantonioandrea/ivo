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

    // Gauss-Legendre nodes.
    
    std::array<Vector<Real>, 2> gauss1(const Natural &, const Real &, const Real &);

}

#endif