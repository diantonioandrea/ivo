/**
 * @file Stiffness.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Problem's stiffness matrix.
 * @date 2024-08-07
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PROBLEM_STIFFNESS
#define PROBLEM_STIFFNESS

#include "./Equation.hpp"

namespace ivo {

    Sparse<Real> stiffness(const Mesh21 &, const Equation &, const Initial &);

}

#endif