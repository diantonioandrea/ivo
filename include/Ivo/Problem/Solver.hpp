/**
 * @file Solver.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Problem's solver.
 * @date 2024-10-07
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PROBLEM_SOLVER
#define PROBLEM_SOLVER

#include "./Stiffness.hpp"
#include "./Forcing.hpp"

namespace ivo {

    Vector<Real> solve(const Mesh21 &, const Sparse<Real> &, const Vector<Real> &, const Initial &);

}

#endif