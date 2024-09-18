/**
 * @file Forcing.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Problem's forcing vector.
 * @date 2024-09-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PROBLEM_FORCING
#define PROBLEM_FORCING

#include "./Equation.hpp"
#include "./Initial.hpp"

namespace ivo {

    Vector<Real> forcing(const Mesh21 &, const Equation &, const Data &, const Initial &);

}

#endif