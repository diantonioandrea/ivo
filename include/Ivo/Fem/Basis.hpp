/**
 * @file Basis.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Basis functions.
 * @date 2024-09-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef FEM_BASIS
#define FEM_BASIS

#include "./Legendre.hpp"

namespace ivo {

    namespace internal {

        // Maps.

        std::tuple<std::array<Vector<Real>, 2>, Real> reference_to_element(const Mesh21 &, const Natural &, const std::array<Vector<Real>, 2> &);

    }

    // Basis functions.

    std::array<Matrix<Real>, 2> basis_t(const Mesh21 &, const Natural &, const Vector<Real> &);
    std::array<Matrix<Real>, 3> basis_s(const Mesh21 &, const Natural &, const std::array<Vector<Real>, 2> &);

}

#endif