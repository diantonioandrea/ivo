/**
 * @file Geometry_Methods_Point21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Methods/Point21.hpp implementation.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    /**
     * @brief Distance(point, point).
     * 
     * @param p Point.
     * @param q Point.
     * @return Real 
     */
    Real distance(const Point21 &p, const Point21 &q) { return std::sqrt((p(0) - q(0)) * (p(0) - q(0)) + (p(1) - q(1)) * (p(1) - q(1)) + (p(2) - q(2)) * (p(2) - q(2))); }

}