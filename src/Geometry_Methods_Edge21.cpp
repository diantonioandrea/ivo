/**
 * @file Geometry_Methods_Edge21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Methods/Edge21.hpp implementation.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Distances.

    /**
     * @brief Distance(edge, point).
     * 
     * @param ab Edge.
     * @param p Point.
     * @return Real 
     */
    Real distance(const Edge21 &ab, const Point21 &p);

    /**
     * @brief Distance(edge, edge).
     * 
     * @param ab 
     * @param cd 
     * @return Real 
     */
    Real distance(const Edge21 &ab, const Edge21 &cd);

    // Containment.

    /**
     * @brief Contains(edge, point).
     * Evaluated by distances.
     * 
     * @param ab Edge.
     * @param p Point.
     * @return true 
     * @return false 
     */
    bool contains(const Edge21 &ab, const Point21 &p) {
        Real ap = distance(ab(0), p);
        Real bp = distance(ab(1), p);

        return std::abs(distance(ab(0), ab(1)) - (ap + bp)) <= GEOMETRICAL_ZERO;
    }

    /**
     * @brief Contains(edge, edge).
     * 
     * @param ab Edge.
     * @param cd Edge.
     * @return true 
     * @return false 
     */
    bool contains(const Edge21 &ab, const Edge21 &cd) {
        return contains(ab, cd(0)) && contains(ab, cd(1));
    }

}