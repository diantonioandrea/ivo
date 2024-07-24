/**
 * @file Edge21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Edge21 related methods.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY21_METHODS_EDGE21
#define GEOMETRY21_METHODS_EDGE21

#include "../Includes.hpp"
#include "../Edge21.hpp"

namespace ivo {

    // Intersections.

    std::optional<Point21> intersections(const Edge21 &, const Edge21 &);

    // Distances.

    Real distance(const Edge21 &, const Point21 &);

    // Containment.

    bool contains(const Edge21 &, const Point21 &);
    bool contains(const Edge21 &, const Edge21 &);

    // Checks.

    bool spatial(const Edge21 &);

}

#endif