/**
 * @file Line21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Line21 related methods.
 * @date 2024-07-22
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY21_METHODS_LINE21
#define GEOMETRY21_METHODS_LINE21

#include "../Polygon21.hpp"

namespace ivo {

    // Intersections.

    std::vector<Point21> intersections(const Line21 &, const Polygon21 &);
    std::optional<Point21> intersections(const Line21 &, const Edge21 &);
    std::optional<Point21> intersections(const Line21 &, const Line21 &);

    // Distances.

    Real distance(const Line21 &, const Point21 &);
    Real distance(const Line21 &, const Edge21 &);
    Real distance(const Line21 &, const Line21 &);

    // Containment.

    bool contains(const Line21 &, const Point21 &);
    bool contains(const Line21 &, const Edge21 &);

    // Checks.

    bool spatial(const Line21 &);

}

#endif