/**
 * @file Polygon21_2.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Polygon21 related methods. Space only.
 * @date 2024-07-23
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY21_METHODS_POLYGON21_2
#define GEOMETRY21_METHODS_POLYGON21_2

#include "../Includes.hpp"
#include "../Polygon21.hpp"

namespace ivo {

    // Polygon methods.

    std::array<Point21, 2> box(const Polygon21 &);

    // Containment.

    bool contains2(const Polygon21 &, const Point21 &);

}

#endif