/**
 * @file Polygon21.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Polygon21 related methods.
 * @date 2024-07-23
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef GEOMETRY21_METHODS_POLYGON21
#define GEOMETRY21_METHODS_POLYGON21

#include "../Includes.hpp"
#include "../Polygon21.hpp"

namespace ivo {

    // Polygon methods.

    Real area(const Polygon21 &);

    Point21 centre(const Polygon21 &);
    Point21 centroid(const Polygon21 &);

    // Checks.

    bool spatial(const Polygon21 &);

}

#endif