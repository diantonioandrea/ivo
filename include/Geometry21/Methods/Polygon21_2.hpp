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

    // (Meshing) polygon methods.

    Polygon21 reduce2(const Polygon21 &, const Line21 &, const Point21 &);
    
    std::vector<Point21> random2(const Polygon21 &, const Natural &);

    std::vector<Polygon21> voronoi2(const Polygon21 &, const std::vector<Point21> &);
    std::vector<Polygon21> voronoi2(const Polygon21 &, const Natural &);

    // (Meshing) diagram postprocessing.

    void lloyd2(const Polygon21 &, std::vector<Polygon21> &);
    void collapse2(const Polygon21 &, std::vector<Polygon21> &);

    // Polygon methods.

    std::array<Point21, 2> box2(const Polygon21 &);

    // Containment.

    bool contains2(const Polygon21 &, const Point21 &);

}

#endif