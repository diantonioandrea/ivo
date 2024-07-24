/**
 * @file Geometry21_Methods_Polygon21.cpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief Geometry/Methods/Polygon21.hpp implementation.
 * @date 2024-07-23
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <Ivo.hpp>

namespace ivo {

    // Checks.

    /**
     * @brief Spatial(polygon).
     * 
     * @param polygon Polygon.
     * @return true 
     * @return false 
     */
    bool spatial(const Polygon21 &polygon) {
        std::vector<Point21> points = polygon.points();

        for(Natural j = 1; j < points.size(); ++j)
            if(std::abs(points[j](2) - points[0](2)) > GEOMETRICAL_ZERO)
                return false;

        return true;
    }

}